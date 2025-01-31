""" spider.py -- create a CSV of single-station observation metadata. """

import os
import sys
import glob
import pandas as pd
import tqdm
import numpy as np
from astropy.time import Time
from astropy.coordinates import EarthLocation
import h5py
import yaml
from dataclasses import dataclass
import pylab as plt
from datetime import datetime
from loguru import logger

# Setup logger
logger.remove(0)
logger.add(sys.stderr, format="<level>{level}</level> | {message}", level="INFO")

# Set SKA-Low Earth Location
eloc = EarthLocation(x=-2561290.83467119, y=5085918.51537833, z=-2864050.87177975, unit='m')

# Define a dataclass to store metadata in
@dataclass
class ObservationMetadata:
    """ Simple dataclass to store observation metadata. """
    obs_id: str
    station: str=''	
    mode: str=''
    sub_mode: str=''	
    intent: str=''	
    notes: str=''	
    observer: str=''	
    reference: str=''	
    utc_start: str=''	
    obs_duration: str='' 	
    loop_duration: str=''
    loop_rest_interval: str=''	
    lst_start: float=''
    start_channel: int=''	
    n_channel: int=''	
    start_frequency: float=''	
    obs_bandwidth: float=''	
    samples_per_frame: int=''	
    time_resolution: float=''	
    n_timesteps: int=''	
    tracking: str=''	
    source_name: str=''	
    right_ascension: str=''	
    declination: str=''	
    altitude: float=''	
    azimuth: float=''	
    obs_date: str=''	
    n_files: int=0	
    size_mb: float=0.0	
    qa: str=''

# Convert HDF5 name to a mode and submode
obs_types = {
    'stationbeam_integ': ['beamformer', 'power'],
    'channel_burst': ['channel-voltages', 'sweep'],
    'channel_integ': ['antenna-bandpass', ' '],
    'correlation_burst': ['correlator', 'sweep'], # NOTE: Could be fixed submode too
    'raw_burst': ['adc-samples', 'synchronous'],  # NOTE: Could be asynchronous submode too
}

###################
## Metadata readers
###################

def get_md_corr(filelist: list, obs_md: ObservationMetadata) -> ObservationMetadata:
    """ Get correlator metadata from filelist. """
    
    bl = [os.path.basename(f) for f in sorted(filelist)]
    dirname = os.path.dirname(filelist[0])
    
    chan_ids, tsteps = [], []
    modestr1, modestr2, chan_id, dt0, dt1, tstep = bl[0].replace('.hdf5', '').split('_')
    for fn in bl:
        _modestr1, _modestr2, chan_id, dt0, dt1, tstep = fn.replace('.hdf5', '').split('_')
        chan_ids.append(chan_id)
        tsteps.append(tstep)
        if modestr1 != _modestr1:
            logger.warning(f"{obs_md.obs_id} - Mixed HDF5 file types in {dirname}")

    # Get unique channels and sequence IDs (timesteps)
    chans = sorted([int(c) for c in set(chan_ids)])
    tsteps = sorted([int(t) for t in set(tsteps)])
    
    n_chan = len(chans)
    n_step = len(tsteps)
    
    obs_md.n_channel     = n_chan
    obs_md.n_timesteps   = n_step
    obs_md.start_channel = np.min(chans)
    obs_md.start_frequency = 0.78125 * obs_md.start_channel
    obs_md.obs_bandwidth   = 0.78125 * obs_md.n_channel 

    # Figure out if sweep mode (n_step == 1) or fixed mode (n_chan == 1)
    try:
        assert n_chan == 1 or n_step == 1
    except AssertionError:
        raise RuntimeError(f"Cannot determine mode! n_chan {n_chan} n_step {n_step}")

    obs_md.sub_mode = 'fixed' if n_chan == 1 else 'sweep'

    # Find first file, extract metadata
    first_file = f"{modestr1}_{modestr2}_{chans[0]}_{dt0}_{dt1}_{tsteps[0]}.hdf5"
    with h5py.File(os.path.join(dirname, first_file), mode='r') as h:
        t0 = Time(h['sample_timestamps']['data'][0, 0], format='unix')
        obs_md.time_resolution = np.round( h['root'].attrs['tsamp'], 11)
        obs_md.utc_start = t0.iso
        obs_md.n_timesteps *= h['sample_timestamps']['data'].shape[0]

    # Get last timestamp: in sweep mode, use last channel
    try:
        if obs_md.sub_mode == 'sweep':
            last_file = f"{modestr1}_{modestr2}_{chans[-1]}_{dt0}_{dt1}_{tsteps[0]}.hdf5"
            
            with h5py.File(os.path.join(dirname, last_file), mode='r') as h:
                t1 = Time(h['sample_timestamps']['data'][0, 0] + obs_md.time_resolution, format='unix')
                
        # Get last timestamp: in fixed mode, this will be last file in sequence
        else:
            if n_step > 1:
                last_file = f"{modestr1}_{modestr2}_{chans[0]}_{dt0}_{dt1}_{tsteps[-1]}.hdf5"
                
                with h5py.File(os.path.join(dirname, last_file), mode='r') as h:
                    t1 = Time(h['sample_timestamps']['data'][-1, -1] + obs_md.time_resolution, format='unix')
            else:
                with h5py.File(os.path.join(dirname, first_file), mode='r') as h:
                    rr = h['root'].attrs
                    t1 = Time(rr['ts_start'] + rr['tsamp'], format='unix')
        
        obs_md.obs_duration = np.round((t1 - t0).sec, 11)

    except FileNotFoundError:
        logger.warning(f"{obs_md.obs_id} - file not found: {last_file} (Mixed file types?)")
            
    return obs_md


def get_md_adc(filelist: list, obs_md: ObservationMetadata) -> ObservationMetadata:
    """ Get ADC sample metadata from filelist. """
    with h5py.File(filelist[0], mode='r') as h:
        
        # ADC dataset is always 32768 in size. Here we check if data > 4096 are zeros,
        # which indicates synchronous mode was used
        dcount = np.sum(h['raw_']['data'][:][4096:])
        obs_md.sub_mode = 'synchronous' if dcount == 0 else 'asynchronous'
        obs_md.time_resolution = np.round(1 / 800e6, 11)
        n_samp = 4096 if dcount == 0 else 32768
        obs_md.obs_duration = n_samp * obs_md.time_resolution
        obs_md.n_timesteps  = n_samp
        
    return obs_md

def get_md_power_beam(filelist: list, obs_md: ObservationMetadata) -> ObservationMetadata:
    """ Get station beamformer (power) metadata from filelist. """
    # This will sort files in time order
    # Files appear to follow stationbeam_integ_0_20240918_11242_0.hdf5 (X and SEQ zero)
    filelist = sorted(filelist)
    dirname = os.path.dirname(filelist[0])

    
    with h5py.File(filelist[0], mode='r') as h:
        obs_md.time_resolution = h['root'].attrs['tsamp']
        obs_md.n_timesteps  = h['root'].attrs['n_blocks']
        obs_md.n_channel    = h['root'].attrs['n_chans']
        
        t0 = h['root'].attrs['ts_start']
        t1 = h['root'].attrs['ts_end']
        
    if len(filelist) > 1:
        with h5py.File(filelist[-1], mode='r') as h:
            t1 = h['root'].attrs['ts_end']
    obs_md.obs_duration = np.round(t1 - t0, 11)
    obs_md.n_timesteps *= len(filelist)
        
    return obs_md


##############
## Main loop
##############

def run_spider(datapath: str, outdir: str='db'):
    """ Run the spider script to create a CSV file.

    Loops through execution block directories in datapath
    """
    logger.info(f"Spidering {datapath}...")
    obslist = sorted(glob.glob(f"{datapath}/eb-t*"))
    
    obs_table = []
    for obs in tqdm.tqdm(obslist):
        subdirlist = os.listdir(f"{obs}/ska-low-mccs")
        n_subdir = len(subdirlist)
        obs_id = os.path.basename(obs)
        
        if n_subdir >= 1:
            obs_md = ObservationMetadata(obs_id)
            if n_subdir >= 2:
                logger.warning(f"{obs_id} - multiple scan directories") 
            for subdir in subdirlist:
                obspath = f"{obs}/ska-low-mccs/{subdir}"
                h5list = sorted(glob.glob(f"{obspath}/*.hdf5"))
                data_size_MB = 0
                    
                if len(h5list) > 0:
                    # Compute data volume
                    for h5 in h5list:
                        data_size_MB += os.path.getsize(h5) / 1e6
                    obs_md.size_mb = np.round(data_size_MB, 2)
                    obs_md.n_files = len(h5list)
                    
                    # Get metadata from HDF5 files
                    with h5py.File(h5list[0], 'r') as fh:
                        try:
                            t = Time(fh['root'].attrs['ts_start'], format='unix', location=eloc)
                            obs_md.utc_start = t.iso
                            #dstr = t.strftime("%Y-%m-%d")
                            #ut = t.unix
                            obs_md.lst_start = np.round(t.sidereal_time('apparent').value, 3)
                            
                            if fh['root'].attrs['station_id'] != 0:
                                station_id = str(fh['root'].attrs['station_id'])
    
                                # As of Jan 2025 some stations have integer IDs. This fixes 'em
                                sid_map = {'345': 's8-1', '350': 's8-6', '352': 's9-2', '431': 's10-3'}
                                station_id = sid_map.get(station_id, station_id)
                                obs_md.station = station_id
                            else:
                                # Station ID is sometimes stored in description field
                                # e.g. "s10-3, sun SFT" or "s9-2 CAL TEST"
                                description = fh['observation_info'].attrs['description'].split(', ')
                                description = description[0].split(' ')
                                if description[0].lower().startswith('s'):
                                    obs_md.station = description[0].upper()
                                if len(description) > 1:
                                    obs_md.intent = description[1]
                        except KeyError:
                            logger.warning(f"Cannot read required keys from {obs_id} {h5list[0]}")
    
                    # Get metadata from YAML
                    yaml_path = f"{obspath}/obs_metadata.yaml"
                    if os.path.exists(yaml_path):
                        with open(yaml_path, 'r') as fh:
                            log = yaml.safe_load(fh)
                            obs_md.station         = log.get('station', '')
                            obs_md.intent          = log.get('intent', '')
                            obs_md.observer        = log.get('observer', '')
                            obs_md.reference       = log.get('reference', '')
                            obs_md.notes           = log.get('notes', '')
                            obs_md.qa              = log.get('qa', '')
                            obs_md.tracking        = log.get('tracking', '')
                            obs_md.source_name     = log.get('source_name', '')
                            obs_md.right_ascension = log.get('ra', '')
                            obs_md.declination     = log.get('dec', '')
                            obs_md.altitude        = log.get('alt', '')
                            obs_md.azimuth         = log.get('az', '')
                            
                    # Get mode and submode from HDF5 filename
                    obstype = "_".join(os.path.basename(h5list[0]).split('_')[:2])
                    mode, sub_mode = obs_types[obstype]
    
                    obs_md.mode = mode
                    obs_md.sub_mode = sub_mode
    
                    try:
                        if mode == 'correlator':
                            obs_md = get_md_corr(h5list, obs_md)
                        elif mode == 'adc-samples':
                            obs_md = get_md_adc(h5list, obs_md)
                        elif mode == 'beamformer' and sub_mode == 'power':
                            obs_md = get_md_power_beam(h5list, obs_md)
                        obs_table.append(obs_md)
                    except OSError:
                        logger.warning(f"OSError when opening: {obs_id} {h5list[0]}")

    
    df = pd.DataFrame(obs_table)

    # Convert floats into string (helps searching)
    df['lst_start']    = df['lst_start'].round(1).astype('str')
    df['obs_duration'] = df['obs_duration'].round(3).astype('str')

    # Display-friendly names
    col_names = {
        'lst_start': 'LST start (hr)',
        'obs_duration': 'Duration (s)',
        'obs_id': 'Observation ID',
        'mode': 'Mode',
        'station': 'Station ID',
        'sub_mode': 'Sub-mode',
        'utc_start': 'UTC Start',
        'qa': 'QA',
        'bandwidth': 'Bandwidth (MHz)',
        'observer': 'Observer',
    }
    df = df.rename(columns=col_names)

    filename = datetime.now().strftime(f"{outdir}/%Y-%m-%d.csv")
    logger.info(f"Saving to {filename}")
    df.to_csv(filename, index=False)
    # Save to latest.csv also
    df.to_csv(f"{outdir}/latest.csv", index=False)

if __name__ == "__main__":
    dpath = "/home/jovyan/daq-data"
    run_spider(dpath, outdir='db')
    