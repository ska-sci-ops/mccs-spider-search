""" spider.py -- create a CSV of single-station observation metadata. """

import os
import sys
import re
import glob
import pandas as pd
import tqdm
import numpy as np
from astropy.time import Time
from astropy.coordinates import EarthLocation
import h5py
import yaml
from dataclasses import dataclass
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
from loguru import logger

# Map of all station IDs
from spider.station_ids import sid_map

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
    'channel_cont': ['channel-voltages', 'fixed'],
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


######################
## Per-observation work
######################

def _process_obs(obs: str) -> list:
    """ Build ObservationMetadata record(s) for a single eb-* folder.

    Pulled out of the main loop so it can be dispatched to a thread pool --
    the work here is almost entirely I/O (HDF5 reads, YAML reads, stat
    calls), so running several observations concurrently overlaps that I/O
    instead of waiting on it one folder at a time.

    Returns a list (usually 0 or 1 entries, but preserves the original
    behaviour of appending one row per valid scan subdirectory in the rare
    case a folder has more than one).
    """
    obs_id = os.path.basename(obs)
    subdirlist = os.listdir(f"{obs}/ska-low-mccs")
    n_subdir = len(subdirlist)

    if n_subdir < 1:
        return []

    obs_md = ObservationMetadata(obs_id)
    if n_subdir >= 2:
        logger.warning(f"{obs_id} - multiple scan directories")

    results = []
    for subdir in subdirlist:
        obspath = f"{obs}/ska-low-mccs/{subdir}"

        # A single scandir() pass gets us both the file list and its size
        # (DirEntry.stat() is cached from the directory read), instead of
        # glob() for names followed by a separate os.path.getsize() stat
        # call per file.
        with os.scandir(obspath) as it:
            h5_entries = sorted(
                (e for e in it if e.name.endswith('.hdf5')),
                key=lambda e: e.name,
            )
        h5list = [e.path for e in h5_entries]
        data_size_MB = sum(e.stat().st_size for e in h5_entries) / 1e6

        if len(h5list) > 0:
            obs_md.size_mb = np.round(data_size_MB, 2)
            obs_md.n_files = len(h5list)

            # Get metadata from HDF5 files
            try:
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
                            # station_id in
                            # https://gitlab.com/ska-telescope/ska-telmodel-data/-/blob/main/tmdata/instrument/ska1_low/layout/data.json
                            station_id = sid_map.get(str(station_id), station_id)
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
                    except OSError:
                        logger.warning(f"OSError encountered on {obs_id} {h5list[0]} (reading metadata)")
                    except KeyError:
                        logger.warning(f"Cannot read required keys from {obs_id} {h5list[0]}")

                    try:
                        # Sometimes the source RA/DEC or name is in the description field.
                        description = fh['observation_info'].attrs['description']
                        paren_match = re.search(r"\(([^)]+)\)", description)
                        if paren_match:
                            paren_value = paren_match.group(1)
                            if paren_value == "NAMED":
                                target_match = re.search(r"pointing at ([A-Z0-9_-]+) \(NAMED\)", description)
                                target = target_match.group(1) if target_match else None
                                obs_md.source_name = target
                            else:
                                radec_match = re.search(r"\[([^\]]+)\]", description)
                                radec = radec_match.group(1) if radec_match else None
                                ra, dec = [float(x.strip()) for x in radec.split(",")]
                                if paren_value == "ALT_AZ":
                                    obs_md.altitude = ra
                                    obs_md.azimuth = dec
                                else:
                                    obs_md.right_ascension = ra
                                    obs_md.declination = dec
                                obs_md.source_name = f"{paren_value} [{ra:.3f}, {dec:.3f}]"
                    except OSError:
                        logger.warning(f"OSError encountered on {obs_id} {h5list[0]} (reading observation_info)")
                    except KeyError:
                        logger.warning(f"Cannot read required key 'description' from {obs_id} {h5list[0]}")

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
                try:
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
                        results.append(obs_md)
                    except OSError:
                        logger.warning(f"OSError when opening: {obs_id} {h5list[0]}")
                    except KeyError:
                        logger.warning(f"Cannot read required keys from {obs_id} {h5list[0]}")
                except KeyError:
                    logger.warning(f"KeyError for obstype: {obstype} {h5list[0]}")
            except PermissionError:
                logger.warning(f"Permission error reading {h5list[0]}")
            except OSError:
                logger.warning(f"OSError when opening: {obs_id} {h5list[0]}")

    return results


##############
## Main loop
##############

def run_spider(datapath: str, outdir: str='db', skip_existing: bool=True, max_workers: int=8):
    """ Run the spider script to create a CSV file.

    Loops through execution block directories in datapath.

    Args:
        datapath: Path to the directory containing eb-* observation folders.
        outdir: Output directory for CSV files.
        skip_existing: If True, read latest.csv and skip any folders whose
                       obs_id is already present, appending new rows to the
                       existing records. Defaults to True.
        max_workers: Number of observation folders to process concurrently.
                     The per-observation work is I/O-bound (HDF5/YAML reads,
                     stat calls), so a thread pool overlaps that I/O instead
                     of processing one folder at a time. Defaults to 8.
    """
    logger.info(f"Spidering {datapath}...")
    # Match both the old eb-t<date>-<seq> folders and the newer
    # short-alphanumeric eb-<id> folders (e.g. eb-18jjem39q2h9).
    obslist = sorted(glob.glob(f"{datapath}/eb-*"))

    # Display-friendly names. Defined here (not just later, near the rename
    # of the newly-spidered df) so it can also be used to normalize whatever
    # column names happen to already be in latest.csv when we load it below.
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

    # Load existing records and build a set of already-spidered obs_ids
    existing_df = None
    existing_ids: set = set()
    if skip_existing:
        latest_csv = os.path.join(outdir, "latest.csv")
        if os.path.exists(latest_csv):
            try:
                existing_df = pd.read_csv(latest_csv)
                # latest.csv may have been written by an older version of this
                # script (raw dataclass field names) or a newer one (display
                # names). Normalize to display names either way -- rename()
                # silently leaves any column not in col_names untouched, so
                # this is a no-op for a file that's already in the new format.
                existing_df = existing_df.rename(columns=col_names)
                id_col = 'Observation ID' if 'Observation ID' in existing_df.columns else 'obs_id'
                existing_ids = set(existing_df[id_col].dropna().astype(str))
                logger.info(f"Loaded {len(existing_ids)} existing records from {latest_csv}; skipping those folders.")
            except Exception as exc:
                logger.warning(f"Could not read {latest_csv}: {exc}. Spidering all folders.")

    # Filter out already-spidered folders up front (this is what actually
    # keeps a daily incremental run fast -- everything below only ever runs
    # against the handful of new folders, not the whole archive).
    todo = [obs for obs in obslist if os.path.basename(obs) not in existing_ids]
    n_skipped = len(obslist) - len(todo)
    if n_skipped:
        logger.debug(f"Skipping {n_skipped} already-spidered folders")

    obs_table = []
    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        futures = {pool.submit(_process_obs, obs): obs for obs in todo}
        for future in tqdm.tqdm(as_completed(futures), total=len(futures)):
            obs = futures[future]
            try:
                obs_mds = future.result()
            except Exception as exc:
                logger.warning(f"Unhandled error processing {os.path.basename(obs)}: {exc}")
                continue
            obs_table.extend(obs_mds)

    df = pd.DataFrame(obs_table)

    # Convert floats into string (helps searching). Guard against an empty
    # DataFrame (e.g. every folder was already in latest.csv, or none had
    # readable metadata) -- in that case pd.DataFrame([]) has no columns at
    # all, so df['lst_start'] would raise a KeyError.
    if not df.empty:
        df['lst_start']    = pd.to_numeric(df['lst_start'], errors='coerce').round(1).astype('str')
        df['obs_duration'] = pd.to_numeric(df['obs_duration'], errors='coerce').round(3).astype('str')

    # Rename to display-friendly names (col_names is defined earlier, where
    # it's also used to normalize existing_df loaded from latest.csv).
    df = df.rename(columns=col_names)

    # Merge newly spidered rows with any existing records from latest.csv
    if skip_existing and existing_df is not None and not df.empty:
        df = pd.concat([existing_df, df], ignore_index=True)
    elif skip_existing and existing_df is not None and df.empty:
        logger.info("No new observations found; latest.csv is already up to date.")
        df = existing_df

    if df is None or df.empty:
        logger.warning("No observations found (new or existing); nothing to save.")
        return

    filename = datetime.now().strftime(f"{outdir}/%Y-%m-%d.csv")
    logger.info(f"Saving to {filename}")
    df.to_csv(filename, index=False)
    # Save to latest.csv also
    df.to_csv(f"{outdir}/latest.csv", index=False)

if __name__ == "__main__":
    dpath = "/home/jovyan/daq-data"
    now_str = datetime.now().strftime("%Y-%m-%d")
    logger.add(f"db/loguru_{now_str}.log", format="<level>{level}</level> | {message}", level="INFO")
    run_spider(dpath, outdir='db')

    print("Updating database on acacia")
    os.system("/home/jovyan/shared/Danny/rclone/rclone copy db/ SKAO:/aa05/mccs-spider-search/db")
