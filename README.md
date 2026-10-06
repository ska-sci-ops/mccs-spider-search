# mccs-spider-search

* spider/ (`mccs-spider`): Runs on MCCS server, extracts HDF5 file metadata and creates database.
* spider/analyse.py (`mccs-analyse`): Builds the analysis workbook (sweep QA, coverage histograms, calibration sets) from `db/latest.csv`.
* app.py: Runs on scs server, simple Panel app to view the MCCS observation database.

### Usage

To update the database, run `mccs-spider` on the MCCS server. 

To do so: 
* login to the Juptyerhub at https://k8s.mccs.low.internal.skao.int/jupyterhub/
* open a new terminal session
* cd `/home/jovyan/shared/Danny/mccs-spider-search`
* run `mccs-spider` (after `uv pip install .`; add `--no-skip-existing` to re-spider everything)

This will spider the files and create the database (a CSV file). This file is uploaded to 
acacia via rclone:

```
rclone copy db/latest.csv SKAO:/aa05/mccs-spider-search/db/
```

To build the analysis workbook from the database:

```
mccs-analyse [db/latest.csv] [out.xlsx] [calibration.xlsx]
```

Sweeps are QA'd on file count: 385 (64-448), or the daytime ranges 321 (128-448), 298 (128-425) and 362 (64-425).

To start the web app, run `start.sh` on the SCS01 server:
* Login via `infra login`
* SSH to `ssh 10.150.0.200`
* cd to `cd /shared/dancpr/mccs-spider-search`
* Activate conda `conda activate aa05`
* Run the start script `./start.sh`

To update the DB on scs01, run the `update_db.sh` script.

