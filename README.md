# mccs-spider-search

* spider.py: Runs on MCCS server, extracts HDF5 file metadata and creates database.
* app.py: Runs on scs server, simple Panel app to view the MCCS observation database.

### Usage

To update the database, run `spider.py` on the MCCS server. 

To do so: 
* login to the Juptyerhub at https://k8s.mccs.low.internal.skao.int/jupyterhub/
* open a new terminal session
* cd `/home/jovyan/shared/Danny/mccs-spider-search`
* run `python spider.py`

This will spider the files and create the database (a CSV file). This file is uploaded to 
acacia via rclone:

```
rclone copy db/latest.csv SKAO:/aa05/mccs-spider-search/db/
```

To start the web app, run `start.sh` on the SCS01 server:
* Login via `infra login`
* SSH to `ssh 10.150.0.200`
* cd to `cd /shared/dancpr/mccs-spider-search`
* Activate conda `conda activate aa05`
* Run the start script `./start.sh`

To update the DB on scs01, run the `update_db.sh` script.

