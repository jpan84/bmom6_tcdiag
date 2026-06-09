# **Fig. Generating Scripts Documentation**
For now, I’ll treat my Github repo located at ~jpan/aquaptc/bmom6\_tcdiag/paper1\_post (on NCAR glade) as the root for all paths below. Make sure to set the correct paths within /paths.py: CAMGR should point to the SCRIP file for the atmospheric mesh, and ARCHRT should contain the directories known as $DOUT\_S\_ROOT in CESM parlance.

## **PP1: Preprocessing TC trajectories in TempestExtremes**
**Script Locations:** `/TC_preprocess/`

* **`par-track_driver.py`**: runs `par-track-aqua.derecho_worker.sh`, which generates TempestExtremes TC trajectories using DetectNodes and StitchNodes for each experiment using .h1i. files
* **`seedlog_merge.py`**: scrapes log files of un/seeding interventions to create tables of the vortex params for each intervention
* **`nff_driver.py`**: runs `nff_gen_TC_masks.sh`, which generates (using the .h1i. files and the output of `par-track_driver.py`) binary masks covering the region surrounding each node where wind speed exceeds X m/s

## **Figs. 1–2: UNSEED/MSEED restarts**
**Script Location:** `/plt_restart_mod.py`

### **Execution Instructions**
Run the script after configuring the following variables at the top of the file:

* **`DIRO`**: Output directory where the generated figures will be saved.
* **`RSDIR`**: Parent directory that contains the date-specific CESM restart folder.
* **`DTSTR`**: Datetime string of the restart. The NetCDF files containing the model restart fields are expected to be located in `{RSDIR}/{DTSTR}`.
* **`EVNTS`**: Path to the Parquet file containing the un/seeding interventions data.
* **`CLAT`, `CLON`**: Central latitude and longitude coordinates of the un/seeding intervention you want to plot. These must be manually specified and can be obtained from the CSV or Parquet file containing the interventions list.
