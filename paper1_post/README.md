# **Fig. Generating Scripts Documentation**

## **Figs. 1–2: UNSEED/MSEED restarts**
* **Script Location:** `/paper1_post/plt_restart_mod.py`

### **Execution Instructions**
Run the script after configuring the following variables at the top of the file:

* **`DIRO`**: Output directory where the generated figures will be saved.
* **`RSDIR`**: Parent directory that contains the date-specific CESM restart folder.
* **`DTSTR`**: Datetime string of the restart. The NetCDF files containing the model restart fields are expected to be located in `{RSDIR}/{DTSTR}`.
* **`EVNTS`**: Path to the Apache Parquet file containing the un/seeding interventions data.
* **`CLAT`, `CLON`**: Central latitude and longitude coordinates of the un/seeding intervention you want to plot. These must be manually specified and can be obtained from the CSV or Parquet file containing the interventions list.
