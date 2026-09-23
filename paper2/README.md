# **Figure and Table Generating Scripts Documentation**
For now, I’ll treat my Github repo located at ~jpan/aquaptc/bmom6\_tcdiag/paper2 (on NCAR glade) as the root for all paths below. Make sure to set the correct paths within /../paper1\_post/paths.py: CAMGR should point to the SCRIP file for the atmospheric mesh, and ARCHRT should contain the directories known as $DOUT\_S\_ROOT in CESM parlance.

## **Shared .py variables and routines**
See paper1 repo for details on the following: `paths.py`, `consts.py`, `sznl_funcs.py`

## **Preprocessing tools**
See paper1 repo for the following categories of tools:

* PP1: Preprocessing TC trajectories in TempestExtremes
* PP1A: Preprocessing TC trajectory statistics
* PP2: Preprocessing coordinate aggregations/interpolations

PP2A: new for paper 2 because some high-freqency variables for masking analyses (e.g., surface latent heat flux hflso) are only available from ocean
* `/remap_TC_masks_mom.sh` maps TC r4 (or other binary masks) from the CAM grid to the MOM grid using an existing map file. The map file should by source-to-destination (stod) in order to maintain binarity.
* `/TC_mask_mom_var.py` applies TC r4 masks (that have already been interpolated to the ocean grid) to selected ocean vars and takes the zonal mean of the masked and unmasked fields

## **Figs. 1, 3: precip and evap absolute and difference line plots, contribution fractions**
**Script Location:** `/plt_TC_masked_frac.py`

### **Preprocessing requirements**
* **PP1**: run `par-track_driver.py` and `trajSN_to_df.py` to obtain TC trajectories. Run `nff\_driver.py` with `PRECT` as a filtvar to generate TC masks and masked precip.
* **PP2**: run `zonmean_driver.py` with `VARS = 'PRECT'` and `TAPE = 'atm/hist/*.h1i.*.nc'`. Run the script again with `VARS = 'TC_R4,PRECT'` and `TAPE = 'atm/nff_4mps/*.h1i.*.nff_4mps'`
* **PP2A**: run both scripts

### **Execution Instructions**
Run the script after pointing `TOTP` to the zonal-mean precip, `TCSP` to the zonal-mean TC masks and masked precip, `TOTE` to the zonal-mean ocean evaporation, and `TCSE` to the zonal-mean masked evaporation.

## **Figs. 2, 4: 850 hPa eddy moisture flux absolute and difference line plots, contribution fractions**
**Script Location:** `/vpqp_decomp.py`

### **Preprocessing requirements**
* **PP1**: run `par-track_driver.py` and `trajSN_to_df.py` to obtain TC trajectories. Run `nff_driver.py` with `V850`,`Q850` as filtvars. Then run `nff_driver.py` again with invert=True so that the inverted mask is applied to those fields.
* **PP2**: run `zonmean_driver.py` on `VARS = 'V850,Q850,V850.Q850'` for 1) the full fields, 2) the TC-masked fields, and 3) the fields that have had the inverted masks applied.

### **Execution Instructions**
Run the script after pointing totfil, tcsfil, and bkgfil to the respective zonal-mean files generated in preprocessing.

