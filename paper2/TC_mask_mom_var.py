import glob
import os
import sys
import dask
import xarray as xr
import numpy as np

# Target variable (defaults to 'hflso' if not specified)
VAR_NAME = sys.argv[1] if len(sys.argv) > 1 else "hflso"


def get_zonal_dim(da):
    """Detect the correct longitude/zonal dimension for the dataset."""
    possible_lon_dims = ["xh", "lon", "nlon", "i", "x"]
    for dim in possible_lon_dims:
        if dim in da.dims:
            return dim
    raise KeyError(
        f"Could not find a recognized zonal/longitude dimension in {list(da.dims)}. "
        f"Expected one of {possible_lon_dims}."
    )


def process_experiment(dr_path, ocn_base, var_name):
    dr = os.path.basename(dr_path)
    mask_dir = os.path.join(dr_path, "nff_4mps_sx0.66av1")
    ocn_dir = os.path.join(ocn_base, dr, "ocn", "hist")

    if not (os.path.isdir(mask_dir) and os.path.isdir(ocn_dir)):
        print(f"Skipping {dr}: Missing required directories.")
        return

    print(f"Processing experiment: {dr} for variable: {var_name}")

    # 1. Open mask dataset lazily
    mask_files = sorted(glob.glob(f"{mask_dir}/*000[6-7]*.nff_4mps"))
    ds_mask = xr.open_mfdataset(
        mask_files, chunks={"time": 30}, combine="nested", concat_dim="time"
    )

    # 2. Open ocean files lazily (selecting only target variable)
    ocn_files = sorted(glob.glob(f"{ocn_dir}/*.sfc.000[6-7]*.nc"))
    ds_ocn = xr.open_mfdataset(
        ocn_files,
        chunks={"time": 30},
        combine="nested",
        concat_dim="time",
        data_vars=[var_name],
    )

    # 3. Direct in-memory time alignment
    print(ds_mask['time'])
    print(ds_ocn['hflso']['time'])
    shared_time = np.intersect1d(ds_mask.time, ds_ocn.time)
    var_subset = ds_ocn[var_name].sel(time=shared_time)
    ds_mask = ds_mask.sel(time=shared_time)

    # 4. Element-wise multiplication with mask
    # Note: Xarray handles broadcasting automatically even if mask uses (lat, lon)
    # and ocean variable uses (yh, xh) as long as dimensions align or coordinates match.
    # If dimensions are named differently between mask and var, rename mask dimensions:
    mask_da = ds_mask["TC_R4"]
    if "yh" in var_subset.dims and "lat" in mask_da.dims:
        mask_da = mask_da.rename({"lat": "yh", "lon": "xh"})

    masked_var = var_subset * mask_da

    # 5. Dynamically detect zonal dimension ('xh', 'lon', 'nlon', etc.)
    zonal_dim = get_zonal_dim(var_subset)
    print(f"  Averaging across zonal dimension: '{zonal_dim}'")

    # Compute zonal means
    ocn_zonal = var_subset.mean(dim=zonal_dim, skipna=True)
    masked_zonal = masked_var.mean(dim=zonal_dim, skipna=True)

    # 6. Dynamically named output files
    out_masked = os.path.join(dr_path, f"{var_name}_masked_4mps.nc")
    out_ocn_zonal = os.path.join(dr_path, f"{var_name}_subset_zonal_mean.nc")
    out_masked_zonal = os.path.join(dr_path, f"{var_name}_masked_4mps_zonal_mean.nc")

    # Save outputs to disk in parallel
    print(f"  Writing output files for {var_name}...")
    xr.save_mfdataset(
        [
            masked_var.to_dataset(name=f"{var_name}_masked"),
            ocn_zonal.to_dataset(name=var_name),
            masked_zonal.to_dataset(name=f"{var_name}_masked"),
        ],
        [out_masked, out_ocn_zonal, out_masked_zonal],
    )
    print(f"  Finished {dr}")


if __name__ == "__main__":
    dask.config.set(scheduler="threads")

    MASK_BASE = "/glade/derecho/scratch/jpan/archive/nff_output"
    OCN_BASE = "/glade/campaign/univ/upsu0032/jpan_aquaptc"

    for dr_path in sorted(glob.glob(f"{MASK_BASE}/b.e23.*0702*")):
        if os.path.isdir(dr_path):
            process_experiment(dr_path, OCN_BASE, VAR_NAME)
