from netCDF4 import Dataset
import numpy as np
import os

src = (
    "b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.251229_seedmatch."
    "mom6.sfc.0016-08-17.nc.oribadtime"
)

dst = (
    "b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.251229_seedmatch."
    "mom6.sfc.0016-08-17.nc.repaired"
)

# Variables whose HDF5 payloads are known to be corrupted.
# time is recreated below.
SKIP = {"time", "tos", "tauy_bot"}

# Recovered time values for 0016-08-17.
TIME_VALUES = np.array(
    [492760800., 492782400., 492804000., 492825600.],
    dtype=np.float64,
)

if not os.path.exists(src):
    raise FileNotFoundError(src)

if os.path.exists(dst):
    raise FileExistsError(
        f"Refusing to overwrite existing file:\n{dst}"
    )

print(f"Source: {src}")
print(f"Output: {dst}")
print()

with Dataset(src, "r") as old, Dataset(dst, "w", format="NETCDF4") as new:

    # ================================================================
    # Dimensions
    # ================================================================

    print("Creating dimensions...")

    for name, dim in old.dimensions.items():
        new.createDimension(
            name,
            None if dim.isunlimited() else len(dim)
        )

        print(
            f"  {name}: "
            f"{'UNLIMITED' if dim.isunlimited() else len(dim)}"
        )

    # ================================================================
    # Global attributes
    # ================================================================

    print("\nCopying global attributes...")

    for attr in old.ncattrs():
        new.setncattr(attr, old.getncattr(attr))

    # ================================================================
    # Variables
    # ================================================================

    print("\nCopying variables...\n")

    for name, var in old.variables.items():

        # ------------------------------------------------------------
        # Skip known-corrupted variables.
        # ------------------------------------------------------------

        if name in SKIP:
            print(f"SKIPPING {name}")
            continue

        print(
            f"Copying {name}: "
            f"{var.dimensions} {var.dtype}"
        )

        # ------------------------------------------------------------
        # Get _FillValue before creating variable.
        # ------------------------------------------------------------

        fill_value = None

        if "_FillValue" in var.ncattrs():
            fill_value = var.getncattr("_FillValue")

        # ------------------------------------------------------------
        # Preserve compression/chunking where possible.
        # ------------------------------------------------------------

        kwargs = {}

        if fill_value is not None:
            kwargs["fill_value"] = fill_value

        try:
            chunking = var.chunking()

            if chunking != "contiguous":
                kwargs["chunksizes"] = chunking

        except Exception:
            pass

        try:
            filters = var.filters()

            if filters is not None:

                if filters.get("zlib", False):
                    kwargs["zlib"] = True
                    kwargs["complevel"] = filters.get(
                        "complevel", 1
                    )

                if filters.get("shuffle", False):
                    kwargs["shuffle"] = True

                if filters.get("fletcher32", False):
                    kwargs["fletcher32"] = True

        except Exception:
            pass

        # ------------------------------------------------------------
        # Create destination variable.
        # ------------------------------------------------------------

        out = new.createVariable(
            name,
            var.dtype,
            var.dimensions,
            **kwargs
        )

        # ------------------------------------------------------------
        # Copy variable attributes.
        #
        # _FillValue must NOT be copied here because it was already
        # supplied to createVariable().
        # ------------------------------------------------------------

        for attr in var.ncattrs():

            if attr == "_FillValue":
                continue

            out.setncattr(
                attr,
                var.getncattr(attr)
            )

        # ------------------------------------------------------------
        # Copy data.
        #
        # For variables with time as the first dimension, copy one
        # timestep at a time. This avoids loading huge variables into
        # memory.
        # ------------------------------------------------------------

        if (
            "time" in var.dimensions
            and var.ndim > 0
            and var.dimensions[0] == "time"
        ):

            for i in range(var.shape[0]):

                print(
                    f"    timestep {i + 1}/{var.shape[0]}"
                )

                index = (
                    (i,)
                    + (slice(None),) * (var.ndim - 1)
                )

                out[index] = var[index]

        else:

            print("    copying complete variable")

            out[:] = var[:]

    # ================================================================
    # Recreate time
    # ================================================================

    print("\nCreating repaired time variable...")

    time = new.createVariable(
        "time",
        "f8",
        ("time",),
        zlib=True,
        complevel=1,
        shuffle=True,
        chunksizes=(512,),
    )

    time[:] = TIME_VALUES

    time.setncattr(
        "units",
        "seconds since 0001-01-01 00:00:00"
    )
    time.setncattr(
        "long_name",
        "time"
    )
    time.setncattr(
        "axis",
        "T"
    )
    time.setncattr(
        "calendar",
        "noleap"
    )

    print("  values:", TIME_VALUES)

print()
print("================================================")
print("Finished.")
print("================================================")
print()
print(f"Repaired file:")
print(dst)
