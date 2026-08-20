from netCDF4 import Dataset

f = "b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.251229_seedmatch.mom6.sfc.0016-08-17.nc.oribadtime"

with Dataset(f) as nc:
    for name, var in nc.variables.items():
        try:
            # Read one time slice
            if "time" in var.dimensions:
                i = var.dimensions.index("time")
                key = [slice(None)] * var.ndim
                key[i] = 0
                var[tuple(key)]
            else:
                var[:]

            print("OK  ", name)

        except Exception as e:
            print("BAD ", name, repr(e))

from netCDF4 import Dataset

f = "b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.251229_seedmatch.mom6.sfc.0016-08-17.nc.oribadtime"

with Dataset(f) as nc:
    for name in ["tauy_bot", "tos", "time"]:
        v = nc.variables[name]

        print("\n", name, v.shape, v.dtype)

        for i in range(v.shape[0]):
            try:
                x = v[i]
                print("  timestep", i, "OK")
            except Exception as e:
                print("  timestep", i, "BAD:", repr(e))

from netCDF4 import Dataset

f = "b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.251229_seedmatch.mom6.sfc.0016-08-16.nc"

with Dataset(f) as nc:
    for name in ["time", "tauy_bot", "tos", "taux_bot", "tauvo"]:
        v = nc.variables[name]
        print("\n", name)
        print(" shape:", v.shape)
        print(" dtype:", v.dtype)
        print(" dims :", v.dimensions)
        print(" attrs:", {a: v.getncattr(a) for a in v.ncattrs()})
        try:
            print(" chunking:", v.chunking())
            print(" compression:", v.filters())
        except Exception as e:
            print(" layout unavailable:", e)
