
from netCDF4 import Dataset

f = "b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.251229_seedmatch.mom6.sfc.0016-08-17.nc.oribadtime"

with Dataset(f) as nc:
    t = nc.variables["time"]

    print(t)
    print("dtype:", t.dtype)
    print("dimensions:", t.dimensions)
    print("shape:", t.shape)
    print("attributes:", {a: t.getncattr(a) for a in t.ncattrs()})

    try:
        print("values:", t[:])
    except Exception as e:
        print("READ FAILED:", repr(e))

from netCDF4 import Dataset

files = [
    "b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.251229_seedmatch.mom6.sfc.0016-08-16.nc",
    "b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.251229_seedmatch.mom6.sfc.0016-08-17.nc.oribadtime",
    "b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.251229_seedmatch.mom6.sfc.0016-08-18.nc",
]

for f in files:
    print("\n", f)
    with Dataset(f) as nc:
        t = nc.variables["time"]
        print("shape:", t.shape)
        print("units:", t.units)
        print("calendar:", t.calendar)
        try:
            print("values:", t[:])
        except Exception as e:
            print("READ FAILED:", e)
