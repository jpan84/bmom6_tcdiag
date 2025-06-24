#learn to compute conservative global integrals with hybrid-level data
import uxarray as ux

TESTFI = '/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl/atm/hist/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl.cam.h0a.0014-04.nc'
camgrid = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'
P0 = 1e5

ds = ux.open_dataset(camgrid, TESTFI)
print(ds['hyai'])
print(ds['hybi'])
print(ds['hyam'])
print(ds['PS'])
#print(ds['P0'])

aterm = (ds['hyai'].isel(ilev=slice(1, None)) - ds['hyai'].isel(ilev=slice(None, -1)).data) * P0
print('\naterm')
print(aterm)

bterm = (ds['hybi'].isel(ilev=slice(1, None)) - ds['hybi'].isel(ilev=slice(None, -1)).data) * ds['PS']
print('\nbterm')
print(bterm)

dp3d = aterm + bterm
dp3d = dp3d.rename(dict(ilev='lev')).assign_coords(lev=ds['lev'])
print('\ndp3d')
print(dp3d)
