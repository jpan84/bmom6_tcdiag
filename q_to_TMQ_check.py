#Joshua Pan

CASES = ['250415_unseed', '250417_ctrl', '250416_seed1x1']
FILI = '/glade/campaign/univ/upsu0032/jpan_aquaptc/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.%s/atm/hist/*.h0a.*.nc'
GRIDFN = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'

import os
import numpy as np
import uxarray as ux
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import sznl_funcs

P0 = 1e5
g = 9.81
FLD = 'Q'

LNCLR = ['blue', 'orange']
YSCL = lambda lat: np.sin(np.deg2rad(lat))
YLAB = np.array([-90, -60, -45, -30, -15, 0, 15, 30, 45, 60, 90]).astype(np.int_)
YLOC = YSCL(YLAB)
ZLIM = [(-1.5, 1.5), (0, 80), (-16, 16)]


def main():
   #if not os.path.exists(OUTDIR):
   #   os.makedirs(OUTDIR)

   dss = [ux.open_mfdataset(GRIDFN, FILI % cs) for cs in CASES]
   ds = dss[0]

   print('\nObtaining dp3d...')
   aterm = (ds['hyai'].isel(ilev=slice(1, None)) - ds['hyai'].isel(ilev=slice(None, -1)).data) * P0
   bterm = (ds['hybi'].isel(ilev=slice(1, None)) - ds['hybi'].isel(ilev=slice(None, -1)).data) * ds['PS']
   dp3d = aterm + bterm
   dp3d = dp3d.rename(dict(ilev='lev')).assign_coords(lev=ds['lev'])
   print(dp3d)
   dss = [dd.assign(variables=dict(dp3d=dp3d, coslat=np.cos(np.deg2rad(ds['lat'])))) for dd in dss]

   colint = [(dd[FLD] * dd['dp3d']).sum(dim='lev') / g for dd in dss] 
   zm = [da.groupby('time.month').mean().zonal_mean(lat=(-90, 90, 1.5)) for da in colint]
   sznl = [sznl_funcs.stack_hemi_sznl(sznl_funcs.monthly2sznl(da)) for da in zm]

   plt.rc('font', size=16)
   plt.rcParams['figure.figsize'] = (30, 6)
   subplot_kw = dict(xlim=(-1, 1))
   fig, axes = plt.subplot_mosaic([['(a)', '(b)', '(c)']], layout='constrained', sharey=False, subplot_kw=subplot_kw)

   for ii, (lbl, ax) in enumerate(axes.items()):
      toplt = sznl[ii]
      if ii in {0, 2}:
         toplt = sznl[ii] - sznl[1]
      for jj, szn in enumerate(toplt.season):
         ax.plot(YSCL(toplt['latitudes']), toplt.sel(season=szn), color=LNCLR[jj], linestyle='dashdot')
         ax.axhline(y=0, linestyle='dotted', color='gray')
         ax.set_xticks(YLOC, YLAB)
         ax.set_ylim(*ZLIM[ii])
         ax.set_title(lbl, loc='left')

   plt.savefig('TMQ_from_q3d.png')
   plt.close()


if __name__ == '__main__':
   main()
