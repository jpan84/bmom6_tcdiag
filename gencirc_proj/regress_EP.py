import glob
import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as colors

VAR = 'EHF500'
PCDIR = '/glade/derecho/scratch/jpan/jpan_tcfields/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl/hist_0012-0014_h1i/EHF_EMF_500_warmNH_0012-0014_eofs'
PCGL = '*pcs*'
LATNM, PNM = 'lat', 'plev'
TDEV = '/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.251010_ctrlbr/atm/0012_JJASON_onpres_0.25_EPF_anom.nc'

ofld = xr.open_dataset(TDEV).squeeze() #other cesm hist fields to regress (preprocessed into zonal mean temporal anom)
globs = sorted(glob.glob(os.path.join(PCDIR, PCGL)))
dss = [xr.open_dataset(gl) for gl in globs]
pcds = xr.concat(dss, dim='mode').squeeze()
#print(pcds)

#normalize PCs to stdev=1
stdev = pcds[VAR].std(dim='time')
pcnorm = pcds[VAR] / stdev
reg_on_pc = lambda vec: (pcnorm * vec).sum(dim='time') / (pcnorm**2).sum(dim='time') #anom of physical quantity corresponding to +1std anom in PC
#print(pcnorm)

EPy = ofld['EPy_EMF'] + ofld['EPy_adv']
EPz = ofld['EPz_EHF'] + ofld['EPz_adv']
EPd = ofld['EPy_EMF_d'] + ofld['EPy_adv_d'] + ofld['EPz_EHF_d'] + ofld['EPz_adv_d']

plt.rcParams['figure.figsize'] = (10, 6)
contourkwargs = {'colors': 'black', 'levels': .02 * 2.**np.arange(-1, 7, 1)}
contourkwargs['levels'] = np.concatenate((-contourkwargs['levels'][::-1], contourkwargs['levels']))
contourfkwargs = {'cmap': 'RdYlBu_r', 'levels': contourkwargs['levels'], 'norm': colors.SymLogNorm(0.01), 'extend': 'both'}
islc = dict(lat=slice(None, None, 16), plev=slice(None, -2))
for lag in np.arange(-8, 10, 2):
   yreg = reg_on_pc(EPy.shift(time=lag * 4))
   zreg = reg_on_pc(EPz.shift(time=lag * 4))
   divplt = reg_on_pc(EPd.isel(**islc).shift(time=lag * 4)) * 86400
   plty, pltz = yreg.isel(**islc) / 1e2, zreg.isel(**islc)
   pltlat, pltp = yreg[LATNM].isel(lat=islc[LATNM]), yreg[PNM].isel(plev=islc[PNM])
   for md in yreg['mode'][:2]:
      #plt.contourf(varreg[LATNM], varreg[PNM], varreg.sel(mode=md), levels=np.arange(-1, 1.01, 0.05), cmap='seismic')
      csf = plt.contourf(pltlat, pltp, divplt.sel(mode=md), **contourfkwargs)
      qv = plt.quiver(pltlat, pltp, plty.sel(mode=md), pltz.sel(mode=md), scale=5e5, pivot='mid')
      plt.ylim(1e3, 1e2)
      plt.yscale('log')
      plt.colorbar(csf)
      #plt.show()
      plt.title('mode %d, lag %d d' % (md + 1, lag))
      plt.savefig('%s_reg_mode%d_d%d.png' % ('EPF', md + 1, lag), bbox_inches='tight')
      plt.close()
