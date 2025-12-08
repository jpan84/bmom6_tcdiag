import glob
import os
from datetime import timedelta as tdel
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as colors

VAR = 'eehf_ix_std'
IXFIL = './0012-0013_JJASON_EEHF_ix_std.nc'
EVFIL = './ehf_events_sep3_sig1.0.nc'
LATNM, PNM = 'lat', 'plev'
TDEV = '/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.251010_ctrlbr/atm/0012-0013_JJASON_onpres_1deg_EPF_anom.nc'
laghrs = np.arange(-336, 337, 48.)
laghrs = np.arange(-96, 97, 6.)
lagtdels = [tdel(hours=hh) for hh in laghrs]

ofld = xr.open_dataset(TDEV).squeeze() #other cesm hist fields to regress (preprocessed into zonal mean temporal anom)
ixds = xr.open_dataset(IXFIL)
evds = xr.open_dataset(EVFIL)

#print(evds)
#print(evds['peaks'])

#globs = sorted(glob.glob(os.path.join(PCDIR, PCGL)))
#dss = [xr.open_dataset(gl) for gl in globs]
#pcds = xr.concat(dss, dim='mode').squeeze()
#print(pcds)

reg_on_pc = lambda vec: (ixds[VAR] * vec).sum(dim='time') / (ixds[VAR]**2).sum(dim='time') #anom of physical quantity corresponding to +1std anom in index

print(ofld['EPy_EMF'])
print(ixds[VAR])
EPy = ofld['EPy_EMF'] + ofld['EPy_adv']
EPz = ofld['EPz_EHF'] + ofld['EPz_adv']
EPd = ofld['EPy_EMF_d'] + ofld['EPy_adv_d'] + ofld['EPz_EHF_d'] + ofld['EPz_adv_d']

bbox = (ofld['plev'] <= 850) & (ofld['plev'] >= 300) & (ofld['lat'] >= 0) & (ofld['lat'] <= 35)
print(bbox.sum() / bbox.size)

plt.rcParams['figure.figsize'] = (10, 6)
contourkwargs = {'colors': 'black', 'levels': .02 * 2.**np.arange(-1, 7, 1)}
contourkwargs['levels'] = np.concatenate((-contourkwargs['levels'][::-1], contourkwargs['levels']))
contourfkwargs = {'cmap': 'RdYlBu_r', 'levels': contourkwargs['levels'], 'norm': colors.SymLogNorm(0.01), 'extend': 'both'}
islc = dict(lat=slice(None, None, 4), plev=slice(None, -2))
for ii, mlag in enumerate(lagtdels):
   #yreg = reg_on_pc(EPy.shift(time=lag * 4)) #TODO: use datetime here to account for gap in JJASON data
   #zreg = reg_on_pc(EPz.shift(time=lag * 4))
   #divplt = reg_on_pc(EPd.isel(**islc).shift(time=lag * 4)) * 86400
   yreg = reg_on_pc(EPy.assign_coords(time=EPy['time'] + mlag))
   zreg = reg_on_pc(EPz.assign_coords(time=EPz['time'] + mlag))
   divplt = reg_on_pc(EPd.isel(**islc).assign_coords(time=EPd['time'] + mlag)) * 86400
   plty, pltz = yreg.isel(**islc) / 1e2, zreg.isel(**islc)
   pltlat, pltp = yreg[LATNM].isel(lat=islc[LATNM]), yreg[PNM].isel(plev=islc[PNM])

   pltlat = np.arange(-89.5, 89.6, 1.)[islc['lat']]

   #plt.contourf(varreg[LATNM], varreg[PNM], varreg.sel(mode=md), levels=np.arange(-1, 1.01, 0.05), cmap='seismic')
   csf = plt.contourf(pltlat, pltp, divplt, **contourfkwargs)
   qv = plt.quiver(pltlat, pltp, plty, pltz, scale=1e6, pivot='mid')
   plt.contour(pltlat, pltp, bbox.isel(**islc), levels=[0.5], zorder=99)
   plt.ylim(8.5e2, 5e1)
   plt.yscale('log')
   plt.gca().set_yticks([50, 70, 100, 200, 300, 500, 700, 850], [50, 70, 100, 200, 300, 500, 700, 850])
   plt.xlim(-10, 60)
   plt.colorbar(csf)
   #plt.show()
   plt.title('lag %dh' % -laghrs[ii])
   plt.savefig('%s_reg_ix_h%d.png' % ('EPF', -laghrs[ii]), bbox_inches='tight')
   plt.close()

   #print(evds['peaks'].data)
   #print(evds['peaks'].data + mlag)
   compo = []
   comptimes = np.intersect1d((evds['troughs'] + mlag), ofld['time'], assume_unique=True)
   compo.append(EPy.sel(time=comptimes).isel(**islc).mean(dim='time'))
   compo.append(EPz.sel(time=comptimes).isel(**islc).mean(dim='time'))
   compo.append(EPd.sel(time=comptimes).isel(**islc).mean(dim='time'))
   #snaps = [[], [], []]
   #for ev in (evds['peaks'].data + mlag):
   #   if ev in ofld['time']:
   #      snaps[0].append(EPy.sel(time=ev).isel(**islc))
   #      snaps[1].append(EPz.sel(time=ev).isel(**islc))
   #      snaps[2].append(EPd.sel(time=ev).isel(**islc))
   #compo = [sum(ll) / len(ll) for ll in snaps]
   divplt, plty, pltz = compo[2] * 86400, compo[0] / 1e2, compo[1]

   #plt.contourf(varreg[LATNM], varreg[PNM], varreg.sel(mode=md), levels=np.arange(-1, 1.01, 0.05), cmap='seismic')
   csf = plt.contourf(pltlat, pltp, divplt, **contourfkwargs)
   qv = plt.quiver(pltlat, pltp, plty, pltz, scale=1e6, pivot='mid')
   plt.contour(pltlat, pltp, bbox.isel(**islc), levels=[0.5], zorder=99)
   plt.ylim(8.5e2, 5e1)
   plt.yscale('log')
   plt.gca().set_yticks([50, 70, 100, 200, 300, 500, 700, 850], [50, 70, 100, 200, 300, 500, 700, 850])
   plt.xlim(-10, 60)
   plt.colorbar(csf)
   #plt.show()
   plt.title('lag %dh' % laghrs[ii])
   plt.savefig('PEHF_%s_cmp_ix_h%d.png' % ('EPF', laghrs[ii]), bbox_inches='tight')
   plt.close()
