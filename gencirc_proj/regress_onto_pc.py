import glob
import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

PCDIR = '/glade/derecho/scratch/jpan/jpan_tcfields/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl/hist_0012-0014_h1i/EHF_EMF_500_warmNH_0012-0014_eofs'
PCGL = '*pcs*'
ANOMS = '/glade/derecho/scratch/jpan/jpan_tcfields/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl/hist_0012-0014_h1i/EHF_EMF_500_warmNH_0012-0014_anom.nc'
VAR = 'EHF500'
LATNM = 'latitudes'
#MEAN = '/glade/derecho/scratch/jpan/jpan_tcfields/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl/hist_0012-0014_h1i/EHF_EMF_500_warmNH_0012-0014_timmean.nc'
TDEV = '/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.251010_ctrlbr/atm/0012_JJASON_onpres_0.25_zm_timdev.nc'
MEAN = '/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.251010_ctrlbr/atm/0012_JJASON_onpres_0.25_zm_tm.nc'

meds = xr.open_dataset(MEAN).squeeze()
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

#perform/test regression by reconstructing anom snapshots
ands = xr.open_dataset(ANOMS)
lrcoef = (pcnorm * ands[VAR]).sum(dim='time') / (pcnorm**2).sum(dim='time') #anom of physical quantity corresponding to +1std anom in PC
#print(lrcoef)

'''
tsel = 1888
plt.rc('font', size=16)
plt.plot(ands[LATNM], ands[VAR].isel(time=tsel), lw=2.5, c='black')
recon = [(pcnorm.isel(time=tsel, mode=slice(0, ii)) * lrcoef.isel(mode=slice(0, ii))).sum(dim='mode') for ii in range(1, 18, 4)]
#print(recon[1])
[plt.plot(ands[LATNM], rec, label=1 + ii * 4) for ii, rec in enumerate(recon)]
plt.title(ands.time.isel(time=tsel).data)
plt.xlabel('Lat [°]')
plt.ylabel('EHF500 anom [K m/s]')
plt.legend()
plt.savefig('recon_%d.png' % tsel, bbox_inches='tight')
#plt.show()
plt.close()

#replot dimensional EOF loading patterns
plt.plot(meds[LATNM], meds[VAR], c='black', lw=2.5)
[plt.axvline(ll, lw=0.5, c='gray') for ll in np.arange(-5, 35, 5)]
plt.ylabel('Mean EHF [K m/s]')
plt.xlabel('Lat [°]')
ax1 = plt.gca().twinx()
ax1.axhline(0, lw=0.5, c='gray')
[plt.plot(lrcoef[LATNM], lrcoef.isel(mode=ii), label=ii + 1) for ii in range(4)]
plt.ylabel('EOF loading patterns [K m/s]')
#plt.show()
plt.legend()
plt.savefig('EOF500_patterns_dimensional.png', bbox_inches='tight')
plt.close()
'''

var2d = 'SST'
LATNM = 'lat'
for lag in np.arange(-16, 20, 4):
   varreg = reg_on_pc(ofld[var2d].shift(time=lag * 4))
   plt.plot(meds[LATNM], meds[var2d], c='black', lw=2.5)
   [plt.axvline(ll, lw=0.5, c='gray') for ll in np.arange(-5, 35, 5)]
   plt.ylabel('Mean %s [%s]' % (var2d, meds[var2d].attrs['units']))
   plt.xlabel('Lat [°]')
   plt.ylim(290, 310)
   ax1 = plt.gca().twinx()
   ax1.axhline(0, lw=0.5, c='gray')
   [plt.plot(varreg[LATNM], varreg.isel(mode=ii), lw=2, label=ii + 1) for ii in range(2)]
   plt.ylabel('%s regressed on PC [%s]\nlag = %d d' % (var2d, meds[var2d].attrs['units'], lag))
   plt.xlim(-5, 30)
   ax1.set_ylim(-.04, .08)
   #plt.show()
   plt.legend()
   plt.savefig('%s_reg_d%d.png' % (var2d, lag), bbox_inches='tight')
   plt.close()

var3d = 'U'
LATNM = 'lat'
PNM = 'plev'
for lag in np.arange(-8, 10, 2):
   varreg = reg_on_pc(ofld[var3d].shift(time=lag * 4))
   print(varreg)
   for md in varreg['mode'][:2]:
      plt.contourf(varreg[LATNM], varreg[PNM], varreg.sel(mode=md), levels=np.arange(-1, 1.01, 0.05), cmap='seismic')
      plt.ylim(1e5, 1e4)
      plt.yscale('log')
      plt.colorbar()
      #plt.show()
      plt.title('mode %d, lag %d d' % (md + 1, lag))
      plt.savefig('%s_reg_mode%d_d%d.png' % (var3d, md + 1, lag), bbox_inches='tight')
      plt.close()
   continue


   plt.plot(meds[LATNM], meds[var2d], c='black', lw=2.5)
   [plt.axvline(ll, lw=0.5, c='gray') for ll in np.arange(-5, 35, 5)]
   plt.ylabel('Mean %s [%s]' % (var2d, meds[var2d].attrs['units']))
   plt.xlabel('Lat [°]')
   plt.ylim(290, 310)
   ax1 = plt.gca().twinx()
   ax1.axhline(0, lw=0.5, c='gray')
   [plt.plot(varreg[LATNM], varreg.isel(mode=ii), lw=2, label=ii + 1) for ii in range(2)]
   plt.ylabel('%s regressed on PC [%s]\nlag = %d d' % (var2d, meds[var2d].attrs['units'], lag))
   plt.xlim(-5, 30)
   ax1.set_ylim(-.04, .08)
   #plt.show()
   plt.legend()
   plt.savefig('%s_reg_d%d.png' % (var2d, lag), bbox_inches='tight')
   plt.close()
