import glob
import os
import xarray as xr
import matplotlib.pyplot as plt

PCDIR = '/glade/derecho/scratch/jpan/jpan_tcfields/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl/hist_0012-0014_h1i/EHF_EMF_500_warmNH_0012-0014_eofs'
PCGL = '*pcs*'
ANOMS = '/glade/derecho/scratch/jpan/jpan_tcfields/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl/hist_0012-0014_h1i/EHF_EMF_500_warmNH_0012-0014_anom.nc'
VAR = 'EHF500'
LATNM = 'latitudes'

globs = sorted(glob.glob(os.path.join(PCDIR, PCGL)))
dss = [xr.open_dataset(gl) for gl in globs]
pcds = xr.concat(dss, dim='mode').squeeze()
#print(pcds)

#normalize PCs to stdev=1
stdev = pcds[VAR].std(dim='time')
pcnorm = pcds[VAR] / stdev
#print(pcnorm)

#perform/test regression by reconstructing anom snapshots
ands = xr.open_dataset(ANOMS)
lrcoef = (pcnorm * ands[VAR]).sum(dim='time') / (pcnorm**2).sum(dim='time') #anom of physical quantity corresponding to +1std anom in PC
#print(lrcoef)

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
plt.show()
