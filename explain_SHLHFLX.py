import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

Lv = 2.501e6
Rv = 461.5

FILI = 'linevslat_h0a_diff/sznlzm.nc'
SZN = 'SON'

LATBND = (-50, 50)
YLIM = [(-5, 5), (-25, 25)]

def es(T):
   aterm = -6810.5245 / T
   bterm = -5.08984 * np.log(T)
   cterm = 55.2966
   return 100 * np.exp(aterm + bterm + cterm)

ds = xr.open_dataset(FILI).sel(latitudes=slice(*LATBND), season=SZN)

dLH_frac = ds['LHFLX'].sel(case=['UNSEED', 'SEED']) / ds['LHFLX'].sel(case='CTRL') 

#print(dLH_frac)
dU_frac = ds['U10'].sel(case=['UNSEED', 'SEED']) / ds['U10'].sel(case='CTRL')

print(es(273.15))
dTs = ds['TS'].sel(case=['UNSEED', 'SEED'])
qsat = es(ds['TS'].sel(case='CTRL')) / ds['PS'].sel(case='CTRL')

dqsat = Lv / Rv * qsat / ds['TS'].sel(case='CTRL')**2 * dTs
#plt.plot(ds['latitudes'], dTs.sel(case='SEED'))
#plt.show()
dqa = ds['QREFHT'].sel(case=['UNSEED', 'SEED'])
qdef = qsat - ds['QREFHT'].sel(case='CTRL')
#plt.plot(ds['latitudes'], qdef)
#plt.show()

dqsat_frac = dqsat / qdef
dqa_frac = dqa / qdef

plt.rc('font', size=16)
plt.rcParams['figure.figsize'] = (15, 4.5)
fig, axes = plt.subplots(1, 2)

orilh = ds['LHFLX'].sel(case='CTRL')
for ii, cs in enumerate(['UNSEED', 'SEED']):
   axes[ii].plot(ds['latitudes'], orilh * dLH_frac.sel(case=cs), color='black', linewidth=3, label='truth')
   axes[ii].plot(ds['latitudes'], orilh * dU_frac.sel(case=cs), label='wind speed')
   axes[ii].plot(ds['latitudes'], orilh * dqsat_frac.sel(case=cs), label='SST')
   axes[ii].plot(ds['latitudes'], orilh * -dqa_frac.sel(case=cs), label='qa')
   axes[ii].plot(ds['latitudes'], orilh * (dU_frac + dqsat_frac - dqa_frac).sel(case=cs), label='linear estimate', linewidth=3, color='gray')
   axes[ii].axhline(0, linestyle='dotted', color='black')
   [axes[ii].axvline(ll, color='gray', linewidth=0.5) for ll in range(-50, 51, 10)]
   axes[ii].legend(fontsize=12)

   axes[ii].set_title(['(a)', '(b)'][ii], loc='left')
   axes[ii].set_title('%s$-$CTRL' % cs)
   axes[ii].set_xlabel('Latitude')
   axes[ii].set_ylabel('LHFLX [W m$^{-2}$]')
   axes[ii].set_ylim(*YLIM[ii])

fig.tight_layout()
plt.show()
