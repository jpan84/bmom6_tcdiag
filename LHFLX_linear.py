import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

Lv = 2.501e6
Rv = 461.5
YSCL = lambda deglat: np.sin(np.deg2rad(deglat))
MT = YSCL(np.array([10, 20]))
STZ = YSCL(np.array([20, 30]))

FILI = 'linevslat_h0a_diff/sznlzm.nc'
SZN = 'SON'

LATBND = (-50, 50)
YLIM = [(-6, 6), (-30, 30)]

def es(T):
   aterm = -6810.5245 / T
   bterm = -5.08984 * np.log(T)
   cterm = 55.2966
   return 100 * np.exp(aterm + bterm + cterm)

ds = xr.open_dataset(FILI).sel(latitudes=slice(*LATBND), season=SZN)
sinlat = YSCL(ds['latitudes'])

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
   axes[ii].plot(sinlat, orilh * dLH_frac.sel(case=cs), color='black', linewidth=3, label='truth')
   axes[ii].plot(sinlat, orilh * (dU_frac + dqsat_frac - dqa_frac).sel(case=cs), label='linear estimate', linewidth=3, color='gray')
   axes[ii].plot(sinlat, orilh * dU_frac.sel(case=cs), label='wind speed')
   axes[ii].plot(sinlat, orilh * dqsat_frac.sel(case=cs), label='SST')
   axes[ii].plot(sinlat, orilh * -dqa_frac.sel(case=cs), label='qa')

   axes[ii].axhline(0, c='gray', lw=0.5)
   axes[ii].axvspan(*MT, fc='purple', alpha=.08)
   axes[ii].axvspan(*STZ, fc='yellow', alpha=.12)
   [axes[ii].axvline(YSCL(ll), color='gray', linewidth=0.5) for ll in range(-50, 51, 10)]
   axes[ii].set_xticks(YSCL(np.arange(-50, 51, 10)), np.arange(-50, 51, 10))
   axes[ii].legend(fontsize=12)
   axes[ii].tick_params(right=True)

   axes[ii].set_title(['(a)', '(b)'][ii], loc='left')
   axes[ii].set_title('%s$-$CTRL' % cs)
   axes[ii].set_xlabel('Latitude')
   axes[ii].set_ylabel('LHFLX [W m$^{-2}$]')
   axes[ii].set_ylim(*YLIM[ii])

fig.tight_layout()
#plt.show()
plt.savefig('./linevslat_h1i_0012-0014/LH_decomp_%s_h0a.svg' % SZN)
plt.close()

#########################################################################
#Compare total dynamical change to non-TC (background) dynamical change
h1idir = './linevslat_h1i_0012-0014'
spdvar = 'wspd_bot'
maskvar = 'TC_R4'
h1ispd = xr.open_dataset(os.path.join(h1idir, 'botspd_tot_sznlzm.nc'))[spdvar].sel(season=SZN)
nTCspd = xr.open_dataset(os.path.join(h1idir, 'inv_TC_R4_botspd_sznlzm.nc'))[spdvar].sel(season=SZN)
maskzm = xr.open_dataset(os.path.join(h1idir, 'TC_R4_sznlzm.nc'))[maskvar].sel(season=SZN)

#print(h1ispd.sel(case='UNSEED').data)
nTCspd.loc[dict(case=['UNSEED', 'SEED'])] += nTCspd.sel(case='CTRL')
#print(h1ispd.sel(case='UNSEED').data)

#renormalize non-TC wind speed to isolate how the background flow changes
non_TC_wspd_renorm = nTCspd / (1 - maskzm)
non_TC_wspd_renorm.loc[dict(case=['UNSEED', 'SEED'])] -= non_TC_wspd_renorm.sel(case='CTRL')

dU_frac_h1i = h1ispd.sel(case=['UNSEED', 'SEED']) / h1ispd.sel(case='CTRL')
dU_frac_nTC = non_TC_wspd_renorm.sel(case=['UNSEED', 'SEED']) / non_TC_wspd_renorm.sel(case='CTRL')

plt.rcParams['figure.figsize'] = (15, 4.5)
fig, axes = plt.subplots(1, 2)

for ii, cs in enumerate(['UNSEED', 'SEED']):
   axes[ii].plot(sinlat, orilh * dU_frac.sel(case=cs), label='U10 monthly mean')
   #axes[ii].plot(sinlat, orilh * dU_frac_h1i.sel(case=cs), label='$|U|_{bot}$', color='C0', linestyle='dotted')
   axes[ii].plot(sinlat, orilh * dU_frac_nTC.sel(case=cs), label='$|U|_{bot}$ non-TC', color='C8')

   axes[ii].axhline(0, c='gray', lw=0.5)
   axes[ii].axvspan(*MT, fc='purple', alpha=.08)
   axes[ii].axvspan(*STZ, fc='yellow', alpha=.12)
   [axes[ii].axvline(YSCL(ll), color='gray', linewidth=0.5) for ll in range(-50, 51, 10)]
   axes[ii].legend(fontsize=12)
   axes[ii].set_xticks(YSCL(np.arange(-50, 51, 10)), np.arange(-50, 51, 10))
   axes[ii].tick_params(right=True)

   axes[ii].set_title(['(c)', '(d)'][ii], loc='left')
   axes[ii].set_title('%s$-$CTRL' % cs)
   axes[ii].set_xlabel('Latitude')
   axes[ii].set_ylabel('LHFLX [W m$^{-2}$]')
   axes[ii].set_ylim(*YLIM[ii])

fig.tight_layout()
plt.savefig('./ms_rev_plts/LH_decomp_%s_dyn.svg' % SZN)
plt.show()
