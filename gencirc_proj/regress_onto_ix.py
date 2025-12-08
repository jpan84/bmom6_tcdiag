import glob
import os
from datetime import timedelta as tdel
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm

IXFIL = './0012-0013_JJASON_EEHF_ix_std.nc'
EVFIL = './ehf_events_sep3_sig1.0.nc'
IXVAR = 'eehf_ix_std'
FLDDIR = '/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.251010_ctrlbr/atm'
FLDFIL = '0012-0013_JJASON_onpres_1deg_zm1_anom.nc'
VARFIL = '0012-0013_JJASON_onpres_1deg_zonvar_anom.nc'
FLDAVG = '0012-0013_JJASON_onpres_1deg_zm1_tm.nc'
VARAVG = '0012-0013_JJASON_onpres_1deg_zonvar_tm.nc'

ixds = xr.open_dataset(IXFIL).squeeze()
evds = xr.open_dataset(EVFIL)
fldds = xr.open_dataset(os.path.join(FLDDIR, FLDFIL)).squeeze()
vards = xr.open_dataset(os.path.join(FLDDIR, VARFIL)).squeeze()
avfds = xr.open_dataset(os.path.join(FLDDIR, FLDAVG)).squeeze()
avvds = xr.open_dataset(os.path.join(FLDDIR, VARAVG)).squeeze()

reg_on_pc = lambda vec: (ixds[IXVAR] * vec).sum(dim='time') / (ixds[IXVAR]**2).sum(dim='time') #anom of physical quantity corresponding to +1std anom in index

meands = avfds
anomds = fldds
var2d = 'FLUT'
laghrs = np.arange(-336, 337, 48.)
lagtds = [tdel(hours=hh) for hh in laghrs]
cnorm = colors.Normalize(vmin=-laghrs.max(), vmax=-laghrs.min())
cmap = cm.get_cmap('managua_r')
mpbl = cm.ScalarMappable(norm=cnorm, cmap=cmap)
mpbl.set_array([])
plt.rc('font', size=16)
plt.plot(meands['lat'], meands[var2d], color='black', lw=2.5)
plt.xlim(-10, 40)
plt.xlabel('Lat')
plt.ylabel(var2d + ' [%s]' % meands[var2d].attrs['units'])
ax1 = plt.gca().twinx()
for ii, lh in enumerate(laghrs):
   compotimes = np.intersect1d((evds['peaks'] - lagtds[ii]), anomds[var2d]['time'], assume_unique=True)
   compo = anomds[var2d].sel(time=compotimes).mean(dim='time')
   #ax1.plot(anomds['lat'], reg_on_pc(anomds[var2d].assign_coords(time=anomds['time'] + lagtds[ii])),\
   #             label='%dh' % -lh, c=cmap(cnorm(-lh))) #negative sign so that lag < 0 means before EEHF index peaks
   ax1.plot(anomds['lat'], compo, c=cmap(cnorm(-lh)), ls='dashed')
ax1.axhline(0, lw=0.5, c='gray')
ax1.set_ylim(-4, 4)
cb = plt.colorbar(mpbl, ax=plt.gca(), pad=0.1)
cb.set_label('lag [days]')
cb.ax.set_yticks(-laghrs, -(laghrs / 24).astype(np.int_))
plt.savefig('OLRzm_reg_on_ix.png', bbox_inches='tight')
plt.show()
plt.close()

autoreg = []
for ii, lh in enumerate(laghrs):
   autoreg.append(reg_on_pc(ixds[IXVAR].assign_coords(time=ixds['time'] + lagtds[ii])))
plt.plot(-laghrs, autoreg)
plt.xlim(0, 340)
plt.xlabel('Lag [hours]')
plt.ylabel('EEHF autocorrelation')
plt.savefig('EEHF_ix_autocor.png', bbox_inches='tight')
plt.close()
#plt.show()

meands = avfds
anomds = fldds
var2d = 'vIVT'
laghrs = np.arange(-336, 337, 48.)
lagtds = [tdel(hours=hh) for hh in laghrs]
cnorm = colors.Normalize(vmin=-laghrs.max(), vmax=-laghrs.min())
cmap = cm.get_cmap('managua_r')
mpbl = cm.ScalarMappable(norm=cnorm, cmap=cmap)
mpbl.set_array([])
plt.rc('font', size=16)
plt.plot(meands['lat'], meands[var2d], color='black', lw=2.5)
plt.xlim(-10, 40)
plt.xlabel('Lat')
plt.ylabel(var2d + ' [%s]' % meands[var2d].attrs['units'])
ax1 = plt.gca().twinx()
for ii, lh in enumerate(laghrs):
   ax1.plot(anomds['lat'], reg_on_pc(anomds[var2d].assign_coords(time=anomds['time'] + lagtds[ii])),\
                label='%dh' % -lh, c=cmap(cnorm(-lh)))
ax1.axhline(0, lw=0.5, c='gray')
#ax1.set_ylim(-2.5, 2.5)
cb = plt.colorbar(mpbl, ax=plt.gca(), pad=0.1)
cb.set_label('lag [days]')
cb.ax.set_yticks(-laghrs, -(laghrs / 24).astype(np.int_))
plt.savefig('vIVTzm_reg_on_ix.png', bbox_inches='tight')
#plt.show()
plt.close()

meands = avvds
anomds = vards
var2d = 'TMQ'
laghrs = np.arange(-336, 337, 48.)
lagtds = [tdel(hours=hh) for hh in laghrs]
cnorm = colors.Normalize(vmin=-laghrs.max(), vmax=-laghrs.min())
cmap = cm.get_cmap('managua_r')
mpbl = cm.ScalarMappable(norm=cnorm, cmap=cmap)
mpbl.set_array([])
plt.rc('font', size=16)
plt.plot(meands['lat'], meands[var2d], color='black', lw=2.5)
plt.xlim(-10, 40)
plt.xlabel('Lat')
plt.ylabel(var2d + ' [%s]' % meands[var2d].attrs['units'])
ax1 = plt.gca().twinx()
for ii, lh in enumerate(laghrs):
   ax1.plot(anomds['lat'], reg_on_pc(anomds[var2d].assign_coords(time=anomds['time'] + lagtds[ii])),\
                label='%dh' % -lh, c=cmap(cnorm(-lh)))
ax1.axhline(0, lw=0.5, c='gray')
#ax1.set_ylim(-2.5, 2.5)
cb = plt.colorbar(mpbl, ax=plt.gca(), pad=0.1)
cb.set_label('lag [days]')
cb.ax.set_yticks(-laghrs, -(laghrs / 24).astype(np.int_))
plt.savefig('CWVvar_reg_on_ix.png', bbox_inches='tight')
plt.show()
plt.close()
