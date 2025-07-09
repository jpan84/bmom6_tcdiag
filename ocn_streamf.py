import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as colors

TESTFILE = '/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl/ocn/hist/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl.mom6.hm.0008-07.nc'

a = 6.371e6
g = 9.81
dens0 = 1035 #Boussinesq background density

def main():
   ds = xr.open_dataset(TESTFILE)
   Psi_EM = ds['vmo'].sum(dim='xh').cumsum(dim='zl').mean(dim='time')
   dlon = np.deg2rad(ds['xh'][1] - ds['xh'][0])
   print(dlon)

   expo = 10
   contourfkwargs = {'cmap': 'coolwarm', 'levels': 2.**np.arange(-2, 7, 1), 'norm': colors.SymLogNorm(2.**-1)}
   contourfkwargs['levels'] = np.concatenate((-contourfkwargs['levels'][::-1], contourfkwargs['levels']))
   contourkwargs = {'colors': 'black', 'levels': 2.**np.arange(-2, 7, 1)}
   contourkwargs['levels'] = np.concatenate((-contourkwargs['levels'][::-1], contourkwargs['levels']))

   csf = plt.contourf(Psi_EM.yq, Psi_EM.zl, Psi_EM.data / 10**expo, **contourfkwargs)
   plt.contour(Psi_EM.yq, Psi_EM.zl, Psi_EM.data / 10**expo, **contourkwargs)
   plt.yscale('log')
   plt.ylabel('Depth [m]')
   plt.xlabel('lat')
   plt.title('Eulerian mean mass streamfunc [1e10 kg/s]')
   plt.gca().invert_yaxis()
   plt.colorbar(csf)
   #plt.show()
   plt.close()

   coslat = np.cos(np.deg2rad(ds['yq'])) #TODO: check numerical accuracy/validity of shifting yh to yq
   print(ds['vmo'].max())
   vmf = (ds['vmo'] / a / coslat / dlon / ds['h'].data).mean(dim=['xh', 'time'])
   print(vmf.max())
   rho_yq = ds['rhoinsitu'].mean(dim=['xh', 'time']).interp(yh=ds['yq'], method='linear')
   vrho_mean = ds['vo'].mean(dim=['xh', 'time']) * rho_yq.data
   vpbp = -g / dens0 * (vmf - vrho_mean)
   bmean = -g / dens0 * (ds['rhoinsitu'].mean(dim=['xh', 'time']) - dens0)
   bz = -bmean.differentiate('zl', edge_order=2)
   coslat_h = np.cos(np.deg2rad(ds['yh']))
   by = (bmean * coslat_h).differentiate('yh', edge_order=2) * 180 / np.pi / a / coslat_h
   kapzm = ds['diftrblo'].mean(dim=['xh', 'time'])
   #plt.hist(bz.data)
   #plt.show()
   #Psi_qs = (2 * np.pi * a * coslat * dens0 * vpbp).transpose('zl', ...) / bz.data
   Psi_qs = (2 * np.pi * a * coslat_h * dens0 * kapzm * by / bz).transpose('zl', ...)
   print(vpbp.shape)
   print(bz.data.shape)
   print(Psi_qs.shape)
   #Psi_qs -= Psi_qs.isel(zl=-1)

   #print(ds['rhoinsitu'].mean())
   #print(vmf.dims)
   #print(vrho_mean.dims)
   #print(vpbp.dims)
   #print(vpbp)
   #print(bz)
   #print(coslat)
   #print(Psi_qs.zl)
   #print(ds['h'])

   csf = plt.contourf(Psi_qs.yh, Psi_qs.zl, Psi_qs.data / 10**expo, **contourfkwargs)
   plt.contour(Psi_qs.yh, Psi_qs.zl, Psi_qs.data / 10**expo, **contourkwargs)
   plt.yscale('log')
   plt.ylabel('Depth [m]')
   plt.xlabel('lat')
   plt.title('Quasi-Stokes mass streamfunc [1e10 kg/s]')
   plt.gca().invert_yaxis()
   plt.colorbar(csf)
   plt.show()
   plt.close()

   ###levels = 10.**np.arange(-2, 4)
   ###levels = np.concatenate((-levels[::-1], levels))
   ###csf = plt.contourf(Psi_qs.yq, Psi_qs.zl, vrho_mean, levels=levels, norm=colors.SymLogNorm(0.01), cmap='coolwarm')#Psi_qs.data / 10**expo, **contourfkwargs) )#
   ###plt.contour(Psi_qs.yq, Psi_qs.zl, vrho_mean, levels=levels, colors='black')
   ###plt.title('Mass flux by Eulerian mean flow vbar * rhobar [kg m-2 s-1]')
   ###plt.yscale('log')
   ###plt.ylabel('Depth [m]')
   ###plt.xlabel('lat')
   ###plt.gca().invert_yaxis()
   ###plt.colorbar(csf)
   ####plt.show()
   ###plt.close()

   ###levels = 10.**np.arange(-2, 4)
   ###levels = np.concatenate((-levels[::-1], levels))
   ###csf = plt.contourf(Psi_qs.yq, Psi_qs.zl, vmf, levels=levels, norm=colors.SymLogNorm(0.01), cmap='coolwarm')#Psi_qs.data / 10**expo, **contourfkwargs) )#
   ###plt.contour(Psi_qs.yq, Psi_qs.zl, vmf, levels=levels, colors='black')
   ###plt.title('Total mass flux vmo / (dx * dz) [kg m-2 s-1]')
   ###plt.yscale('log')
   ###plt.ylabel('Depth [m]')
   ###plt.xlabel('lat')
   ###plt.gca().invert_yaxis()
   ###plt.colorbar(csf)
   ####plt.show()
   ###plt.close()

   ###csf = plt.contourf(Psi_qs.yq, Psi_qs.zl, vpbp, cmap='coolwarm')#Psi_qs.data / 10**expo, **contourfkwargs) )#
   ####plt.contour(Psi_qs.yq, Psi_qs.zl, vpbp, levels=levels, colors='black')
   ###plt.title('Eddy buoyancy flux v\'b\' [m-2 s-3]')
   ###plt.yscale('log')
   ###plt.ylabel('Depth [m]')
   ###plt.xlabel('lat')
   ###plt.gca().invert_yaxis()
   ###plt.colorbar(csf)
   ####plt.show()
   ###plt.close()

   ###csf = plt.contourf(Psi_qs.yq, Psi_qs.zl, vpbp * dens0 / g, cmap='coolwarm')#Psi_qs.data / 10**expo, **contourfkwargs) )#
   ####plt.contour(Psi_qs.yq, Psi_qs.zl, vpbp, levels=levels, colors='black')
   ###plt.title('Eddy mass flux v\'rho\' [kg m-2 s-1]')
   ###plt.yscale('log')
   ###plt.ylabel('Depth [m]')
   ###plt.xlabel('lat')
   ###plt.gca().invert_yaxis()
   ###plt.colorbar(csf)
   ####plt.show()
   ###plt.close()

   ###zscl = np.vectorize(lambda zl: zl / 2850 if zl <= 2000 else 2000 / 2850 + (zl - 2000) / 2000 * 850 / 2850)
   ###zlab = np.arange(1000, 5000, 1000)
   ###zloc = zscl(zlab)
   ###csf = plt.contourf(ds.yh, zscl(ds.zl), ds['rhoinsitu'].mean(dim=['xh', 'time']), cmap='viridis', levels=np.arange(1020, 1040.1, 0.5))
   ####plt.contour(Psi_qs.yq, Psi_qs.zl, vpbp, levels=levels, colors='black')
   ###plt.title('Zonal mean density [kg m-3]')
   ####plt.yscale('log')
   ###plt.yticks(zloc, zlab)
   ###plt.ylabel('Depth [m]')
   ###plt.xlabel('lat')
   ###plt.gca().invert_yaxis()
   ###plt.colorbar(csf)
   ###plt.show()
   ###plt.close()

   ###csf = plt.contourf(Psi_qs.yq, Psi_qs.zl, bz, levels=np.arange(0, 1.6e-4, 3e-5))
   ####plt.contour(Psi_qs.yq, Psi_qs.zl, vpbp, levels=levels, colors='black')
   ###plt.title('Buoyancy freq squared [s-2]')
   ###plt.yscale('log')
   ###plt.ylabel('Depth [m]')
   ###plt.xlabel('lat')
   ###plt.gca().invert_yaxis()
   ###plt.colorbar(csf)
   ####plt.show()
   ###plt.close()

   ###vbar = ds['vo'].mean(dim=['xh', 'time'])
   ###csf = plt.contourf(Psi_qs.yq, Psi_qs.zl, vbar, cmap='coolwarm', levels=np.arange(-0.2, 0.21, .04))
   ####plt.contour(Psi_qs.yq, Psi_qs.zl, vpbp, levels=levels, colors='black')
   ###plt.title('V velocity [m s-1]')
   ###plt.yscale('log')
   ###plt.ylabel('Depth [m]')
   ###plt.xlabel('lat')
   ###plt.gca().invert_yaxis()
   ###plt.colorbar(csf)
   ####plt.show()
   ###plt.close()

   ###Psi_qs_vel = (vpbp / bz.data).differentiate('zl', edge_order=2)
   ###csf = plt.contourf(Psi_qs.yq, Psi_qs.zl, Psi_qs_vel, cmap='coolwarm', levels=np.arange(-0.2, 0.21, .04))
   ####plt.contour(Psi_qs.yq, Psi_qs.zl, vpbp, levels=levels, colors='black')
   ###plt.title('V QS velocity [m s-1]')
   ###plt.yscale('log')
   ###plt.ylabel('Depth [m]')
   ###plt.xlabel('lat')
   ###plt.gca().invert_yaxis()
   ###plt.colorbar(csf)
   ####plt.show()
   ###plt.close()

   Psi_resid = Psi_EM + Psi_qs.interp(yh=Psi_EM['yq'].data).data
   csf = plt.contourf(Psi_EM.yq, Psi_EM.zl, Psi_resid.data / 10**expo, **contourfkwargs)
   plt.contour(Psi_EM.yq, Psi_EM.zl, Psi_resid.data / 10**expo, **contourkwargs)
   plt.yscale('log')
   plt.ylabel('Depth [m]')
   plt.xlabel('lat')
   plt.title('Residual mass streamfunc [1e10 kg/s]')
   plt.gca().invert_yaxis()
   plt.colorbar(csf)
   plt.show()
   plt.close()

   csf = plt.contourf(Psi_EM.yq, Psi_EM.zl, Psi_EM.data / 10**expo, **contourfkwargs)
   plt.contour(Psi_EM.yq, Psi_EM.zl, Psi_EM.data / 10**expo, **contourkwargs)
   plt.yscale('log')
   plt.ylabel('Depth [m]')
   plt.xlabel('lat')
   plt.title('Eulerian mass streamfunc [1e10 kg/s]')
   plt.gca().invert_yaxis()
   plt.colorbar(csf)
   plt.show()
   plt.close()


if __name__ == '__main__':
   main()
