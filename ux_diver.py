import uxarray as ux
import xarray as xr
import numpy as np

FILI = '/glade/derecho/scratch/jpan/jpan_tcfields1/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl.cam.h1i.0014-12-26-00000.nc'
FILO = '/glade/derecho/scratch/jpan/jpan_tcfields1/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl.cam.h1i.0014-12-26-00000.nc.div200.nc'
UGRD = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'

UVAR, VVAR = 'U200', 'V200'

apply_grad = lambda uxda: uxda.gradient()

def main():
   uxds = ux.open_dataset(UGRD, FILI)
   print(uxds[UVAR])
   gradds = uxds['T500'].isel(time=0).gradient()

   print(np.sum(~np.isnan(gradds['meridional_gradient'])))
   print(uxds[UVAR].isel(time=0).divergence(uxds[VVAR].isel(time=0)))
   exit()


   # --- CRITICAL STEP ---
   # Manually trigger the area and centroid calculation once.
   # This ensures the 'uxgrid' object inside uxds is fully populated.
   _ = uxds.uxgrid.face_areas
   _, _ = uxds.uxgrid.face_x, uxds.uxgrid.face_y
   # ---------------------
   grad_tslc = [uxds[UVAR].sel(time=tt).gradient() for tt in uxds['time']]
   grad_ds = xr.concat(grad_tslc, dim=uxds['time'])
   print(grad_ds)
   print(grad_ds.compute())
   exit()

   divu = uxds[UVAR].gradient()['zonal_gradient']
   divv = uxds[VVAR].gradient()['meridional_gradient']

   outds = ux.UxDataset(data_vars={UVAR + '_div': divu, VVAR + '_div': divv})
   outds.to_netcdf(FILO)

def ux_gradient_multdim(uxda, uxgrid):
   global apply_grad
   return xr.apply_ufunc(apply_grad, uxda, input_core_dims=[['n_face']], output_core_dims=[['n_face']],\
               output_dtypes=[uxda.dtype], vectorize=True, dask='parallelized')#, kwargs=dict(uxgrid=uxgrid))

if __name__ == '__main__':
   main()
