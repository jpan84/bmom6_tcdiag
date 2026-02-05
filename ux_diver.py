import uxarray as ux
import xarray as xr

FILI = '/glade/derecho/scratch/jpan/jpan_tcfields1/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl.cam.h1i.0014-12-26-00000.nc'
FILO = '/glade/derecho/scratch/jpan/jpan_tcfields1/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl.cam.h1i.0014-12-26-00000.nc.div200.nc'
UGRD = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'

UVAR, VVAR = 'U200', 'V200'

apply_grad = lambda uxda: uxda.gradient()

def main():
   uxds = ux.open_mfdataset(UGRD, FILI)
   #out = xr.apply_ufunc(ux.UxDataArray.gradient, uxds[UVAR], input_core_dims=[['n_face']], output_core_dims=[['n_face']],\
   #            output_dtypes=[uxds[UVAR].dtype], vectorize=True, dask='parallelized', kwargs=dict(uxgrid=uxds.uxgrid))
   print(ux_gradient_multdim(uxds[UVAR], uxds.uxgrid).compute())
   #u_gradient = uxds[UVAR].groupby('time').map(lambda x: x.gradient())
   #print(u_gradient)
   #print(u_gradient.compute())
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
