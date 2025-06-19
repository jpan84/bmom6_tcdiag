import uxarray as ux
import matplotlib.pyplot as plt

g = 9.81
cp = 1004
lv = 2500840

testfile = '/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl/atm/hist_0010_h1i/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl.cam.h1i.0010-07-01-21600.nc'
camgrid = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'

def main():
   ds = ux.open_dataset(camgrid, testfile)
   print(ds['area'])

   mse850 = g * ds['Z850'] + cp * ds['T850'] + lv * ds['Q850']
   plt.hist(mse850.data)
   plt.show()

   umf500 = (ds['OMEGA500'] < 0) * (-ds['OMEGA500'] * ds['area'] / g)

if __name__ == '__main__':
   main()
