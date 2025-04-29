#Joshua Pan jvp5930@psu.edu 20240119
#Adapted from 20221002 version
#Updated to compute zonal accel from EP flux divergence.
#Use data on pressure levels, interpolated using hy2pres.py.
#Compute the transformed Eulerian mean (TEM) streamfunctions based on Eq. 10 from White et al. 2020 (10.1175/JCLI-D-19-0697.1).
#qcmd -q casper -l walltime=00:15:00 -l select=1:ncpus=36:mem=128GB -A UPSU0063 python3 streamf.py 0

#import paths as pt
OUTDIR = './streamf/'

import sys
sys.path.append('/glade/u/home/jpan/aquaptc/aquapgrid/')
import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker as mticker
import consts as c
import linevslat

#GLOBS = [dr[1] for dr in pt.HDIRS(0)]

print(sys.argv)
#DIRIDX = int(sys.argv[1])
DIFF = bool(int(sys.argv[1])) #False if len(sys.argv) < 3 else bool(sys.argv[2])
#DIRIDX2 = None if not DIFF else int(sys.argv[3])
DIR1 = '/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250415_unseed/atm/hist_regrid_onpres/*h0a*.nc'
DIR2 = '/glade/derecho/scratch/jpan/QPC5-ne30np4-aquap10-seed3x3_regrid_onpres/*h0*' if DIFF else None
pres_name = 'plev'

H = 7e3 #m
dens0 = c.p0 / c.g / H
LATLAB = linevslat.LATLAB

def main():
   if not os.path.exists(OUTDIR):
      os.makedirs(OUTDIR)

   print('Reading history files...')
   print(DIR1)
   HIST_DS1 = xr.open_mfdataset(DIR1) #open all h1 files into one dask-chunked array
   streamf1 = comppsi(HIST_DS1)
   #print([(s.lat.min(),s.lat.max()) for s in streamf1])
   plotfields = streamf1
   expo = 10

   if DIFF:
      HIST_DS2 = xr.open_mfdataset(DIR2)
      streamf2 = comppsi(HIST_DS2)
      #print([(s.lat.min(),s.lat.max()) for s in streamf2])
      plotfields = [streamf2[i] - sf1 for i, sf1 in enumerate(streamf1)]
      expo = 9
   
   print('Plotting streamfunctions...')
   plt.rc('font', size=16)
   plt.rcParams['figure.figsize'] = (30, 6)
   contourfkwargs = {'cmap': 'bwr' if DIFF else 'coolwarm', 'levels': 2.**np.arange(-2, 7, 1), 'norm': colors.SymLogNorm(2.**-1)}
   contourfkwargs['levels'] = np.concatenate((-contourfkwargs['levels'][::-1], contourfkwargs['levels']))
   contourkwargs = {'colors': 'black', 'levels': 2.**np.arange(-2, 7, 1)}
   contourkwargs['levels'] = np.concatenate((-contourkwargs['levels'][::-1], contourkwargs['levels']))
   clabelkwargs = {'inline': 1, 'fontsize': 10, 'colors': 'black', 'fmt': '%.1f'}
   subplot_kw = dict(xlim=(-1, 1), ylim=(100, 1000), yscale='log')
   fig, axes = plt.subplots(1, 3, sharey=True, subplot_kw=subplot_kw)
   plotfields = [field / 10**expo for field in plotfields]
   plottitles = ['Eulerian Mean Streamfunction', 'Stokes Drift', 'Residual Streamfunction $\Psi^*$']

   for i, ax in enumerate(axes):
      CSF = ax.contourf(np.sin(np.deg2rad(HIST_DS1.lat)), HIST_DS1[pres_name], plotfields[i].values, **contourfkwargs)
      CS1 = ax.contour(np.sin(np.deg2rad(HIST_DS1.lat)), HIST_DS1[pres_name], plotfields[i].values, **contourkwargs) 
      if i == 0:
         ax.yaxis.set_minor_formatter(mticker.ScalarFormatter())
         ax.yaxis.set_major_formatter(mticker.ScalarFormatter())
         ax.set_ylabel('Pressure [hPa]')
         ax.set_yticks(np.arange(100, 1001, 100))
         plt.colorbar(CSF, ax=axes)
      ax.set_xticks(np.sin(np.deg2rad(LATLAB)), labels=LATLAB)
      ax.set_xlim(np.sin(np.deg2rad(-15)), np.sin(np.deg2rad(90)))
      plt.gca().invert_yaxis()
      #ax.clabel(CS1, **clabelkwargs)
      ax.set_title('%s (10$^{%d}$ kg s$^{-1}$)' % (plottitles[i], expo))
      ax.set_xlabel('Latitude [°]')

   plt.savefig(os.path.join(OUTDIR, '%s_%s_TEMstreamf.png' % (DIR1.split('/')[-4], DIR2.split('/')[-4] if DIFF else '')), bbox_inches='tight')
   plt.close()

   #print('Saving to .nc...')
   #outds = xr.Dataset(data_vars = {'PSI_EM': PSI_EM, 'PSI_stokes': PSI_bal, 'PSI_resid': PSI_resid, 'EPy': Fy, 'EPz': Fz, 'EPylp': Fylp, 'EPzlp': Fzlp, 'EPyhp': Fyhp, 'EPzhp': Fzhp, 'EPdiv': dF})
   #outds.to_netcdf(path = '%s_%s_TEM.nc' % plotfileargs)

   print('Done.')

def comppsi(HIST_DS):
   print('Computing Eulerian mean streamfunction...')
   #compute the streamfunction assoc w/ Eulerian-mean V
   V_ZM = HIST_DS.V.mean(dim='lon')
   V_ZMEM = V_ZM.mean(dim='time')
   coslat = np.cos(HIST_DS.lat * np.pi / 180)
   latcirc = 2 * np.pi * c.a * coslat #circumference at each latitude
   integrand = V_ZMEM * latcirc / c.g
   PSI_EM = trapint(integrand, HIST_DS[pres_name]) * 100

   print('Computing Stokes drift...')
   VTH_ZMEM = thta(HIST_DS.VT, HIST_DS[pres_name]).mean(dim=['time', 'lon']) 
   TH_MEAN = thta(HIST_DS.T, HIST_DS[pres_name]).mean(dim=['time', 'lon'])
   V_MEAN = HIST_DS.V.mean(dim=['time', 'lon'])
   VTH_MEAN = V_MEAN * TH_MEAN
   EHF = VTH_ZMEM - VTH_MEAN #- VTH_STN

   dTHTA_dp = TH_MEAN.differentiate(pres_name, edge_order=2) / 100 #K Pa-1
   #plt.boxplot(dTHTA_dp.mean(dim='time').values, showfliers=False)
   #plt.savefig('%s_dTHTA_dp.png' % CASENAME)
   #plt.close()
   dTHTA_dp = dTHTA_dp.clip(max=-1e-4) #takes care of convective instability (or low stability) in TTL
   PSI_bal = (-EHF * latcirc / c.g / dTHTA_dp)#.mean(dim='time')

   PSI_resid = PSI_EM + PSI_bal

   return [PSI_EM, PSI_bal, PSI_resid]

#compute the trapezoidal integral given integrand values and 1D coordinate values
#start at a zero integral value for the 1st coord value
def trapint(igrnd, coords):
   dim = list(coords.coords)[0]
   #sz = coords.shape[0]
   igral = xr.DataArray(np.zeros(igrnd.shape, dtype=np.float64), dims=igrnd.dims, coords=igrnd.coords)
   for i, pt in enumerate(coords[1:]):
      prev = {dim: i} #use [] to index DataArrays
      curr = {dim: pt} #use .loc[] to index DataArrays
      step = pt - coords[prev]
      igral.loc[curr] = igral[prev] + (igrnd.loc[curr] + igrnd[prev]) / 2 * step
   return igral

def thta(T, p):
   """
   Returns potential temperature given temp [K] and pressure [hPa].
   """
   return T * (1000 / p)**(c.Rd/c.cp)


if __name__ == '__main__':
   main()
