import uxarray as ux
import numpy as np
import holoviews as hv
import os
import glob
from paths import ARCHRT, CAMGR

DIRS = [ARCHRT[1], ARCHRT[3]]
FILS = [['0005-08-06-64800', '0005-08-07-00000', '0005-08-07-21600', '0005-08-07-43200'],
        ['0006-09-29-64800', '0006-09-30-00000', '0006-09-30-21600', '0006-09-30-43200']]

CTRS = [(98.79, 18.01), (22.73, 20.78)]
DEGW = 7.5

PLEVS = np.concatenate((np.arange(980, 1005, 5), np.arange(1005, 1015, 1)))

def main():
   pths = [glob.glob(os.path.join(ar, 'atm/hist/', f'*h1i.{fi}*.nc') for fi in FILS[ii]]) for ii, ar in enumerate(DIRS)]

   dss = [ux.open_mfdataset(CAMGR, [pt for pt in pths[ii]]) for ii, ar in enumerate(DIRS)]

   subs = [ds['PSL'].subset.bounding_box((CTRS[ii][0] - DEGW, CTRS[ii][0] + DEGW), (CTRS[ii][1] - DEGW, CTRS[ii][1] + DEGW)) for ii, ds in enumerate(dss)]

   plts = np.array([[(sb.sel(time=tt) / 100).plot(rasterize=True) for tt in sb['time']] for sb in subs])

   grid_layout = hv.Layout(plts.flatten()).cols(4)

   # OPTION A: Save as an interactive HTML file (Best for sharing/viewing in a browser)
   hv.save(grid_layout, 'seed_evo.html', backend='bokeh')
   print("Interactive plot saved to 'seed_evo.html'")

   # OPTION B: Show interactively in a window (Requires an X11 server running if remote)
   # import bokeh.io
   # bokeh.io.show(hv.render(grid_layout, backend='bokeh'))

if __name__ == "__main__":
   main()
