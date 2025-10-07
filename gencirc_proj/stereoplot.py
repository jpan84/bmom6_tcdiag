#Joshua Pan jp872@cornell.edu 20220624
#Updated to take existing stereo axes as kwarg
#Set up North Pole stereographic plotting

import cartopy.crs as ccrs
import matplotlib.path as mpath
import numpy as np
import matplotlib.pyplot as plt
from time import perf_counter_ns

#return axes
#Args: southernmost latitude of plot, ax kwarg should be an existing axes instance with projection ccrs.NorthPolarStereo()
def NPstereoaxes(southbound, ax=None):
   if ax == None:
      ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=0.0), label=str(int(perf_counter_ns())))
   ax.set_extent([-180, 180, southbound, 90], ccrs.PlateCarree())
   # Compute a circle in axes coordinates, which we can use as a boundary
   # for the map. We can pan/zoom as much as we like - the boundary will be
   # permanently circular.
   theta = np.linspace(0, 2*np.pi, 100)
   center, radius = [0.5, 0.5], 0.5
   verts = np.vstack([np.sin(theta), np.cos(theta)]).T
   circle = mpath.Path(verts * radius + center)
   ax.set_boundary(circle, transform=ax.transAxes)
   ax.gridlines(crs=ccrs.PlateCarree(), draw_labels = True, x_inline = True, y_inline = True) #unsure how to properly label gridlines (labels not showing)
   return ax
