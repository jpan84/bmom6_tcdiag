f=addfiles("/glade/u/home/jpan/grids/b.e23.BMOM.ne30np4_sx0.66av1.aqua.production.250505_internal_heat.mom6.static.nc", "r")
lath = f[:]->geolat
lonh = f[:]->geolon

opt = True
opt@NetCDFType = "netcdf4"
curvilinear_to_SCRIP("/glade/u/home/jpan/grids/sx0.66av1_hgrid_SCRIP.nc", lath, lonh, opt)
