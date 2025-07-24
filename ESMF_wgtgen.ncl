sgrid = "/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc"
dgrid = "/glade/u/home/jpan/grids/sx0.66av1_hgrid_SCRIP.nc"
wgtf = "/glade/u/home/jpan/grids/map_ne120np4_to_sx0.66av1_hgrid.stod.250724.nc"

opt = True
opt@InterpMethod = "neareststod"
opt@DstRegional = True
opt@NetCDFType = "netcdf4"
opt@Debug = True

ESMF_regrid_gen_weights(sgrid, dgrid, wgtf, opt)
