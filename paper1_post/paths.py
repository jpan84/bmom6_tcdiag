import os

ARCHRT = ['/glade/campaign/univ/upsu0032/jpan_aquaptc/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250702_unseed2hPa6m/', '/glade/campaign/univ/upsu0032/jpan_aquaptc/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250415_unseed/', '/glade/campaign/univ/upsu0032/jpan_aquaptc/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl/', '/glade/derecho/scratch/jpan/archive/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.251229_seedmatch/', '/glade/campaign/univ/upsu0032/jpan_aquaptc/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250416_seed1x1/']
CASENAMES = [os.path.basename(os.path.normpath(ar)) for ar in ARCHRT]
ALIA = ['UNSEED_EX', 'UNSEED', 'CTRL', 'SEED', 'SEED_EX']
CTLIX = 2
IXHORS = {0: 1, 1: 2, 2: 3, 3: 0, 4: 2, 5: 4}
CAMGR = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'
OCN_CUST = [True, True, False, True, True]
