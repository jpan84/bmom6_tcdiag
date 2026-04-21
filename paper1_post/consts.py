#from atm.log
g = 9.797640000
OM = 7.292115000e-5
a_e = 6.371000000e6
cp = 1004.64000000000
MWDRY = 0.2896623325e-1
MWH2O = 0.1801618113e-1

P0 = 1e5

#need to double check below
lv = 2500840.
rho_w = 1000.
Runi = 8.31432

Rd = Runi / MWDRY
Rv = Runi / MWH2O
kapd = Rd / cp
