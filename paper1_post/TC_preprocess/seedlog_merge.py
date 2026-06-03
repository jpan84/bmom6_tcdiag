#from seedmatch_stats import dfcols as mscols
#from seed_stats import dfcols as xscols
#from unseed_stats import dfcols as uscols
#
#print(set(mscols).union(set(xscols), set(uscols)))
#{'clon', 'sstval', 'exppr', 'tid', 'rmw_target', 'zp', 'rp', 'psamb', 'minsep', 'dt', 'psmin', 'dp', 'sstlat', 'dp [hPa]', 'clat', 'rmw_final'}

#construct a pandas dataframe of un/seed events
import os
import sys
import glob
import re
import numpy as np
import pandas as pd
from datetime import datetime as dt

alias = sys.argv[1] #'251229_seed_match'
MSTRCOLS = ['dt', 'sstlat', 'sstval', 'tid', 'clat', 'clon', 'minsep', 'psamb', 'psmin', 'dp', 'rmw_target', 'rmw_final', 'rp', 'zp', 'exppr']
STNG = ['rp', 'dp', 'zp', 'exppr']
logpaths = sorted(glob.glob('/glade/u/home/jpan/work/MOM6_CASEDIRS/%s.out*' % alias))

ps_pattern_1000 = r"^\[1\d{5}\."
ps_pattern_900 = r"^9\d{4}\."

def main():
   seedlog = list(open(logpaths[0], 'r'))
   for pt in logpaths[1:]:
      seedlog.extend(list(open(pt, 'r')))

   outdf = pd.DataFrame(columns=MSTRCOLS)
   data = {col: np.nan for col in MSTRCOLS}

   mkplt, prevln, dtstr, settings = False, '', None, None

   for ln in seedlog:
      spl = ln.strip('\n').split(' ')
      if not len(spl[0]):
         continue
      if 'not in same hemisphere' in ln and ln != prevln:
         data['psmin'], settings = np.nan, []
         mkplt = True
      if spl[0] == 'find-sst-max.py:':
         data['sstlat'] = float(spl[4])
         data['sstval'] = float(spl[-1])
      if spl[0] in ['sedding', 'lat/lon']:
         #211 sedding lat/lon: -8.347259 48.750000
         data['clon'] = float(spl[-1])
         data['clat'] = float(spl[-2])
      elif spl[:2] == ['Debugging', 'restart']:
         #212 Debugging restart file /glade/derecho/scratch/jpan/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250129_unseed_all/run/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250129_unseed_all.cam.r.0001-02-02-00000.nc at time 5407747
         data['tid'] = int(spl[-1])
      elif spl[0][0] == "'":
         #200 '/glade/derecho/scratch/jpan/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250129_unseed_all/run/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250129_unseed_all.cam.r.0001-02-02-00000.nc' -> '/glade/derecho/scratch/jpan/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250129_unseed_all/run/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250129_unseed_all.cam.r.0001-02-02-00000.nc.ORIG.nc'
         dtstr = re.search(r"(\d{4}-\d{2}-\d{2}-\d{5})", spl[0]).group()
         print(dtstr)
      elif spl[:2] == ['minimum', 'separation']:
         data['minsep'] = float(spl[-1])
      elif spl[0][0] == '[':
         if re.match(ps_pattern_1000, spl[0]): 
            data['psmin'] = float(spl[0].strip('[')) / 100
         if spl[0] == '[' and re.match(ps_pattern_900, spl[1]):
            data['psmin'] = float(spl[1]) / 100
      elif spl[0] == 'ambient':
         data['psamb'] = float(spl[-1])
      elif spl[0] == 'iter:':
         #print(spl)
         data['rmw_target'] = float(spl[-4])
      elif spl[:2] == ['random', 'minp:']:
         data['dp'] = float(spl[-1])
      elif spl[:2] == ['random', 'zp:']:
         data['zp'] = float(spl[-1])
      elif spl[:2] == ['random', 'exppr:']:
         data['exppr'] = float(spl[-1])
      elif spl[0] == 'rp:':
         data['rp'] = float(spl[1])
         #mkplt = True
      elif prevln[:5] == 'iter:' and spl[0] != 'iter:':
         data['rmw_final'] = float(prevln.strip('\n').split(' ')[7])

      if 'Final BEST SETTINGS' in ln:
         settings = [float(num) for num in re.compile(r'\d+\.\d+').findall(ln)]
         for jj, param in enumerate(STNG):
            data[param] = settings[jj]
         mkplt = True #make plot now that we have all the info for 1 storm
      if spl[:3] == ['Average', 'PSDRY', 'out']:
         mkplt = True

      #TODO: dtstr to dtobj?
      if mkplt:
         data['dt'] = dt.strptime('-'.join(dtstr.split('-')[:-1]), '%Y-%m-%d')
         outdf.loc[len(outdf)] = data

         mkplt, dtstr = False, None
         data = {col: np.nan for col in MSTRCOLS}
         settings = None if settings == None else []

      prevln = ln

   # Save to dynamic files out
   os.makedirs('seed_stats', exist_ok=True)
   outdf.to_parquet(f'seed_stats/{alias}_events.parquet')
   outdf.to_csv(f'seed_stats/{alias}_events.csv')
   print(f"Processed {alias} into master schema shape: {outdf.shape}")

if __name__ == '__main__':
    main()
