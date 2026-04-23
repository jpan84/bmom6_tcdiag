import pandas as pd
import numpy as np
import sys
import cftime
import datetime as dt

SEDFIL = sys.argv[1]
PARFIL = sys.argv[2]
LAT1 = float(sys.argv[3]) #lower bound (positive, inclusive)
LAT2 = float(sys.argv[4]) #upper bound (positive, exclusive)
dt_fmt = '%Y-%m-%d %H:%M:%S'
print(sys.argv, '\n')

SELLAG = 12.
SELRNG = 3.0
lags = np.arange(6, 19, 6.)
rngs = np.arange(2, 4.1, 0.5)
#lags = np.arange(6, 7, 6.)
#rngs = np.arange(1, 3, 1.)
tab1 = np.zeros((lags.size, rngs.size)) #count exactly one match
tab2 = np.zeros(tab1.shape) #more than one match

def main():
   sdf = pd.read_csv(SEDFIL)
   ori_pdf = pd.read_csv(PARFIL)
   ori_pdf['is_seeded'] = False
   #print(sdf.head())
   #print(pdf.head())

   sdf['dt'] = sdf['dt'].apply(cftime.datetime.strptime, args=(dt_fmt,), calendar='noleap')
   sdf = sdf[(sdf['clat'] >= LAT1) & (sdf['clat'] < LAT2) | (sdf['clat'] <= -LAT1) & (sdf['clat'] > -LAT2)]

   pdf = ori_pdf[ori_pdf['isgen']]
   try:
      pdf['dt'] = pdf['dt'].apply(cftime.datetime.strptime, args=(dt_fmt,), calendar='noleap')
   except ValueError: #tc_stats post-processing accidentally interpolated leap days
      pdf['dt'] = pdf['dt'].apply(cftime.datetime.strptime, args=(dt_fmt,), calendar='Gregorian')
      #print(pdf['dt'].apply(lambda x: x.day))
      pdf = pdf[~pdf['dt'].apply(lambda x: x.month == 2 and x.day == 29)]
      pdf['dt'] = pdf['dt'].apply(lambda d: cftime.datetime(d.year, d.month, d.day, d.hour, calendar='noleap'))
   pdf['dt'] = pdf['dt'] - dt.timedelta(days=2000 * 365)

   #for ix, rw in sdf.iterrows(): 
   #   #print(type(rw['dt']))
   #   dtmatch = pdf[pdf['dt'] == rw['dt']]
   #   #print(np.deg2rad(dtmatch['lon']))
   #   print(rw['clon'], rw['clat'])
   #   print(dtmatch)
   #   print(gcd_deg(rw['clon'], rw['clat'], dtmatch['lon'], dtmatch['lat']))
   #   break

   for ii, ll in enumerate(lags):
      for jj, rg in enumerate(rngs):
         print('Working on lag, range', ll, rg)
         for ix, rw in sdf.iterrows():
            if rw['dt'].day == 1 and rw['dt'].hour == 0:
               print('\t', rw['dt'])
            dtmatch = pdf[(pdf['dt'] - rw['dt'] <= dt.timedelta(hours=ll)) & (pdf['dt'] - rw['dt'] >= dt.timedelta(hours=lags[0]))]
            #print(rw['dt'])
            #print(dtmatch)
            inrng = dtmatch[gcd_deg(rw['clon'], rw['clat'], dtmatch['lon'], dtmatch['lat']) <= rg]
            if inrng.shape[0] == 1:
               tab1[ii, jj] += 1
            elif inrng.shape[0] > 1:
               tab2[ii, jj] += 1

            if ll == SELLAG and rg == SELRNG:
               if inrng.shape[0] == 1:               
                  ori_pdf.loc[ori_pdf['stmnum'].isin(inrng['stmnum']), 'is_seeded'] = True
               elif inrng.shape[0] > 1:
                  import warnings
                  warnings.warn('More than one storm has been linked to a single seeding event.')
                  exit()

   print(tab1 / sdf.shape[0])
   print('\n')
   print(tab2 / sdf.shape[0])

   out_path = PARFIL.replace('.csv', '.seedflagged.')
   ori_pdf.to_csv(out_path + 'csv', index=False)
   ori_pdf.to_parquet(out_path + 'parquet', index=False)

def gcd_deg(clon, clat, ser_lon, ser_lat):
   clonr, clatr = np.deg2rad(clon), np.deg2rad(clat)
   ser_lonr, ser_latr = np.deg2rad(ser_lon), np.deg2rad(ser_lat)

   dlat = ser_latr - clatr
   dlon = ser_lonr - clonr
   #dlon = (dlon + np.pi) % (2 * np.pi) - np.pi

   hav = np.sin(dlat / 2)**2 + np.cos(clatr) * np.cos(ser_latr) * np.sin(dlon / 2)**2
   hav = np.clip(hav, 0, 1)
   gcdr = 2 * np.arctan2(np.sqrt(hav), np.sqrt(1 - hav))
   return np.rad2deg(gcdr)

if __name__ == '__main__':
   main()
