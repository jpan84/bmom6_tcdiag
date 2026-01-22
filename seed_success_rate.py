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

lags = np.arange(6, 13, 6.)
rngs = np.arange(1, 3.1, 0.5)
#lags = np.arange(6, 7, 6.)
#rngs = np.arange(1, 3, 1.)
tab1 = np.zeros((lags.size, rngs.size)) #count exactly one match
tab2 = np.zeros(tab1.shape) #more than one match

def main():
   sdf = pd.read_csv(SEDFIL)
   pdf = pd.read_csv(PARFIL)
   #print(sdf.head())
   #print(pdf.head())

   sdf['dt'] = sdf['dt'].apply(cftime.datetime.strptime, args=(dt_fmt,), calendar='noleap')
   sdf = sdf[(sdf['clat'] >= LAT1) & (sdf['clat'] < LAT2) | (sdf['clat'] <= LAT1) & (sdf['clat'] > LAT2)]

   pdf = pdf[pdf['isgen']]
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
            print('\t', rw['dt'])
            dtmatch = pdf[(pdf['dt'] - rw['dt'] <= dt.timedelta(hours=ll)) & (pdf['dt'] - rw['dt'] >= dt.timedelta(hours=lags[0]))]
            #print(rw['dt'])
            #print(dtmatch)
            inrng = dtmatch[gcd_deg(rw['clon'], rw['clat'], dtmatch['lon'], dtmatch['lat']) <= rg]
            if inrng.shape[0] == 1:
               tab1[ii, jj] += 1
            elif inrng.shape[0] > 1:
               tab2[ii, jj] += 1

   print(tab1 / sdf.shape[0])
   print('\n')
   print(tab2 / sdf.shape[0])

def gcd_deg(clon, clat, ser_lon, ser_lat):
   clonr, clatr = np.deg2rad(clon), np.deg2rad(clat)
   ser_lonr, ser_latr = np.deg2rad(ser_lon), np.deg2rad(ser_lat)

   dlat = ser_latr - clatr
   dlon = ser_lonr - clonr

   hav = np.sin(dlat / 2)**2 + np.cos(clatr) * np.cos(ser_latr) * np.sin(dlon / 2)**2
   gcdr = 2 * np.arctan2(np.sqrt(hav), np.sqrt(1 - hav))
   return np.rad2deg(gcdr)

if __name__ == '__main__':
   main()
