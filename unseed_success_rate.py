import pandas as pd
import numpy as np
import sys
import cftime
import datetime as dt

SEDFIL = sys.argv[1]
PARFIL = sys.argv[2]
#LAT1 = float(sys.argv[3]) #lower bound (positive, inclusive)
#LAT2 = float(sys.argv[4]) #upper bound (positive, exclusive)
dt_fmt = '%Y-%m-%d %H:%M:%S'
print(sys.argv, '\n')

MAXLAG, MINLAG = 0, 0 #12, -12
RNG = 0.01 #6.0
#rngs = np.arange(2, 3.1, 0.5)
#lags = np.arange(6, 7, 6.)
#rngs = np.arange(1, 3, 1.)
#tab1 = np.zeros((lags.size, rngs.size)) #count exactly one match
#tab2 = np.zeros(tab1.shape) #more than one match
frac1, frac2 = 0.0, 0.0

def main():
   global frac1, frac2
   sdf = pd.read_csv(SEDFIL)
   pdf = pd.read_csv(PARFIL)
   #print(sdf.head())
   #print(pdf.head())

   sdf['dt'] = sdf['dt'].apply(cftime.datetime.strptime, args=(dt_fmt,), calendar='noleap')
   sdf = sdf[sdf['psmin'].notna()]
   print(sdf.head())
   #sdf = sdf[(sdf['clat'] >= LAT1) & (sdf['clat'] < LAT2) | (sdf['clat'] <= -LAT1) & (sdf['clat'] > -LAT2)]

   try:
      pdf['dt_obj'] = pdf['dt'].apply(cftime.datetime.strptime, args=(dt_fmt,), calendar='noleap')
   except ValueError: #tc_stats post-processing accidentally interpolated leap days
      pdf['dt_obj'] = pdf['dt'].apply(cftime.datetime.strptime, args=(dt_fmt,), calendar='Gregorian')
      #print(pdf['dt'].apply(lambda x: x.day))
      pdf = pdf[~pdf['dt_obj'].apply(lambda x: x.month == 2 and x.day == 29)]
      pdf['dt_obj'] = pdf['dt_obj'].apply(lambda d: cftime.datetime(d.year, d.month, d.day, d.hour, calendar='noleap'))
   pdf['dt_obj'] = pdf['dt_obj'] - dt.timedelta(days=2000 * 365)
   #pdf['istouched'] = False #flag for whether a discrete storm has been touched by unseed
   pdf['unseed_pt'] = False #flag for whether a fix is associated with an unseed event


   for ix, rw in sdf.iterrows():
      if rw['dt'].day == 1 and rw['dt'].hour == 0:
         print('\t', rw['dt'])
      dtmatch = pdf[(pdf['dt_obj'] - rw['dt'] <= dt.timedelta(hours=MAXLAG)) & (pdf['dt_obj'] - rw['dt'] >= dt.timedelta(hours=MINLAG))] #allow for SN trajectory gap
      #print(rw['dt'])
      #print(dtmatch)
      inrng = dtmatch[gcd_deg(rw['clon'], rw['clat'], dtmatch['lon'], dtmatch['lat']) <= RNG]
      #inrng = inrng.drop_duplicates(subset=['stmnum'], keep='first')
      if inrng.shape[0] == 1:
         frac1 += 1
      elif inrng.shape[0] > 1:
         print(inrng)
         frac2 += 1

      #flag storms and fixes that were touched by unseed
      #pdf.loc[pdf['stmnum'].isin(inrng['stmnum']), 'istouched'] = True
      pdf.loc[inrng.index, 'unseed_pt'] = True

   #count the number of unseed events applied to each trajectory
   attempt_cnt = pdf[pdf['unseed_pt'] == True].groupby('stmnum').size()
   pdf['istouched'] = pdf['stmnum'].map(attempt_cnt).fillna(0).astype(int)

   print(sdf.shape)
   print('Exactly one match:', frac1 / sdf.shape[0])
   print('\n')
   print('More than one match:', frac2 / sdf.shape[0], '\n')

   #tchd = pdf[pdf['istouched']].drop_duplicates(subset=['stmnum'], keep='first').shape[0]
   #nstm = pdf.drop_duplicates(subset=['stmnum'], keep='first').shape[0]
   #print('storms touched by unseed:', tchd, '/', nstm)

   pdf = pdf.drop(columns=['dt_obj'])
   pdf.to_csv(PARFIL + '.flagged.csv')
   pdf.to_parquet(PARFIL + '.flagged.parquet')

def gcd_deg(clon, clat, ser_lon, ser_lat):
   clonr, clatr = np.deg2rad(clon), np.deg2rad(clat)
   ser_lonr, ser_latr = np.deg2rad(ser_lon), np.deg2rad(ser_lat)

   dlat = ser_latr - clatr
   dlon = ser_lonr - clonr
   dlon = (dlon + np.pi) % (2 * np.pi) - np.pi

   hav = np.sin(dlat / 2)**2 + np.cos(clatr) * np.cos(ser_latr) * np.sin(dlon / 2)**2
   hav = np.clip(hav, 0, 1)
   gcdr = 2 * np.arctan2(np.sqrt(hav), np.sqrt(1 - hav))
   return np.rad2deg(gcdr)

if __name__ == '__main__':
   main()
