import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

FILI = '/glade/u/home/jpan/aquaptc/tempest/260415_density_5exp/sznl_climo.csv' #from track_dens.py
COLS = ['mo_in_szn', 'varnm', 'case', 'NH', 'SH', 'totyrs']
nh_warm_mo = [6, 9]
sh_warm_mo = [2, 3]

DIRE = '~/aquaptc/tempest/'
EVNTS = [('250702_unseed_2hPa6m_unseed_events.parquet', 'us'), ('250415_unseed_production_unseed_events.parquet', 'us'), None, ('251229_seed_match_seed_events.parquet', 'sd'), ('250416_seed1x1_production_seed_events.parquet', 'sd')]
ORDER = ['unseed2', 'unseed', 'ctrl', 'mseed', 'seed']
PIELBL = [['(a)', '(b)'], ['(c)', '(d)']]

def main():
   climo = pd.read_csv(FILI, header=None, names=COLS).drop_duplicates()
   climo[['NH', 'SH']] = climo[['NH', 'SH']].div(climo['totyrs'], axis=0)
   print(climo)

   evdfs = [pd.read_parquet(os.path.join(DIRE, ev[0])) if ev is not None else None for ev in EVNTS]
   #TODO: filter unseed dfs by (~dp.isnan())
   evdfs = [df[~df['rp'].isna()] if df is not None else None for df in evdfs]
   print([df.shape if df is not None else None for df in evdfs])

   #print(tot_nTCs(climo))
   for ii, df in enumerate(evdfs):
      print('\nWorking on case', ORDER[ii])
      if EVNTS[ii] is None:
         continue
      selcase = climo[climo['case'] == ORDER[ii]]
      ev_pyr = df.shape[0] / selcase['totyrs'].iloc[0] #un/seed events per year
      tccnt = selcase[selcase['varnm'] == 'genesis points (storm count)']
      nTCs = (tccnt['NH'] + tccnt['SH']).sum()
      dn = nTCs - tot_nTCs(climo[climo['case'] == 'ctrl']).iloc[0]
      #print(dn)

      if EVNTS[ii][1] == 'sd':
         seed_gen = selcase[selcase['varnm'] == 'seeded genesis points']
         n_converted = (seed_gen['NH'] + seed_gen['SH']).sum()
         print('The seed conversion rate for %s is %.1f / %.1f = %.3f' % (ORDER[ii], n_converted, ev_pyr, n_converted / ev_pyr))
         print('The efficacy of seeds for %s is %.1f / %.1f = %.3f' % (ORDER[ii], dn, n_converted, dn / n_converted))

      if EVNTS[ii][1] == 'us':
         uscnts = [selcase[selcase['varnm'] == 'stms with %d unseeds' % nn] for nn in range(4)]
         nTCs_by_attempts = [(subdf['NH'] + subdf['SH']).sum() for subdf in uscnts]
         n_inelig = uscnts[0][~uscnts[0]['mo_in_szn'].isin(nh_warm_mo)]['NH'].sum() + uscnts[0][~uscnts[0]['mo_in_szn'].isin(sh_warm_mo)]['SH'].sum() #approx count of TCs ineligible for unseeding due to undergoing lysis in wrong season
         #print(n_inelig)
         missed = selcase[selcase['varnm'] == 'TCs missed by unseed']
         nTCs_by_attempts[0] = (missed['NH'] + missed['SH']).sum()
         print(nTCs_by_attempts)
         pct_of_TCs = 100. * np.array(nTCs_by_attempts) / sum(nTCs_by_attempts) #percent of TCs in the warm-season climo binned by number of attempts
         print(pct_of_TCs)
         n_admitted = np.dot(np.arange(4), nTCs_by_attempts) #number of unseed events spent on TCs that do show up in offline warm-season climo
         pct_admitted = 100. * n_admitted / ev_pyr
         print(n_admitted, pct_admitted)

         plt.rcParams['figure.figsize'] = (10, 4)
         plt.rc('font', size=13)
         fig, axes = plt.subplots(1, 2)

         #try0
         #axes[0].pie([n_admitted, ev_pyr - n_admitted], labels=['TCs', 'denied genesis'], autopct='%1.1f%%', explode=[0.1, 0])#labels=['Unseeding targeted at TCs', 'Unseeding preventing TCs']
         #axes[0].set_title('Unseeding events')

         #try0
         #axes[1].pie(np.arange(4) * nTCs_by_attempts, labels=np.arange(4), autopct='%1.1f%%') #labels=['Unseeding attempts consumed to remove TCs that took %d events' % nn for nn in range(4)],
         #axes[1].set_title('For unseeding events targeted at TCs,\nhow many events did it take to remove a TC?')

         #try0.5
         #axes[1].pie(nTCs_by_attempts, labels=np.arange(4), autopct='%1.1f%%') #labels=['Unseeding attempts consumed to remove TCs that took %d events' % nn for nn in range(4)],
         #axes[1].set_title('How many TCs experienced\n$n$ unseeding events?')

         ##try1
         #axes[0].pie([ev_pyr - n_admitted] + list(nTCs_by_attempts[1:] * np.arange(1, 4)), labels=['denied genesis'] + list(range(1, 4)), autopct=pie_pct_fmt(ev_pyr), explode=[0.1, 0, 0, 0])#labels=['Unseeding targeted at TCs', 'Unseeding preventing TCs']
         #axes[0].set_title('Unseeding interventions (%.1f annually)' % ev_pyr)

         ##try1
         #axes[1].pie(nTCs_by_attempts, labels=np.arange(4), autopct=pie_pct_fmt(nTCs - n_inelig)) #labels=['Unseeding attempts consumed to remove TCs that took %d events' % nn for nn in range(4)],
         #axes[1].set_title('Warm-season TCs (%.1f annually)\nby number of unseeding attempts' % (nTCs - n_inelig))

         #try1.5
         axes[0].pie([ev_pyr - n_admitted, np.dot(nTCs_by_attempts[1:], np.arange(1, 4))], labels=['denied\ngenesis', 'offline-tracked\nTCs'],
                      colors=['C8', 'C6'], autopct=pie_pct_fmt(ev_pyr), explode=[0, 0])#labels=['Unseeding targeted at TCs', 'Unseeding preventing TCs']
         axes[0].set_title('Unseeding interventions (%.1f annually)' % ev_pyr)
         axes[0].set_xlabel(PIELBL[ii][0])
         axes[0].set_ylabel(ORDER[ii], fontsize=16, labelpad=40)

         #try1.5
         axes[1].pie(nTCs_by_attempts, labels=np.arange(4), autopct=pie_pct_fmt(nTCs - n_inelig), colors=['C0', 'C2', 'C1', 'C3']) #labels=['Unseeding attempts consumed to remove TCs that took %d events' % nn for nn in range(4)],
         axes[1].set_title('Warm-season TCs (%.1f annually)\nby number of unseeding attempts' % (nTCs - n_inelig))
         axes[1].set_xlabel(PIELBL[ii][1])

         fig.tight_layout(w_pad=3)
         plt.savefig(os.path.join(os.path.dirname(FILI), '%s_unseed_pie.svg' % ORDER[ii]))
         plt.show()


#get total number of TCs by case
def tot_nTCs(df):
   tccnt = df[df['varnm'] == 'genesis points (storm count)']
   return (tccnt['NH'] + tccnt['SH']).groupby(df['case'], sort=False).sum()

def pie_pct_fmt(ntot):
   def autopct(pct):
       val = pct * ntot / 100.0
       return ('%1.1f%%\n(%.1f)' % (pct, val)) if pct > 10. else ''
   return autopct

if __name__ == '__main__':
   main()
