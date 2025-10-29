import numpy as np
import matplotlib.pyplot as plt

ROWS = ['JJA', 'SON']
COLS = ['CTRL$-$UNSEED', 'SEED$-$CTRL']
YLIMS = [(-2.01e8, 1.0e8), (-1.05e9, 5.0e8)]
LETTERS = [['(a)', '(b)'], ['(c)', '(d)']]

mP_TC = np.array([[1.89e8, -9.55e8], [1.75e8, -9.04e8]])
mP_zm = np.array([[4.66e7, -3.32e8], [3.17e7, -3.00e8]])

vqc_zm = np.array([[-3.53e7, 3.58e8], [-6.38e6, 4.53e8]])

E_TC = np.array([[-8.98e7, 3.69e8], [-8.14e7, 3.69e8]])
E_zm = np.array([[-3.02e7, -2.72e7], [-8.24e6, -1.57e8]])

#reverse the sign of UNSEED for easier comparison with SEED
BARS = [[('-P', -mP_zm), ('+E', -E_zm), ('$-\partial_y (vq)$', -vqc_zm)], [('-P', mP_zm), ('+E', E_zm), ('$-\partial_y (vq)$', vqc_zm)]]

fig, axes = plt.subplots(2, 2)

for ii, szn in enumerate(ROWS):
   for jj, cs in enumerate(COLS):
      axes[ii][jj].bar([bb[0] for bb in BARS[jj]], [bb[1][ii][jj] for bb in BARS[jj]])
      axes[ii][jj].bar('sum', sum([bb[1][ii][jj] for bb in BARS[jj]]), color='gray')
      print(sum([bb[1][ii][jj] for bb in BARS[jj]]))
      sgn = -1 if jj == 0 else 1
      axes[ii][jj].scatter(0, sgn * mP_TC[ii][jj], color='black')
      axes[ii][jj].scatter(1, sgn * E_TC[ii][jj], color='black')
      axes[ii][jj].axhline(0, color='gray', linewidth=0.5)
      axes[ii][jj].set_ylim(*YLIMS[jj])
      axes[ii][jj].set_title(LETTERS[ii][jj], loc='left')

      axes[ii][jj].set_ylabel(('$q$') + ' tendency [kg s$^{-1}$]') # if sgn == 1 else '$-q$'
      if ii == 0:
         axes[ii][jj].set_title(COLS[jj], fontsize=14)

#fig.suptitle('20-30°N water vapor budget', fontsize=16)
fig.suptitle('JJA SON')
fig.tight_layout()

plt.savefig('iwv_20-30_bar.png', bbox_inches='tight')
plt.show()
