from LT.datafile import dfile
import LT.box as B
import numpy as np

f = dfile('pm580_laget_bc_corr_set1.txt')

th_nq = np.array(f['xb'])
pm = np.array(f['yb'])

a = np.array(f['pwiaRC_dataXsec'])
sig_a = np.array(f['pwiaRC_dataXsec_err'])
b = np.array(f['fsiRC_dataXsec'])
sig_b = np.array(f['fsiRC_dataXsec_err'])

    
R =  a / b

dR_da = 1./b
dR_db = -a/b**2

dR = np.sqrt( (dR_da)**2 * (sig_a)**2 + (dR_db)**2 * (sig_b)**2 )

B.plot_exp(pm, R, dR)
B.pl.ylim(-10, 10)
B.pl.show()
