import LT.box as B
import numpy as np
from LT.datafile import dfile

f1 = dfile('redXsec_combined_original.txt')
f2 = dfile('redXsec_combined.txt')

thnq1 = f1['xb']
p_miss1 = (f1['pm_avg'])[thnq1==35]
norm_syst_tot1 = (f1['norm_syst_tot'])[thnq1==35]
y1 = 0. * p_miss1

thnq2 = f2['xb']
p_miss2 = (f2['pm_avg'])[thnq2==35]
norm_syst_tot2_35 = (f2['norm_syst_tot'])[thnq2==35]
norm_syst_tot2_45 = (f2['norm_syst_tot'])[thnq2==45]
norm_syst_tot2_75 = (f2['norm_syst_tot'])[thnq2==75]

y2 = 0. * p_miss2


print('norm_rel_err_35 = ',norm_syst_tot2_35)
print('norm_rel_err_45 = ',norm_syst_tot2_45)
print('norm_rel_err_75 = ',norm_syst_tot2_75)


#B.plot_exp(p_miss2, y2, norm_syst_tot2*100, marker='s', color='r', label='with tgt. wall / spec. acc systematics')
#B.plot_exp(p_miss1, y1, norm_syst_tot1*100, marker='o', color='b', label='original')

#B.pl.legend()
#B.pl.show('same')
