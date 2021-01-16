#Python script to plot cross-sections
import LT.box as B
from LT.datafile import dfile
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.legend_handler import HandlerTuple

#Use latex commands (e.g. \textit ot \textbf)
rc('text', usetex=True)
#Set default font to times new roman
font = {'family' : 'Times New Roman',
        'weight' : 'normal',
        'size'   : 12
}
plt.rc('font', **font)


dir_name = "./bin_centering_corrections/all_thnq_Q2_4to5/"

fname_80 = dir_name + "pm80_laget_bc_corr.txt"
fname_580_set1 = dir_name + "pm580_laget_bc_corr_set1.txt"
fname_580_set2 = dir_name + "pm580_laget_bc_corr_set2.txt"

fname_750_set1 = dir_name + "pm750_laget_bc_corr_set1.txt"
fname_750_set2 = dir_name + "pm750_laget_bc_corr_set2.txt"
fname_750_set3 = dir_name + "pm750_laget_bc_corr_set3.txt"

kin_80 = dfile(fname_80)

kin_580_s1 = dfile(fname_580_set1)
kin_580_s2 = dfile(fname_580_set2)

kin_750_s1 = dfile(fname_750_set1)
kin_750_s2 = dfile(fname_750_set2)
kin_750_s3 = dfile(fname_750_set3)


#--Get 80 MeV/c setting--

#determine rel stats error
dataXsec_80         = np.array(kin_80['fsiRC_dataXsec_fsibc_corr'])    #data Xsec (corrected for rad. corr and bin-center corr. using  Laget FSI model
dataXsec_80_err     = np.array(kin_80['fsiRC_dataXsec_fsibc_corr_err'])
stats_80            = dataXsec_80_err / dataXsec_80                    #relative stats. error

#require stats. error <= 50 %
thnq_80             = np.array(kin_80['xb'])[stats_80<=0.50]
pm_avg_80           = np.array(kin_80['yb'])[stats_80<=0.50]  #GeV/c

dataXsec_80         = np.array(kin_80['fsiRC_dataXsec_fsibc_corr'])[stats_80<=0.50]   
dataXsec_80_err     = np.array(kin_80['fsiRC_dataXsec_fsibc_corr_err'])[stats_80<=0.50]

#make plots
#fig1 = B.pl.figure(figsize=(8,6))
#thnq0 = 35
#B.plot_exp(pm_avg_80[thnq_80==thnq0], dataXsec_80[thnq_80==thnq0], dataXsec_80_err[thnq_80==thnq0], marker='s', linestyle='none', color='b', mfc='b', ecolor='b', capsize=0, logy='true', label=r'80 MeV/c setting' )
#B.pl.show()


#--Get 580 MeV/c setting--

#determine rel stats error
dataXsec_580_s1         = np.array(kin_580_s1['fsiRC_dataXsec_fsibc_corr'])    
dataXsec_580_err_s1     = np.array(kin_580_s1['fsiRC_dataXsec_fsibc_corr_err'])
stats_580_s1            = dataXsec_580_err_s1 / dataXsec_580_s1

dataXsec_580_s2         = np.array(kin_580_s2['fsiRC_dataXsec_fsibc_corr'])    
dataXsec_580_err_s2     = np.array(kin_580_s2['fsiRC_dataXsec_fsibc_corr_err'])
stats_580_s2            = dataXsec_580_err_s2 / dataXsec_580_s2

#require stats. error <= 50 %
thnq_580_s1             = np.array(kin_580_s1['xb'])[stats_580_s1<=0.50]
pm_avg_580_s1           = np.array(kin_580_s1['yb'])[stats_580_s1<=0.50]  #GeV/c
dataXsec_580_s1         = np.array(kin_580_s1['fsiRC_dataXsec_fsibc_corr'])[stats_580_s1<=0.50]   
dataXsec_580_err_s1     = np.array(kin_580_s1['fsiRC_dataXsec_fsibc_corr_err'])[stats_580_s1<=0.50]

thnq_580_s2             = np.array(kin_580_s2['xb'])[stats_580_s2<=0.50]
pm_avg_580_s2           = np.array(kin_580_s2['yb'])[stats_580_s2<=0.50]  #GeV/c
dataXsec_580_s2         = np.array(kin_580_s2['fsiRC_dataXsec_fsibc_corr'])[stats_580_s2<=0.50]   
dataXsec_580_err_s2     = np.array(kin_580_s2['fsiRC_dataXsec_fsibc_corr_err'])[stats_580_s2<=0.50]


#--Get 750 MeV/c setting--

#determine rel stats error
dataXsec_750_s1         = np.array(kin_750_s1['fsiRC_dataXsec_fsibc_corr'])    
dataXsec_750_err_s1     = np.array(kin_750_s1['fsiRC_dataXsec_fsibc_corr_err'])
stats_750_s1            = dataXsec_750_err_s1 / dataXsec_750_s1

dataXsec_750_s2         = np.array(kin_750_s2['fsiRC_dataXsec_fsibc_corr'])    
dataXsec_750_err_s2     = np.array(kin_750_s2['fsiRC_dataXsec_fsibc_corr_err'])
stats_750_s2            = dataXsec_750_err_s2 / dataXsec_750_s2

dataXsec_750_s3         = np.array(kin_750_s3['fsiRC_dataXsec_fsibc_corr'])    
dataXsec_750_err_s3     = np.array(kin_750_s3['fsiRC_dataXsec_fsibc_corr_err'])
stats_750_s3            = dataXsec_750_err_s3 / dataXsec_750_s3


#require stats. error <= 50 %
thnq_750_s1             = np.array(kin_750_s1['xb'])[stats_750_s1<=0.50]
pm_avg_750_s1           = np.array(kin_750_s1['yb'])[stats_750_s1<=0.50]  #GeV/c
dataXsec_750_s1         = np.array(kin_750_s1['fsiRC_dataXsec_fsibc_corr'])[stats_750_s1<=0.50]   
dataXsec_750_err_s1     = np.array(kin_750_s1['fsiRC_dataXsec_fsibc_corr_err'])[stats_750_s1<=0.50]

thnq_750_s2             = np.array(kin_750_s2['xb'])[stats_750_s2<=0.50]
pm_avg_750_s2           = np.array(kin_750_s2['yb'])[stats_750_s2<=0.50]  #GeV/c
dataXsec_750_s2         = np.array(kin_750_s2['fsiRC_dataXsec_fsibc_corr'])[stats_750_s2<=0.50]   
dataXsec_750_err_s2     = np.array(kin_750_s2['fsiRC_dataXsec_fsibc_corr_err'])[stats_750_s2<=0.50]

thnq_750_s3             = np.array(kin_750_s3['xb'])[stats_750_s3<=0.50]
pm_avg_750_s3           = np.array(kin_750_s3['yb'])[stats_750_s3<=0.50]  #GeV/c
dataXsec_750_s3         = np.array(kin_750_s3['fsiRC_dataXsec_fsibc_corr'])[stats_750_s3<=0.50]   
dataXsec_750_err_s3     = np.array(kin_750_s3['fsiRC_dataXsec_fsibc_corr_err'])[stats_750_s3<=0.50]

#make plots for higher missing momentum
fig1 = B.pl.figure(figsize=(8,6))
thnq0 = 45

#B.plot_exp(pm_avg_80[thnq_80==thnq0], dataXsec_80[thnq_80==thnq0], dataXsec_80_err[thnq_80==thnq0], marker='D', linestyle='none', color='b', mfc='b', ecolor='b', capsize=0, logy='true', label=r'80 MeV/c setting' )

B.plot_exp(pm_avg_580_s1[thnq_580_s1==thnq0], dataXsec_580_s1[thnq_580_s1==thnq0], dataXsec_580_err_s1[thnq_580_s1==thnq0], marker='s', linestyle='none', color='k', mfc='white', ecolor='k', capsize=0, logy='true', label=r'580 MeV/c setting (set1)' )
B.plot_exp(pm_avg_580_s2[thnq_580_s2==thnq0]+0.001, dataXsec_580_s2[thnq_580_s2==thnq0], dataXsec_580_err_s2[thnq_580_s2==thnq0], marker='o', linestyle='none', color='r', mfc='white', ecolor='r', capsize=0, logy='true', label=r'580 MeV/c setting (set2)' )


B.plot_exp(pm_avg_750_s1[thnq_750_s1==thnq0]+0.002, dataXsec_750_s1[thnq_750_s1==thnq0], dataXsec_750_err_s1[thnq_750_s1==thnq0], marker='^', linestyle='none', color='g', mfc='white', ecolor='g', capsize=0, logy='true', label=r'750 MeV/c setting (set1)' )
B.plot_exp(pm_avg_750_s2[thnq_750_s2==thnq0]+0.003, dataXsec_750_s2[thnq_750_s2==thnq0], dataXsec_750_err_s2[thnq_750_s2==thnq0], marker='v', linestyle='none', color='b', mfc='white', ecolor='b', capsize=0, logy='true', label=r'750 MeV/c setting (set2)' )
#B.plot_exp(pm_avg_750_s3[thnq_750_s3==thnq0]+0.004, dataXsec_750_s3[thnq_750_s3==thnq0], dataXsec_750_err_s3[thnq_750_s3==thnq0], marker='<', linestyle='none', color='m', mfc='white', ecolor='m', capsize=0, logy='true', label=r'750 MeV/c setting (set3)' )


#Plot Labels/Titles
B.pl.title(r'Cross Section, $\theta_{nq}=45\pm5^{\circ}$', fontsize=20)
B.pl.xlabel('$p_{\mathrm{r}}$ [GeV/c]', fontsize=20)
B.pl.ylabel('$d^{5}\sigma/d\Omega_{e}d\Omega_{p}d\omega$ $\mathrm{[\mu b \cdot sr^{-2} \cdot MeV^{-1}]}$ ', fontsize=20)

plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
B.pl.legend(fontsize=18)

B.pl.show()



'''
stats_80            = np.array(kin_80['tot_stats_err'])  #relative stats error
thnq_80             = np.array(kin_80['thnq_bin'])[stats_80<=0.50]
pm_avg_80           = np.array(kin_80['pm_avg'])[stats_80<=0.50]  #GeV/c
red_Xsec_80         = np.array(kin_80['red_dataXsec_avg'])[stats_80<=0.50]  #avg. red Xsec [fm^3]
red_Xsec_sys_err_80 = np.array(kin_80['red_dataXsec_avg_syst_err'])[stats_80<=0.50]  #avg. red Xsec systematic error [fm^3]
red_Xsec_tot_err_80 = np.array(kin_80['red_dataXsec_avg_tot_err'])[stats_80<=0.50]  #avg. red Xsec tot error [fm^3]
rel_syst_err_80 = np.array(kin_80['tot_syst_err'])[stats_80<=0.50]
rel_tot_err_80 = np.array(kin_80['tot_err'])[stats_80<=0.50]


fig1 = B.pl.figure(figsize=(8,6))

#B.plot_exp(pm_avg_80[thnq_80==35], red_Xsec_80[thnq_80==35], red_Xsec_tot_err_80[thnq_80==35], marker='s', linestyle='none', color='b', mfc='b', ecolor='b', capsize=0, logy='true', label=r'80 MeV/c setting' )
#B.plot_exp(pm_avg_80[thnq_80==35], red_Xsec_80[thnq_80==35], red_Xsec_sys_err_80[thnq_80==35], marker='s', linestyle='none', color='b', mfc='b', ecolor='lightgray', elinewidth=3, capsize=0, logy='true')

B.plot_exp(pm_avg_580[thnq_580==35], red_Xsec_580[thnq_580==35], red_Xsec_tot_err_580[thnq_580==35], marker='o', linestyle='none', color='r', mfc='r', ecolor='r', capsize=0, logy='true', label=r'580 MeV/c setting' )
B.plot_exp(pm_avg_580[thnq_580==35], red_Xsec_580[thnq_580==35], red_Xsec_sys_err_580[thnq_580==35], marker='o', linestyle='none', color='r', mfc='r', ecolor='lightgray', elinewidth=3, capsize=0, logy='true')

B.plot_exp(pm_avg_750[thnq_750==35]+0.01, red_Xsec_750[thnq_750==35], red_Xsec_tot_err_750[thnq_750==35], marker='^', linestyle='none', color='magenta', mfc='magenta', ecolor='magenta', capsize=0, logy='true', label=r'750 MeV/c setting' )
B.plot_exp(pm_avg_750[thnq_750==35]+0.01, red_Xsec_750[thnq_750==35], red_Xsec_sys_err_750[thnq_750==35], marker='^', linestyle='none', color='magenta', mfc='magenta', ecolor='lightgray', elinewidth=3, capsize=0, logy='true')

plt.subplots_adjust(left=0.12, right=0.96, top=0.93, bottom=0.11)

#Set X-limits
B.pl.xlim(0.4, 1.2)
#Set Y-limits
B.pl.ylim(2e-6, 2e-4)

#Plot Labels/Titles
B.pl.title(r'Reduced Cross Section, $\theta_{nq}=35\pm5^{\circ}$', fontsize=20)
B.pl.xlabel('$p_{\mathrm{r}}$ [GeV/c]', fontsize=20)
B.pl.ylabel('$\sigma_{\mathrm{red}}$ [fm$^{3}$]', fontsize=20)

plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
B.pl.legend(fontsize=18)

#B.pl.show()
B.pl.savefig('redXsec_35.pdf')
#-----
'''
