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


fname_80 = "./red_dataXsec_pm80.txt"
fname_580 = "./red_dataXsec_pm580.txt"
fname_750 = "./red_dataXsec_pm750.txt"

kin_80 = dfile(fname_80)
kin_580 = dfile(fname_580)
kin_750 = dfile(fname_750)

#Get 80 MeV/c setting
stats_80            = np.array(kin_80['tot_stats_err'])  #relative stats error
thnq_80             = np.array(kin_80['thnq_bin'])[stats_80<=0.50]
pm_avg_80           = np.array(kin_80['pm_avg'])[stats_80<=0.50]  #GeV/c
red_Xsec_80         = np.array(kin_80['red_dataXsec_avg'])[stats_80<=0.50]  #avg. red Xsec [fm^3]
red_Xsec_sys_err_80 = np.array(kin_80['red_dataXsec_avg_syst_err'])[stats_80<=0.50]  #avg. red Xsec systematic error [fm^3]
red_Xsec_tot_err_80 = np.array(kin_80['red_dataXsec_avg_tot_err'])[stats_80<=0.50]  #avg. red Xsec tot error [fm^3]
rel_syst_err_80 = np.array(kin_80['tot_syst_err'])[stats_80<=0.50]
rel_tot_err_80 = np.array(kin_80['tot_err'])[stats_80<=0.50]

#Get 580 MeV/c setting
stats_580            = np.array(kin_580['tot_stats_err'])  #relative stats error
thnq_580             = np.array(kin_580['thnq_bin'])[stats_580<=0.50]
pm_avg_580           = np.array(kin_580['pm_avg'])[stats_580<=0.50]  #GeV/c
red_Xsec_580         = np.array(kin_580['red_dataXsec_avg'])[stats_580<=0.50]  #avg. red Xsec [fm^3]
red_Xsec_sys_err_580 = np.array(kin_580['red_dataXsec_avg_syst_err'])[stats_580<=0.50]  #avg. red Xsec systematic error [fm^3]
red_Xsec_tot_err_580 = np.array(kin_580['red_dataXsec_avg_tot_err'])[stats_580<=0.50]  #avg. red Xsec tot error [fm^3]
rel_syst_err_580 = np.array(kin_580['tot_syst_err'])[stats_580<=0.50]
rel_tot_err_580 = np.array(kin_580['tot_err'])[stats_580<=0.50]


#Get 750 MeV/c setting
stats_750            = np.array(kin_750['tot_stats_err'])  #relative stats error
thnq_750             = np.array(kin_750['thnq_bin'])[stats_750<=0.50]
pm_avg_750           = np.array(kin_750['pm_avg'])[stats_750<=0.50]  #GeV/c
red_Xsec_750         = np.array(kin_750['red_dataXsec_avg'])[stats_750<=0.50]  #avg. red Xsec [fm^3]
red_Xsec_sys_err_750 = np.array(kin_750['red_dataXsec_avg_syst_err'])[stats_750<=0.50]  #avg. red Xsec systematic error [fm^3]
red_Xsec_tot_err_750 = np.array(kin_750['red_dataXsec_avg_tot_err'])[stats_750<=0.50]  #avg. red Xsec tot error [fm^3]
rel_syst_err_750 = np.array(kin_750['tot_syst_err'])[stats_750<=0.50]
rel_tot_err_750 = np.array(kin_750['tot_err'])[stats_750<=0.50]

'''
fig1, (ax1, ax2) = B.pl.subplots(2, sharex=True)

ax1.errorbar(pm_avg_80[thnq_80==35], red_Xsec_80[thnq_80==35], red_Xsec_tot_err_80[thnq_80==35], marker='s', linestyle='none', color='b', mfc='b', ecolor='b', capsize=0, label=r'80 MeV/c setting' )
ax1.errorbar(pm_avg_80[thnq_80==35], red_Xsec_80[thnq_80==35], red_Xsec_sys_err_80[thnq_80==35], marker='s', linestyle='none', color='b', mfc='b', ecolor='lightgray', elinewidth=3, capsize=0)

ax1.errorbar(pm_avg_580[thnq_580==35], red_Xsec_580[thnq_580==35], red_Xsec_tot_err_580[thnq_580==35], marker='o', linestyle='none', color='r', mfc='r', ecolor='r', capsize=0, label=r'580 MeV/c setting' )
ax1.errorbar(pm_avg_580[thnq_580==35], red_Xsec_580[thnq_580==35], red_Xsec_sys_err_580[thnq_580==35], marker='o', linestyle='none', color='r', mfc='r', ecolor='lightgray', elinewidth=3, capsize=0)

ax1.errorbar(pm_avg_750[thnq_750==35], red_Xsec_750[thnq_750==35], red_Xsec_tot_err_750[thnq_750==35], marker='^', linestyle='none', color='magenta', mfc='magenta', ecolor='magenta', capsize=0, label=r'750 MeV/c setting' )
ax1.errorbar(pm_avg_750[thnq_750==35], red_Xsec_750[thnq_750==35], red_Xsec_sys_err_750[thnq_750==35], marker='^', linestyle='none', color='magenta', mfc='magenta', ecolor='lightgray', elinewidth=3, capsize=0)

ax1.set_yscale('log')

#Plot Labels/Titles
ax1.set_title(r'Reduced Cross Section, $\theta_{nq}=35\pm5^{\circ}$', fontsize=20)
ax2.set_xlabel('$p_{\mathrm{r}}$ [GeV/c]', fontsize=18)
ax1.set_ylabel('$\sigma_{\mathrm{red}}$ [fm$^{3}$]', fontsize=18)
ax1.legend()

plt.xticks(fontsize=16)
plt.yticks(fontsize=16)

#---Plot relative error bars
y_dummy = np.zeros(len(red_Xsec_80[thnq_80==35]))
l1, = ax2.plot([], [], c='b')
l2, = ax2.plot([], [], c='r')
l3, = ax2.plot([], [], c='magenta')
l4, = ax2.plot([], [], c='gray', label='Systematic Uncerainty')

ax2.errorbar(pm_avg_80[thnq_80==35], y_dummy, rel_tot_err_80[thnq_80==35], marker='s', linestyle='none', color='b', mfc='white', ecolor='b', capsize=0)
ax2.errorbar(pm_avg_80[thnq_80==35], y_dummy, rel_syst_err_80[thnq_80==35], marker='s', linestyle='none', color='b', mfc='white', ecolor='gray', elinewidth=3, capsize=0)

y_dummy = np.zeros(len(red_Xsec_580[thnq_580==35]))
ax2.errorbar(pm_avg_580[thnq_580==35], y_dummy, rel_tot_err_580[thnq_580==35], marker='o', linestyle='none', color='r', mfc='white', ecolor='r', capsize=0)
ax2.errorbar(pm_avg_580[thnq_580==35], y_dummy, rel_syst_err_580[thnq_580==35], marker='o', linestyle='none', color='r', mfc='white', ecolor='gray', elinewidth=3, capsize=0)

y_dummy = np.zeros(len(red_Xsec_750[thnq_750==35]))
ax2.errorbar(pm_avg_750[thnq_750==35], y_dummy, rel_tot_err_750[thnq_750==35], marker='^', linestyle='none', color='magenta', mfc='white', ecolor='magenta', capsize=0)
ax2.errorbar(pm_avg_750[thnq_750==35], y_dummy, rel_syst_err_750[thnq_750==35], marker='^', linestyle='none', color='magenta', mfc='white', ecolor='gray', elinewidth=3, capsize=0)

ax2.set_ylim([-0.7,0.7])

#Plot Labels/Titles
ax2.set_ylabel('Relative Uncertainty', fontsize=18)

#ax2.legend()
ax2.legend([(l1,l2, l3), l4], ['Total Uncertainty', 'Systematic Uncertainty'], handler_map={tuple: HandlerTuple(ndivide=None)}, loc='upper left')
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)

B.pl.show()
'''


fig1 = B.pl.figure(figsize=(8,6))

#B.plot_exp(pm_avg_80[thnq_80==35], red_Xsec_80[thnq_80==35], red_Xsec_tot_err_80[thnq_80==35], marker='s', linestyle='none', color='b', mfc='b', ecolor='b', capsize=0, logy='true', label=r'80 MeV/c setting' )
#B.plot_exp(pm_avg_80[thnq_80==35], red_Xsec_80[thnq_80==35], red_Xsec_sys_err_80[thnq_80==35], marker='s', linestyle='none', color='b', mfc='b', ecolor='lightgray', elinewidth=3, capsize=0, logy='true')

B.plot_exp(pm_avg_580[thnq_580==35], red_Xsec_580[thnq_580==35], red_Xsec_tot_err_580[thnq_580==35], marker='o', linestyle='none', color='r', mfc='r', ecolor='r', capsize=0, logy='true', label=r'580 MeV/c setting' )
B.plot_exp(pm_avg_580[thnq_580==35], red_Xsec_580[thnq_580==35], red_Xsec_sys_err_580[thnq_580==35], marker='o', linestyle='none', color='r', mfc='r', ecolor='lightgray', elinewidth=3, capsize=0, logy='true')

B.plot_exp(pm_avg_750[thnq_750==35], red_Xsec_750[thnq_750==35], red_Xsec_tot_err_750[thnq_750==35], marker='^', linestyle='none', color='magenta', mfc='magenta', ecolor='magenta', capsize=0, logy='true', label=r'750 MeV/c setting' )
B.plot_exp(pm_avg_750[thnq_750==35], red_Xsec_750[thnq_750==35], red_Xsec_sys_err_750[thnq_750==35], marker='^', linestyle='none', color='magenta', mfc='magenta', ecolor='lightgray', elinewidth=3, capsize=0, logy='true')

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
B.pl.legend(fontsize=16)

#B.pl.show()
B.pl.savefig('redXsec_35.pdf')
#-----

fig2 = B.pl.figure(figsize=(8,6))

#B.plot_exp(pm_avg_80[thnq_80==45], red_Xsec_80[thnq_80==45], red_Xsec_tot_err_80[thnq_80==45], marker='s', color='b', mfc='white', ecolor='b', capsize=0, logy='true', label=r'80 MeV/c setting' )
#B.plot_exp(pm_avg_80[thnq_80==45], red_Xsec_80[thnq_80==45], red_Xsec_sys_err_80[thnq_80==45], marker='s', color='b', mfc='white', ecolor='lightgray', elinewidth=3, capsize=0, logy='true')

B.plot_exp(pm_avg_580[thnq_580==45], red_Xsec_580[thnq_580==45], red_Xsec_tot_err_580[thnq_580==45], marker='o', color='r', ecolor='r', capsize=0, logy='true', label=r'580 MeV/c setting' )
B.plot_exp(pm_avg_580[thnq_580==45], red_Xsec_580[thnq_580==45], red_Xsec_sys_err_580[thnq_580==45], marker='o', color='r', ecolor='lightgray', elinewidth=3, capsize=0, logy='true')

B.plot_exp(pm_avg_750[thnq_750==45], red_Xsec_750[thnq_750==45], red_Xsec_tot_err_750[thnq_750==45], marker='^', color='magenta', ecolor='magenta', capsize=0, logy='true', label=r'750 MeV/c setting' )
B.plot_exp(pm_avg_750[thnq_750==45], red_Xsec_750[thnq_750==45], red_Xsec_sys_err_750[thnq_750==45], marker='^', color='magenta', ecolor='lightgray', elinewidth=3, capsize=0, logy='true')

plt.subplots_adjust(left=0.12, right=0.96, top=0.93, bottom=0.11)

#Set X-limits
B.pl.xlim(0.4, 1.2)
#Set Y-limits
B.pl.ylim(2e-6, 2e-4)

#Plot Labels/Titles
B.pl.title(r'Reduced Cross Section, $\theta_{nq}=45\pm5^{\circ}$', fontsize=20)
B.pl.xlabel('$p_{\mathrm{r}}$ [GeV/c]', fontsize=20)
B.pl.ylabel('$\sigma_{\mathrm{red}}$ [fm$^{3}$]', fontsize=20)

plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
B.pl.legend(fontsize=16)

#B.pl.show()
B.pl.savefig('redXsec_45.pdf')

#------

#-----

fig3 = B.pl.figure()

B.plot_exp(pm_avg_80[thnq_80==75], red_Xsec_80[thnq_80==75], red_Xsec_tot_err_80[thnq_80==75], marker='s', color='b', mfc='white', ecolor='b', capsize=0, logy='true', label=r'80 MeV/c setting' )
B.plot_exp(pm_avg_80[thnq_80==75], red_Xsec_80[thnq_80==75], red_Xsec_sys_err_80[thnq_80==75], marker='s', color='b', mfc='white', ecolor='lightgray', elinewidth=3, capsize=0, logy='true')

#B.plot_exp(pm_avg_580[thnq_580==75], red_Xsec_580[thnq_580==75], red_Xsec_tot_err_580[thnq_580==75], marker='o', color='r', ecolor='r', capsize=0, logy='true', label=r'580 MeV/c setting' )
#B.plot_exp(pm_avg_580[thnq_580==75], red_Xsec_580[thnq_580==75], red_Xsec_sys_err_580[thnq_580==75], marker='o', color='r', ecolor='lightgray', elinewidth=3, capsize=0, logy='true')

#B.plot_exp(pm_avg_750[thnq_750==75], red_Xsec_750[thnq_750==75], red_Xsec_tot_err_750[thnq_750==75], marker='^', color='magenta', ecolor='magenta', capsize=0, logy='true', label=r'750 MeV/c setting' )
#B.plot_exp(pm_avg_750[thnq_750==75], red_Xsec_750[thnq_750==75], red_Xsec_sys_err_750[thnq_750==75], marker='^', color='magenta', ecolor='lightgray', elinewidth=3, capsize=0, logy='true')

#Plot Labels/Titles
B.pl.title(r'Reduced Cross Section, $\theta_{nq}=75\pm5^{\circ}$', fontsize=20)
B.pl.xlabel('$p_{r}$ [GeV/c]', fontsize=18)
B.pl.ylabel('$\sigma_{red}$ [fm$^{3}$]', fontsize=18)

plt.xticks(fontsize=16)
plt.yticks(fontsize=16)

B.pl.legend()
B.pl.show()
