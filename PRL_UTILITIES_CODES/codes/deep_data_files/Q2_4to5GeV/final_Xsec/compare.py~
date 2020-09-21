import LT.box as B
from LT.datafile import dfile
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy.ma as ma
import sys                                     
import os                                                                                                       
from sys import argv  
import matplotlib
from matplotlib import rc
from matplotlib import *
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, InsetPosition, mark_inset)

MeV2fm = 197.3**3 
GeV2fm = 0.1973**3 
#1 fm = (1/197.3) MeV^-1


#READ WB pm_35 deg
#Conversion factor: 1 fm = 1 / (0.1973 GeV)  or  1 fm^-1 = 0.1973 GeV
fminv2GeV = 1./(1./0.1973)
GeV2fm = 0.1973**3  #inverse fermi^3 to GeV^3

fname_wb = "pm_distribution_results_q4_35.00.data"
kin = dfile(fname_wb)

p_miss = np.array(kin['p_miss'])/1000.   # in GeV
th_e = np.array(kin['th_e'])                #in deg
Ei = np.array(kin['Ei'])              #in MeV
omega = np.array(kin['omega'])        #in MeV
qlab = np.array(kin['qlab'])         #in MeV
pf = np.array(kin['pf'])             #in MeV
sig_red_exp = np.array(kin['sig_red_exp']) 
y = np.zeros(len(th_e))

#Get average kinematics
def get_avg_kinWB(header, ithnq):
    fname = 'av_kin_theta1/theta_pm_%i.00_fsi_norad_avgkin.data'%(ithnq)
    kin = dfile(fname)
    arr = (np.array(kin[header]))
    return arr

#-------READ CURRENT FILE-------
#Get average kinematics
def get_avg_kin(header, pm, data_set, ithnq ):
    fname = 'pm%i_fsi_norad_avgkin_set%i.txt'%(pm, data_set)
    kin = dfile(fname)
    thnq = np.array(kin['xb'])
    arr = (np.array(kin[header]))[thnq==ithnq]
    return arr

def get_redXsec(pm, data_set, ithnq ):
    fname = 'pm%i_laget_bc_corr_set%i.txt'%(pm, data_set)
    kin = dfile(fname)
    thnq = np.array(kin['xb'])
    arr = (np.array(kin['red_dataXsec']))[thnq==ithnq]* MeV2fm
    return arr  #in fm^3


#get_avg_kinWB(header, ithnq)


#R = (get_redXsec(80, 1, 35) - sig_red_exp[0])/get_redXsec(80, 1, 35) * 100.
#B.plot_exp(get_avg_kin('pm', 80, 1, 35)/1000., R, marker='^', color='r', label=r'Current 80 MeV $\sigma_{red}$ fm$^{3}$ ', logy=False)
#B.pl.show()

#B.plot_exp(p_miss, sig_red_exp, marker='o', color='black', label=r'ORIGINAL  $\sigma_{red}$, units$ ', logy=True)
#B.plot_exp(get_avg_kin('pm', 80, 1, 35)/1000.,  get_redXsec(80, 1, 35), marker='^', color='r', label=r'Current 80 MeV $\sigma_{red}$ fm$^{3}$ ', logy=True)
#B.plot_exp(get_avg_kin('pm', 580, 1, 35)/1000.,  get_redXsec(580, 1, 35), marker='^', color='r', label=r'Current 580 (set1) $\sigma_{red}$ fm$^{3}$ ', logy=True)
#B.plot_exp(get_avg_kin('pm', 580, 2, 35)/1000.,  get_redXsec(580, 2, 35), marker='^', color='r', label=r'Current 580 (set2) $\sigma_{red}$ fm$^{3}$ ', logy=True)
#B.plot_exp(get_avg_kin('pm', 750, 1, 35)/1000.,  get_redXsec(750, 1, 35), marker='^', color='r', label=r'Current 750 (set1) $\sigma_{red}$ fm$^{3}$ ', logy=True)
#B.plot_exp(get_avg_kin('pm', 750, 2, 35)/1000.,  get_redXsec(750, 2, 35), marker='^', color='r', label=r'Current 750 (set2) $\sigma_{red}$ fm$^{3}$ ', logy=True)
#B.plot_exp(get_avg_kin('pm', 750, 3, 35)/1000.,  get_redXsec(750, 3, 35), marker='^', color='r', label=r'Current 750 (set3) $\sigma_{red}$ fm$^{3}$ ', logy=True)

#B.pl.xlabel('Averaged Missing Momentum [GeV/c]')
#B.pl.ylabel('Reduced Cross Sections')
#B.pl.legend()
#B.pl.show()


#B.plot_exp(p_miss, th_e, marker='o', color='black', label=r'ORIGINAL$ ', logy=False)
B.plot_exp(get_avg_kinWB('pm', 75)/1000.,  get_avg_kinWB('th_e', 75), marker='o', color='k', label=r'Original: 75 deg ', logy=False)

B.plot_exp(get_avg_kin('pm', 80, 1, 75)/1000.,  get_avg_kin('th_e', 80, 1, 75), marker='^', color='r', label=r'Current 80 MeV ', logy=False)
B.plot_exp(get_avg_kin('pm', 580, 1, 75)/1000.,  get_avg_kin('th_e', 580, 1, 75), marker='^', color='r', label=r'Current 580 (set1) ', logy=False)
B.plot_exp(get_avg_kin('pm', 580, 2, 75)/1000.,  get_avg_kin('th_e', 580, 2, 75), marker='^', color='r', label=r'Current 580 (set2) ', logy=False)
B.plot_exp(get_avg_kin('pm', 750, 1, 75)/1000.,  get_avg_kin('th_e', 750, 1, 75), marker='^', color='r', label=r'Current 750 (set1) ', logy=False)
B.plot_exp(get_avg_kin('pm', 750, 2, 75)/1000.,  get_avg_kin('th_e', 750, 2, 75), marker='^', color='r', label=r'Current 750 (set2) ', logy=False)
B.plot_exp(get_avg_kin('pm', 750, 3, 75)/1000.,  get_avg_kin('th_e', 750, 3, 75), marker='^', color='r', label=r'Current 750 (set3) ', logy=False)


B.pl.xlabel('Averaged Missing Momentum [GeV/c]')
B.pl.ylabel('Averaged Electron Angle')
B.pl.legend()
B.pl.show()

#-------------

#B.plot_exp(p_miss, Ei, marker='o', color='black', label=r'ORIGINAL: 75 deg ', logy=False)
B.plot_exp(get_avg_kinWB('pm', 75)/1000.,  get_avg_kinWB('Ei', 75), marker='o', color='k', label=r'Original: 75 deg ', logy=False)

B.plot_exp(get_avg_kin('pm', 80, 1, 75)/1000.,  get_avg_kin('Ei', 80, 1, 75), marker='^', color='r', label=r'Current 80 MeV ', logy=False)
B.plot_exp(get_avg_kin('pm', 580, 1, 75)/1000.,  get_avg_kin('Ei', 580, 1, 75), marker='^', color='r', label=r'Current 580 (set1) ', logy=False)
B.plot_exp(get_avg_kin('pm', 580, 2, 75)/1000.,  get_avg_kin('Ei', 580, 2, 75), marker='^', color='r', label=r'Current 580 (set2) ', logy=False)
B.plot_exp(get_avg_kin('pm', 750, 1, 75)/1000.,  get_avg_kin('Ei', 750, 1, 75), marker='^', color='r', label=r'Current 750 (set1) ', logy=False)
B.plot_exp(get_avg_kin('pm', 750, 2, 75)/1000.,  get_avg_kin('Ei', 750, 2, 75), marker='^', color='r', label=r'Current 750 (set2) ', logy=False)
B.plot_exp(get_avg_kin('pm', 750, 3, 75)/1000.,  get_avg_kin('Ei', 750, 3, 75), marker='^', color='r', label=r'Current 750 (set3) ', logy=False)


B.pl.xlabel('Averaged Missing Momentum [GeV/c]')
B.pl.ylabel('Averaged Electron Beam Energy')
B.pl.legend()
B.pl.show()

#-------------
#B.plot_exp(p_miss, omega, marker='o', color='black', label=r'ORIGINAL$ ', logy=False)
B.plot_exp(get_avg_kinWB('pm', 75)/1000.,  get_avg_kinWB('omega', 75), marker='o', color='k', label=r'Original: 75 deg ', logy=False)

B.plot_exp(get_avg_kin('pm', 80, 1, 75)/1000.,  get_avg_kin('omega', 80, 1, 75), marker='^', color='r', label=r'Current 80 MeV ', logy=False)
B.plot_exp(get_avg_kin('pm', 580, 1, 75)/1000.,  get_avg_kin('omega', 580, 1, 75), marker='^', color='r', label=r'Current 580 (set1) ', logy=False)
B.plot_exp(get_avg_kin('pm', 580, 2, 75)/1000.,  get_avg_kin('omega', 580, 2, 75), marker='^', color='r', label=r'Current 580 (set2) ', logy=False)
B.plot_exp(get_avg_kin('pm', 750, 1, 75)/1000.,  get_avg_kin('omega', 750, 1, 75), marker='^', color='r', label=r'Current 750 (set1) ', logy=False)
B.plot_exp(get_avg_kin('pm', 750, 2, 75)/1000.,  get_avg_kin('omega', 750, 2, 75), marker='^', color='r', label=r'Current 750 (set2) ', logy=False)
B.plot_exp(get_avg_kin('pm', 750, 3, 75)/1000.,  get_avg_kin('omega', 750, 3, 75), marker='^', color='r', label=r'Current 750 (set3) ', logy=False)


B.pl.xlabel('Averaged Missing Momentum [GeV/c]')
B.pl.ylabel('Averaged Energy Transfer')
B.pl.legend()
B.pl.show()


#---------
#B.plot_exp(p_miss, qlab, marker='o', color='black', label=r'ORIGINAL$ ', logy=False)
B.plot_exp(get_avg_kinWB('pm', 75)/1000.,  get_avg_kinWB('q_lab', 75), marker='o', color='k', label=r'Original: 75 deg ', logy=False)

B.plot_exp(get_avg_kin('pm', 80, 1, 75)/1000.,  get_avg_kin('q_lab', 80, 1, 75), marker='^', color='r', label=r'Current 80 MeV ', logy=False)
B.plot_exp(get_avg_kin('pm', 580, 1, 75)/1000.,  get_avg_kin('q_lab', 580, 1, 75), marker='^', color='r', label=r'Current 580 (set1) ', logy=False)
B.plot_exp(get_avg_kin('pm', 580, 2, 75)/1000.,  get_avg_kin('q_lab', 580, 2, 75), marker='^', color='r', label=r'Current 580 (set2) ', logy=False)
B.plot_exp(get_avg_kin('pm', 750, 1, 75)/1000.,  get_avg_kin('q_lab', 750, 1, 75), marker='^', color='r', label=r'Current 750 (set1) ', logy=False)
B.plot_exp(get_avg_kin('pm', 750, 2, 75)/1000.,  get_avg_kin('q_lab', 750, 2, 75), marker='^', color='r', label=r'Current 750 (set2) ', logy=False)
B.plot_exp(get_avg_kin('pm', 750, 3, 75)/1000.,  get_avg_kin('q_lab', 750, 3, 75), marker='^', color='r', label=r'Current 750 (set3) ', logy=False)


B.pl.xlabel('Averaged Missing Momentum [GeV/c]')
B.pl.ylabel('Averaged 3-momentum Transfer')
B.pl.legend()
B.pl.show()

#---------
#B.plot_exp(p_miss, pf, marker='o', color='black', label=r'ORIGINAL$ ', logy=False)
B.plot_exp(get_avg_kinWB('pm', 75)/1000.,  get_avg_kinWB('pf', 75), marker='o', color='k', label=r'Original: 75 deg ', logy=False)

B.plot_exp(get_avg_kin('pm', 80, 1, 75)/1000.,  get_avg_kin('pf', 80, 1, 75), marker='^', color='r', label=r'Current 80 MeV ', logy=False)
B.plot_exp(get_avg_kin('pm', 580, 1, 75)/1000.,  get_avg_kin('pf', 580, 1, 75), marker='^', color='r', label=r'Current 580 (set1) ', logy=False)
B.plot_exp(get_avg_kin('pm', 580, 2, 75)/1000.,  get_avg_kin('pf', 580, 2, 75), marker='^', color='r', label=r'Current 580 (set2) ', logy=False)
B.plot_exp(get_avg_kin('pm', 750, 1, 75)/1000.,  get_avg_kin('pf', 750, 1, 75), marker='^', color='r', label=r'Current 750 (set1) ', logy=False)
B.plot_exp(get_avg_kin('pm', 750, 2, 75)/1000.,  get_avg_kin('pf', 750, 2, 75), marker='^', color='r', label=r'Current 750 (set2) ', logy=False)
B.plot_exp(get_avg_kin('pm', 750, 3, 75)/1000.,  get_avg_kin('pf', 750, 3, 75), marker='^', color='r', label=r'Current 750 (set3) ', logy=False)


B.pl.xlabel('Averaged Missing Momentum [GeV/c]')
B.pl.ylabel('Averaged Final Proton Momentum')
B.pl.legend()
B.pl.show()

