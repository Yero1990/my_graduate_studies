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
#from prl_utilities import *
from prl_utilities_test import *


#Use latex commands (e.g. \textit ot \textbf)
rc('text', usetex=True)
#Set default font to times new roman
font = {'family' : 'Times New Roman',
        'weight' : 'normal',
        'size'   : 12
}
plt.rc('font', **font)

#Set font
csfont = {'fontname':'Times New Roman'}

def make_prl_plots(plot2_inset=0):
  
    #Read Hall C (E12-10-003) experimental data (requires <=50% statistical uncertainty, units: fm^3, GeV/c)
    pm_avg35, red_dataXsec_avg_masked35, red_dataXsec_avg_tot_err35 = read_hallc_data(35)
    pm_avg45, red_dataXsec_avg_masked45, red_dataXsec_avg_tot_err45 = read_hallc_data(45)
    pm_avg75, red_dataXsec_avg_masked75, red_dataXsec_avg_tot_err75 = read_hallc_data(75)

    #Read Hall A experimental data
    pm_ha35, red_dataXsec_ha35, red_dataXsec_err_ha35 = read_halla_data(35)
    pm_ha45, red_dataXsec_ha45, red_dataXsec_err_ha45 = read_halla_data(45)  
    pm_ha75, red_dataXsec_ha75, red_dataXsec_err_ha75 = read_halla_data(75) 

    #Read Theoretical calculations
    #PWIA (or PWBA, which assumes the neutron, in addition to a proton can be a plane wave - free particle)
    pm_v18_35_pwia,  f_red_pwiaXsec_V18_35  = read_theoretical_models("V18", "PWIA", 35)
    pm_cd_35_pwia,   f_red_pwiaXsec_CD_35  = read_theoretical_models("CD-Bonn", "PWIA", 35)
    pm_jml_35_pwia,  f_red_pwiaXsec_JML_35  = read_theoretical_models("JML", "PWIA", 35)
    pm_wjc2_35_pwba, f_red_pwbaXsec_WJC2_35 = read_JWVO_theory(35, "WJC2", "GKex05", "PWBA")
    
    pm_v18_45_pwia,  f_red_pwiaXsec_V18_45  = read_theoretical_models("V18", "PWIA", 45)
    pm_cd_45_pwia,   f_red_pwiaXsec_CD_45   = read_theoretical_models("CD-Bonn", "PWIA", 45)
    pm_jml_45_pwia,  f_red_pwiaXsec_JML_45  = read_theoretical_models("JML", "PWIA", 45)
    pm_wjc2_45_pwba, f_red_pwbaXsec_WJC2_45 = read_JWVO_theory(45, "WJC2", "GKex05", "PWBA")

    pm_v18_75_pwia,  f_red_pwiaXsec_V18_75  = read_theoretical_models("V18", "PWIA", 75)
    pm_cd_75_pwia,   f_red_pwiaXsec_CD_75   = read_theoretical_models("CD-Bonn", "PWIA", 75)
    pm_jml_75_pwia,  f_red_pwiaXsec_JML_75  = read_theoretical_models("JML", "PWIA", 75)
    pm_wjc2_75_pwba, f_red_pwbaXsec_WJC2_75 = read_JWVO_theory(75, "WJC2", "GKex05", "PWBA")

    #FSI (or DWBA)
    pm_v18_35_fsi, f_red_fsiXsec_V18_35 = read_theoretical_models("V18", "FSI", 35)
    pm_cd_35_fsi,  f_red_fsiXsec_CD_35  = read_theoretical_models("CD-Bonn", "FSI", 35)
    pm_jml_35_fsi, f_red_fsiXsec_JML_35 = read_theoretical_models("JML", "FSI", 35)
    pm_wjc2_35_fsi, f_red_fsiXsec_WJC2_35 = read_JWVO_theory(35, "WJC2", "GKex05", "DWBA")

    pm_v18_45_fsi, f_red_fsiXsec_V18_45 = read_theoretical_models("V18", "FSI", 45)
    pm_cd_45_fsi,  f_red_fsiXsec_CD_45  = read_theoretical_models("CD-Bonn", "FSI", 45)
    pm_jml_45_fsi, f_red_fsiXsec_JML_45 = read_theoretical_models("JML", "FSI", 45)
    pm_wjc2_45_fsi, f_red_fsiXsec_WJC2_45 = read_JWVO_theory(45, "WJC2", "GKex05", "DWBA")

    pm_v18_75_fsi, f_red_fsiXsec_V18_75 = read_theoretical_models("V18", "FSI", 75)
    pm_cd_75_fsi,  f_red_fsiXsec_CD_75  = read_theoretical_models("CD-Bonn", "FSI", 75)
    pm_jml_75_fsi, f_red_fsiXsec_JML_75 = read_theoretical_models("JML", "FSI", 75)
    pm_wjc2_75_fsi, f_red_fsiXsec_WJC2_75 = read_JWVO_theory(75, "WJC2", "GKex05", "DWBA")

    print('=========MSV18 (35 deg)=========')
    print('pr(pwia) = ',  pm_v18_35_pwia)
    print('redXsec_pwia = ', f_red_pwiaXsec_V18_35( pm_v18_35_pwia) )
    print('pr(fsi) = ',  pm_v18_35_fsi)
    print('redXsec_fsi = ', f_red_fsiXsec_V18_35( pm_v18_35_fsi) )
    print('=========MSV18 (35 deg)=========')

    #======== Define Ratios of Reduced Cross Sections (to be used in PRL PLOT 2)=====
    
    #Data (and all models) ratio to CD-Bonn PWIA

    #=====================
    #= THETA_NQ = 35 DEG =
    #=====================

    R_data35 = red_dataXsec_avg_masked35 / f_red_pwiaXsec_CD_35(pm_avg35)
    R_data_err35 = red_dataXsec_avg_tot_err35 / f_red_pwiaXsec_CD_35(pm_avg35)

    R_HAdata35 = red_dataXsec_ha35 / f_red_pwiaXsec_CD_35(pm_ha35)
    R_HAdata_err35 = red_dataXsec_err_ha35 / f_red_pwiaXsec_CD_35(pm_ha35)

    R_CD_fsi35    = f_red_fsiXsec_CD_35(pm_avg35)    / f_red_pwiaXsec_CD_35(pm_avg35)
    R_JML_pwia35  = f_red_pwiaXsec_JML_35(pm_avg35)  / f_red_pwiaXsec_CD_35(pm_avg35)
    R_JML_fsi35   = f_red_fsiXsec_JML_35(pm_avg35)   / f_red_pwiaXsec_CD_35(pm_avg35)
    R_V18_pwia35  = f_red_pwiaXsec_V18_35(pm_avg35)  / f_red_pwiaXsec_CD_35(pm_avg35)
    R_V18_fsi35   = f_red_fsiXsec_V18_35(pm_avg35)   / f_red_pwiaXsec_CD_35(pm_avg35)
    R_WJC2_pwba35 = f_red_pwbaXsec_WJC2_35(pm_avg35) / f_red_pwiaXsec_CD_35(pm_avg35)
    R_WJC2_fsi35  = f_red_fsiXsec_WJC2_35(pm_avg35)  / f_red_pwiaXsec_CD_35(pm_avg35)


    #=====================
    #= THETA_NQ = 45 DEG =
    #=====================

    R_data45 = red_dataXsec_avg_masked45 / f_red_pwiaXsec_CD_45(pm_avg45)
    R_data_err45 = red_dataXsec_avg_tot_err45 / f_red_pwiaXsec_CD_45(pm_avg45)

    R_HAdata45 = red_dataXsec_ha45 / f_red_pwiaXsec_CD_45(pm_ha45)
    R_HAdata_err45 = red_dataXsec_err_ha45 / f_red_pwiaXsec_CD_45(pm_ha45)

    R_CD_fsi45    = f_red_fsiXsec_CD_45(pm_avg45)    / f_red_pwiaXsec_CD_45(pm_avg45)
    R_JML_pwia45  = f_red_pwiaXsec_JML_45(pm_avg45)  / f_red_pwiaXsec_CD_45(pm_avg45)
    R_JML_fsi45   = f_red_fsiXsec_JML_45(pm_avg45)   / f_red_pwiaXsec_CD_45(pm_avg45)
    R_V18_pwia45  = f_red_pwiaXsec_V18_45(pm_avg45)  / f_red_pwiaXsec_CD_45(pm_avg45)
    R_V18_fsi45   = f_red_fsiXsec_V18_45(pm_avg45)   / f_red_pwiaXsec_CD_45(pm_avg45)
    R_WJC2_pwba45 = f_red_pwbaXsec_WJC2_45(pm_avg45) / f_red_pwiaXsec_CD_45(pm_avg45)
    R_WJC2_fsi45  = f_red_fsiXsec_WJC2_45(pm_avg45)  / f_red_pwiaXsec_CD_45(pm_avg45)

    
    #=====================
    #= THETA_NQ = 75 DEG =
    #=====================

    R_data75 = red_dataXsec_avg_masked75 / f_red_pwiaXsec_CD_75(pm_avg75)
    R_data_err75 = red_dataXsec_avg_tot_err75 / f_red_pwiaXsec_CD_75(pm_avg75)

    R_HAdata75 = red_dataXsec_ha75 / f_red_pwiaXsec_CD_75(pm_ha75)
    R_HAdata_err75 = red_dataXsec_err_ha75 / f_red_pwiaXsec_CD_75(pm_ha75)

    R_CD_fsi75    = f_red_fsiXsec_CD_75(pm_avg75)    / f_red_pwiaXsec_CD_75(pm_avg75)
    R_JML_pwia75  = f_red_pwiaXsec_JML_75(pm_avg75)  / f_red_pwiaXsec_CD_75(pm_avg75)
    R_JML_fsi75   = f_red_fsiXsec_JML_75(pm_avg75)   / f_red_pwiaXsec_CD_75(pm_avg75)
    R_V18_pwia75  = f_red_pwiaXsec_V18_75(pm_avg75)  / f_red_pwiaXsec_CD_75(pm_avg75)
    R_V18_fsi75   = f_red_fsiXsec_V18_75(pm_avg75)   / f_red_pwiaXsec_CD_75(pm_avg75)
    R_WJC2_pwba75 = f_red_pwbaXsec_WJC2_75(pm_avg75) / f_red_pwiaXsec_CD_75(pm_avg75)
    R_WJC2_fsi75  = f_red_fsiXsec_WJC2_75(pm_avg75)  / f_red_pwiaXsec_CD_75(pm_avg75)

    print(">>>>>>>>>>>> ", red_dataXsec_avg_tot_err35)
    #------------------------------------------------------------------------------------
    #--------MAKE PRL PLOT 1 (Reduced Cross Sections vs. recoil momenta)-----------------
    #------------------------------------------------------------------------------------
    font_size = 16
    label_size = 16
    axes_fontsize = 19
    #-----Create Subplots-----
    fig = B.pl.subplots(3, sharex=True, sharey=True, figsize=(14.2, 3.9))   #figsize (width, height), (13.4,6) is in inches;  PRL requirements: width = 8.6 cm (3.38583 in) --> 1 manuscript column (for wider figs, use 1.5 or 2 columns: (17.2 cm or 6.771654 in)
    gs = gridspec.GridSpec(1, 3) 

    #=====================
    #= THETA_NQ = 35 DEG =
    #====================
    ax0 = B.pl.subplot(gs[0])

    B.pl.text(0.25, 0.5e-6, r'$\theta_{nq}=35\pm5^{o}$', fontsize=font_size)
    B.pl.text(1.0, 1e0, r'(a)', fontsize=font_size)

    #Plot Experimental Data (Hall A or Hall C)
    l1 = B.plot_exp(pm_avg35, red_dataXsec_avg_masked35, red_dataXsec_avg_tot_err35, marker='o', markersize=5, color='k', capsize=0, markerfacecolor='k', logy=True, label='This Experiment (Hall C)', zorder=4)
    l2 = B.plot_exp(pm_ha35, red_dataXsec_ha35, red_dataXsec_err_ha35, marker='s',  markersize=5, color='#ff1000', markerfacecolor='white', capsize=0, logy=True,  label='Hall A Data', zorder=3)

    #Plot theoretical calculations
    l3, = B.plot_exp(pm_jml_35_pwia, f_red_pwiaXsec_JML_35(pm_jml_35_pwia), linestyle='--', marker='None', color='#0000ff', logy=True, label='Paris PWIA', zorder=2)
    l4, = B.plot_exp(pm_jml_35_fsi,  f_red_fsiXsec_JML_35(pm_jml_35_fsi), linestyle='-', marker='None', color='#0000ff', logy=True, label='Paris FSI', zorder=2)
    
    l5, = B.plot_exp(pm_v18_35_pwia, f_red_pwiaXsec_V18_35(pm_v18_35_pwia), linestyle='--', marker='None', color='#009000', logy=True, label='AV18 PWIA', zorder=2)   
    l6, = B.plot_exp(pm_v18_35_fsi,  f_red_fsiXsec_V18_35(pm_v18_35_fsi), linestyle='-', marker='None', color='#009000', logy=True, label='AV18 FSI', zorder=2) 
    
    l7, = B.plot_exp(pm_cd_35_pwia, f_red_pwiaXsec_CD_35(pm_cd_35_pwia), linestyle='--', marker='None', color='#ff00ff', logy=True, label='CD-Bonn PWIA', zorder=2)     
    l8, = B.plot_exp(pm_cd_35_fsi,  f_red_fsiXsec_CD_35(pm_cd_35_fsi), linestyle='-', marker='None', color='#ff00ff', logy=True, label='CD-Bonn FSI', zorder=2) 

    l9, = B.plot_exp(pm_wjc2_35_pwba, f_red_pwbaXsec_WJC2_35(pm_wjc2_35_pwba), linestyle='--', marker='None', color='orange', logy=True, label='WJC2 (GKex05) PWBA', zorder=0 )
    l10, = B.plot_exp(pm_wjc2_35_fsi, f_red_fsiXsec_WJC2_35(pm_wjc2_35_fsi), linestyle='-', marker='None', color='orange', logy=True, label='WJC2 (GKex05) DWBA', zorder=0 )
    
    #Set axis limits
    B.pl.xlim(0.0, 1.19)
    B.pl.ylim(1e-7, 10.)
 
    #Set Tick Marks
    ax0.tick_params(axis='both', which='both', direction='in', top=True, bottom=True, left=True, right=True, labelsize=label_size)

    #Set Axes Labels for subplot 0
    B.pl.xlabel('')                                                                                                                                                                                                   
    B.pl.ylabel(r'$\sigma_{\mathrm{red}}$ ($\mathrm{fm}^{3}$) ', fontsize=axes_fontsize,  labelpad=10)                                                                                           
    B.pl.title('')

    #=====================
    #= THETA_NQ = 45 DEG =
    #====================
    ax1 = B.pl.subplot(gs[1], sharey=ax0)

    B.pl.text(0.25, 0.5e-6, r'$\theta_{nq}=45\pm5^{o}$', fontsize=font_size)
    B.pl.text(1.0, 1e0, r'(b)', fontsize=font_size)

    #Remove un-necessary Y tick labels from subplot
    B.pl.setp(ax1.get_yticklabels(), visible=False)
    
    #Plot Experimental Data (Hall A or Hall C)
    l1 = B.plot_exp(pm_avg45, red_dataXsec_avg_masked45, red_dataXsec_avg_tot_err45, marker='o', markersize=5, color='k', capsize=0, markerfacecolor='k', logy=True, label='This Experiment (Hall C)', zorder=4)
    l2 = B.plot_exp(pm_ha45, red_dataXsec_ha45, red_dataXsec_err_ha45, marker='s',  markersize=5, color='#ff1000', markerfacecolor='white', capsize=0, logy=True,  label='Hall A Data', zorder=3)

    #Plot theoretical calculations
    l3, = B.plot_exp(pm_jml_45_pwia, f_red_pwiaXsec_JML_45(pm_jml_45_pwia), linestyle='--', marker='None', color='#0000ff', logy=True, label='Paris PWIA', zorder=2)
    l4, = B.plot_exp(pm_jml_45_fsi,  f_red_fsiXsec_JML_45(pm_jml_45_fsi), linestyle='-', marker='None', color='#0000ff', logy=True, label='Paris FSI', zorder=2)
    
    l5, = B.plot_exp(pm_v18_45_pwia, f_red_pwiaXsec_V18_45(pm_v18_45_pwia), linestyle='--', marker='None', color='#009000', logy=True, label='AV18 PWIA', zorder=2)   
    l6, = B.plot_exp(pm_v18_45_fsi,  f_red_fsiXsec_V18_45(pm_v18_45_fsi), linestyle='-', marker='None', color='#009000', logy=True, label='AV18 FSI', zorder=2) 
    
    l7, = B.plot_exp(pm_cd_45_pwia, f_red_pwiaXsec_CD_45(pm_cd_45_pwia), linestyle='--', marker='None', color='#ff00ff', logy=True, label='CD-Bonn PWIA', zorder=2)     
    l8, = B.plot_exp(pm_cd_45_fsi,  f_red_fsiXsec_CD_45(pm_cd_45_fsi), linestyle='-', marker='None', color='#ff00ff', logy=True, label='CD-Bonn FSI', zorder=2) 

    l9, = B.plot_exp(pm_wjc2_45_pwba, f_red_pwbaXsec_WJC2_45(pm_wjc2_45_pwba), linestyle='--', marker='None', color='orange', logy=True, label='WJC2 (GKex05) PWBA', zorder=0 )
    l10, = B.plot_exp(pm_wjc2_45_fsi, f_red_fsiXsec_WJC2_45(pm_wjc2_45_fsi), linestyle='-', marker='None', color='orange', logy=True, label='WJC2 (GKex05) DWBA', zorder=0 )
    
    #Set axis limits
    B.pl.xlim(0.0, 1.19)
    B.pl.ylim(1e-7, 10.)

   
    #Set Tick Marks
    ax1.tick_params(axis='both', which='both', direction='in', top=True, bottom=True, left=True, right=True, labelsize=label_size)

    #Set Axes Labels for subplot 1
    B.pl.xlabel(r'$p_{\mathrm{r}}$ (GeV/c)', fontsize=axes_fontsize,  labelpad=10)                                                                                                                                                                                               
    B.pl.ylabel(r'', fontsize=axes_fontsize,  labelpad=10)                                                                                           
    B.pl.title('')

    #=====================
    #= THETA_NQ = 75 DEG =
    #====================
    ax2 = B.pl.subplot(gs[2], sharey=ax0)

    B.pl.text(0.7, 0.5e-6, r'$\theta_{nq}=75\pm5^{o}$', fontsize=font_size)
    B.pl.text(1.05, 0.9, r'(c)', fontsize=font_size)

    #Remove un-necessary Y tick labels from subplot
    B.pl.setp(ax2.get_yticklabels(), visible=False)
    
    #Plot Experimental Data (Hall A or Hall C)
    l1 = B.plot_exp(pm_avg75, red_dataXsec_avg_masked75, red_dataXsec_avg_tot_err75, marker='o', markersize=5, color='k', capsize=0, markerfacecolor='k', logy=True, label='This Experiment (Hall C)', zorder=4)
    l2 = B.plot_exp(pm_ha75, red_dataXsec_ha75, red_dataXsec_err_ha75, marker='s',  markersize=5, color='#ff1000', markerfacecolor='white', capsize=0, logy=True,  label='Hall A Data', zorder=3)

    #Plot theoretical calculations
    l3, = B.plot_exp(pm_jml_75_pwia, f_red_pwiaXsec_JML_75(pm_jml_75_pwia), linestyle='--', marker='None', color='#0000ff', logy=True, label='Paris PWIA', zorder=2)
    l4, = B.plot_exp(pm_jml_75_fsi,  f_red_fsiXsec_JML_75(pm_jml_75_fsi), linestyle='-', marker='None', color='#0000ff', logy=True, label='Paris FSI', zorder=2)
    
    l5, = B.plot_exp(pm_v18_75_pwia[0:len(pm_v18_75_pwia)-4], (f_red_pwiaXsec_V18_75(pm_v18_75_pwia))[0:len(pm_v18_75_pwia)-4], linestyle='--', marker='None', color='#009000', logy=True, label='AV18 PWIA', zorder=2)   
    l6, = B.plot_exp(pm_v18_75_fsi[0:len(pm_v18_75_fsi)-4],  (f_red_fsiXsec_V18_75(pm_v18_75_fsi))[0:len(pm_v18_75_fsi)-4], linestyle='-', marker='None', color='#009000', logy=True, label='AV18 FSI', zorder=2) 
    
    l7, = B.plot_exp(pm_cd_75_pwia[0:len(pm_cd_75_pwia)-4], (f_red_pwiaXsec_CD_75(pm_cd_75_pwia))[0:len(pm_cd_75_pwia)-4], linestyle='--', marker='None', color='#ff00ff', logy=True, label='CD-Bonn PWIA', zorder=2)     
    l8, = B.plot_exp(pm_cd_75_fsi[0:len(pm_cd_75_fsi)-4],  (f_red_fsiXsec_CD_75(pm_cd_75_fsi))[0:len(pm_cd_75_fsi)-4], linestyle='-', marker='None', color='#ff00ff', logy=True, label='CD-Bonn FSI', zorder=2) 

    l9, = B.plot_exp(pm_wjc2_75_pwba, f_red_pwbaXsec_WJC2_75(pm_wjc2_75_pwba), linestyle='--', marker='None', color='orange', logy=True, label='WJC2 (GKex05) PWBA', zorder=0 )
    l10, = B.plot_exp(pm_wjc2_75_fsi, f_red_fsiXsec_WJC2_75(pm_wjc2_75_fsi), linestyle='-', marker='None', color='orange', logy=True, label='WJC2 (GKex05) DWBA', zorder=0 )
    
    #Set axis limits
    B.pl.xlim(0.0, 1.19)
    B.pl.ylim(1e-7, 10.)

   
    #Set Tick Marks
    ax2.tick_params(axis='both', which='both', direction='in', top=True, bottom=True, left=True, right=True, labelsize=label_size)

    #Set Axes Labels for subplot 2
    B.pl.xlabel('')                                                                                                                                                                                                   
    B.pl.ylabel(r'', fontsize=axes_fontsize,  labelpad=10)                                                                                           
    B.pl.title('')
    
    #====================================
    #= SUBPLOTS SIZE / SPACING / LEGEND =
    #====================================

    leg_font_size = 9.2
    #Remove spacing between subplots
    plt.subplots_adjust(wspace = 0.0000000001, bottom=0.16, top=0.97, left=0.085, right=0.98) #, hspace = 0.001, wspace = 0.001)
    line_labels=['This Experiment (Hall C)', 'Hall A Data', 'JML Paris PWIA', 'JML Paris FSI', 'MS AV18 PWBA', 'MS AV18 FSI', 'MS CD-Bonn PWBA', 'MS CD-Bonn FSI', 'JVO WJC2 PWBA', 'JVO WJC2 FSI']
    ax2.legend([l1, l2, l3, l4, l5, l6, l7, l8, l9, l10], line_labels, loc='upper right', frameon=False, fontsize=leg_font_size)      #subplot to use for common legend


    #B.pl.show()
    B.pl.savefig('./PRL_plot1.pdf', format='pdf')



    
    #------------------------------------------------------------------------------------------
    #--------MAKE PRL PLOT 2 (Reduced Cross Sections Ratio vs. recoil momenta)-----------------
    #------------------------------------------------------------------------------------------
        
    #-----Create Subplots-----
    fig = B.pl.subplots(3, sharex=True, sharey=True, figsize=(6.7, 5.5))    #ORIGINAL SIZE: (width, height) = (6.7, 10)
    gs = gridspec.GridSpec(3, 1) 

    #=====================
    #= THETA_NQ = 35 DEG =
    #=====================
    
    ax0 = B.pl.subplot(gs[0])

    #B.pl.text(0.1, 4, r'(a)', fontsize=19)  #original txt label
    B.pl.text(1.05, 10, r'(a)', fontsize=font_size)  #with inset txt label

    #Remove un-necessary X tick labels from subplot
    B.pl.setp(ax0.get_xticklabels(), visible=False)

    #Plot the Data (and all models) to CD-Bonn PWIA model
    B.pl.axhline(y=1.0, xmin = 0.0, xmax=1.2, color='#ff00ff', linestyle='--', label='MS CD-Bonn PWBA',zorder=2)
    
    B.plot_exp(pm_ha35, R_HAdata35, R_HAdata_err35, marker='s', color='#ff0000', markerfacecolor='white',  label='Hall A Data', capsize=0, zorder=3)
    B.plot_exp(pm_avg35, R_data35, R_data_err35, marker='o', color='k', label='This Experiment (Hall C)', capsize=0, zorder=4)
    B.plot_exp(pm_avg35, R_CD_fsi35, marker='None', linestyle='-', color='#ff00ff', label='MS CD-Bonn FSI', zorder=1)
    B.plot_exp(pm_avg35, R_JML_pwia35, marker='None', linestyle='--', color='#0000ff', label='JML Paris PWIA', zorder=1)
    B.plot_exp(pm_avg35, R_JML_fsi35, marker='None', linestyle='-', color='#0000ff', label='JML Paris FSI', zorder=1)
    
    B.plot_exp(pm_avg35, R_V18_pwia35, marker='None', linestyle='--', color='#009000', label='MS AV18 PWBA', zorder=1)
    B.plot_exp(pm_avg35, R_V18_fsi35, marker='None', linestyle='-', color='#009000', label='MS AV18 FSI', zorder=1)

    B.plot_exp(pm_avg35, R_WJC2_pwba35, marker='None', linestyle='--', color='orange', label='JVO WJC2 PWBA', zorder=1)
    B.plot_exp(pm_avg35, R_WJC2_fsi35, marker='None', linestyle='-', color='orange', label='JVO WJC2 FSI', zorder=1)
    
    #Set axis limits
    B.pl.xlim(0.0, 1.2)
    B.pl.ylim(-0.5, 5.8)  #ORIGINAL
    B.pl.ylim(-0.5, 16.8)  #INSET
    
    #log settings
    B.pl.yscale('linear')

    #Set Tick Marks
    ax0.tick_params(axis='both', which='both', direction='in', top=True, bottom=True, left=True, right=True, labelsize=label_size)    

    #Set Axes Labels for subplot 0
    B.pl.xlabel('')                                                                                                                                                                                                   
    B.pl.ylabel('')                                                                                                                                                                        
    B.pl.title('') 

    #=================
    # INSET: 35 DEG
    #================
    if(plot2_inset==1):
        #----CRTE INSET PLOT(35 DEG)
        axins35 = inset_axes(ax0, width='60%', height='60%', loc='upper left')
        #B.pl.xlim(-0.1, 1.2)  --ORIGINAL INSET RANGE
        #B.pl.ylim(0.3, 4.5)   --ORIGINAL INSET RANGE

        B.pl.xlim(-0.1, 0.700)  #--INSET RANGE (M.K.Jones)
        B.pl.ylim(0.3, 2.08)   #--INSET RANGE (M.K.Jones)
        
        
        #Plot the Data (and all models) to CD-Bonn PWIA model
        B.pl.axhline(y=1.0, xmin = 0.0, xmax=1.2, color='#ff00ff', linestyle='--', label='MS CD-Bonn PWBA',zorder=2)
       
        B.plot_exp(pm_ha35, R_HAdata35, R_HAdata_err35, marker='s', markersize=4, color='#ff0000', markerfacecolor='white',  label='Hall A Data', capsize=0, zorder=3)
        B.plot_exp(pm_avg35, R_data35, R_data_err35, marker='o', markersize=4, color='k', label='This Experiment (Hall C)', capsize=0, zorder=4)
        B.plot_exp(pm_avg35, R_CD_fsi35, marker='None', linestyle='-', color='#ff00ff', label='MS CD-Bonn FSI', zorder=1)
        B.plot_exp(pm_avg35, R_JML_pwia35, marker='None', linestyle='--', color='#0000ff', label='JML Paris PWIA', zorder=1)
        B.plot_exp(pm_avg35, R_JML_fsi35, marker='None', linestyle='-', color='#0000ff', label='JML Paris FSI', zorder=1)
        
        B.plot_exp(pm_avg35, R_V18_pwia35, marker='None', linestyle='--', color='#009000', label='MS AV18 PWBA', zorder=1)
        B.plot_exp(pm_avg35, R_V18_fsi35, marker='None', linestyle='-', color='#009000', label='MS AV18 FSI', zorder=1)
        
        B.plot_exp(pm_avg35, R_WJC2_pwba35, marker='None', linestyle='--', color='orange', label='JVO WJC2 PWBA', zorder=1)
        B.plot_exp(pm_avg35, R_WJC2_fsi35, marker='None', linestyle='-', color='orange', label='JVO WJC2 FSI', zorder=1)
        
        axins35.yaxis.tick_right()
    
        
        B.pl.xlabel('')
        B.pl.ylabel('')
        B.pl.title('')
    
    #=====================
    #= THETA_NQ = 45 DEG =
    #=====================
    
    ax1 = B.pl.subplot(gs[1])

    #B.pl.text(0.1, 16, r'(b)', fontsize=19)  #original txt label
    B.pl.text(1.05, 18, r'(b)', fontsize=font_size)  #with inset txt label

    #Remove un-necessary X tick labels from subplot
    B.pl.setp(ax1.get_xticklabels(), visible=False)

    #Plot the Data (and all models) to CD-Bonn PWIA model
    B.pl.axhline(y=1.0, xmin = 0.0, xmax=1.2, color='#ff00ff', linestyle='--', label='MS CD-Bonn PWBA',zorder=2)
    
    B.plot_exp(pm_ha45, R_HAdata45, R_HAdata_err45, marker='s', color='#ff0000', markerfacecolor='white',  label='Hall A Data', capsize=0, zorder=3)
    B.plot_exp(pm_avg45, R_data45, R_data_err45, marker='o', color='k', label='This Experiment (Hall C)', capsize=0, zorder=4)
    B.plot_exp(pm_avg45, R_CD_fsi45, marker='None', linestyle='-', color='#ff00ff', label='MS CD-Bonn FSI', zorder=1)
    B.plot_exp(pm_avg45, R_JML_pwia45, marker='None', linestyle='--', color='#0000ff', label='JML Paris PWIA', zorder=1)
    B.plot_exp(pm_avg45, R_JML_fsi45, marker='None', linestyle='-', color='#0000ff', label='JML Paris FSI', zorder=1)
    
    B.plot_exp(pm_avg45, R_V18_pwia45, marker='None', linestyle='--', color='#009000', label='MS AV18 PWBA', zorder=1)
    B.plot_exp(pm_avg45, R_V18_fsi45, marker='None', linestyle='-', color='#009000', label='MS AV18 FSI', zorder=1)

    B.plot_exp(pm_avg45, R_WJC2_pwba45, marker='None', linestyle='--', color='orange', label='JVO WJC2 PWBA', zorder=1)
    B.plot_exp(pm_avg45, R_WJC2_fsi45, marker='None', linestyle='-', color='orange', label='JVO WJC2 FSI', zorder=1)
    
    #Set axis limits
    B.pl.xlim(0.0, 1.2)
    B.pl.ylim(-0.5, 22.9) #ORIGINAL

    #log settings
    #B.pl.ylim(0.5, 2.3)
    B.pl.yscale('linear')

    #Set Tick Marks
    ax1.tick_params(axis='both', which='both', direction='in', top=True, bottom=True, left=True, right=True, labelsize=label_size)    

    #Set Axes Labels for subplot 0
    B.pl.xlabel('')                                                                                                                                                                                                   
    B.pl.ylabel(r'$R = \sigma_{\mathrm{red}} / \sigma^{\textrm{\large CD-Bonn PWIA}}_{\mathrm{red}}$', fontsize=axes_fontsize,  labelpad=10)                                                                                                                                                                        
    B.pl.title('') 

    #=================
    # INSET: 45 DEG
    #================
    if(plot2_inset==1):
        #----CRTE INSET PLOT(45 DEG)
        axins45 = inset_axes(ax1, width='60%', height='60%', loc='upper left')
        #B.pl.xlim(-0.1, 1.2) --ORIGINAL INSET RANGE
        #B.pl.ylim(0.3, 4.5)  --ORIGINAL INSET RANGE

        B.pl.xlim(-0.1, 0.700)  #--INSET RANGE (M.K.Jones)
        B.pl.ylim(0.3, 2.08)   #--INSET RANGE (M.K.Jones)
        
        #Plot the Data (and all models) to CD-Bonn PWIA model
        B.pl.axhline(y=1.0, xmin = 0.0, xmax=1.2, color='#ff00ff', linestyle='--', label='MS CD-Bonn PWBA',zorder=2)
       
        B.plot_exp(pm_ha45, R_HAdata45, R_HAdata_err45, marker='s', markersize=4, color='#ff0000', markerfacecolor='white',  label='Hall A Data', capsize=0, zorder=3)
        B.plot_exp(pm_avg45, R_data45, R_data_err45, marker='o', markersize=4, color='k', label='This Experiment (Hall C)', capsize=0, zorder=4)
        B.plot_exp(pm_avg45, R_CD_fsi45, marker='None', linestyle='-', color='#ff00ff', label='MS CD-Bonn FSI', zorder=1)
        B.plot_exp(pm_avg45, R_JML_pwia45, marker='None', linestyle='--', color='#0000ff', label='JML Paris PWIA', zorder=1)
        B.plot_exp(pm_avg45, R_JML_fsi45, marker='None', linestyle='-', color='#0000ff', label='JML Paris FSI', zorder=1)
        
        B.plot_exp(pm_avg45, R_V18_pwia45, marker='None', linestyle='--', color='#009000', label='MS AV18 PWBA', zorder=1)
        B.plot_exp(pm_avg45, R_V18_fsi45, marker='None', linestyle='-', color='#009000', label='MS AV18 FSI', zorder=1)
        
        B.plot_exp(pm_avg45, R_WJC2_pwba45, marker='None', linestyle='--', color='orange', label='JVO WJC2 PWBA', zorder=1)
        B.plot_exp(pm_avg45, R_WJC2_fsi45, marker='None', linestyle='-', color='orange', label='JVO WJC2 FSI', zorder=1)
        
        axins45.yaxis.tick_right()
    
        
        B.pl.xlabel('')
        B.pl.ylabel('')
        B.pl.title('')
    
    #=====================
    #= THETA_NQ = 75 DEG =
    #=====================
    
    ax2 = B.pl.subplot(gs[2])

    B.pl.text(0.1, 4, r'(c)', fontsize=font_size)  #original txt label
    #B.pl.text(1.1, 15, r'(c)', fontsize=19)  #with inset txt label

    #Remove un-necessary X tick labels from subplot
    B.pl.setp(ax1.get_xticklabels(), visible=False)

    #Plot the Data (and all models) to CD-Bonn PWIA model
    B.pl.axhline(y=1.0, xmin = 0.0, xmax=0.5, color='#ff00ff', linestyle='--', label='MS CD-Bonn PWBA',zorder=2)
    
    B.plot_exp(pm_ha75, R_HAdata75, R_HAdata_err75, marker='s', color='#ff0000', markerfacecolor='white',  label='Hall A Data', capsize=0, zorder=3)
    B.plot_exp(pm_avg75, R_data75, R_data_err75, marker='o', color='k', label='This Experiment (Hall C)', capsize=0, zorder=4)
    B.plot_exp(pm_avg75, R_CD_fsi75, marker='None', linestyle='-', color='#ff00ff', label='MS CD-Bonn FSI', zorder=1)
    B.plot_exp(pm_avg75, R_JML_pwia75, marker='None', linestyle='--', color='#0000ff', label='JML Paris PWIA', zorder=1)
    B.plot_exp(pm_avg75, R_JML_fsi75, marker='None', linestyle='-', color='#0000ff', label='JML Paris FSI', zorder=1)
    
    B.plot_exp(pm_avg75, R_V18_pwia75, marker='None', linestyle='--', color='#009000', label='MS AV18 PWBA', zorder=1)
    B.plot_exp(pm_avg75, R_V18_fsi75, marker='None', linestyle='-', color='#009000', label='MS AV18 FSI', zorder=1)

    B.plot_exp(pm_avg75, R_WJC2_pwba75, marker='None', linestyle='--', color='orange', label='JVO WJC2 PWBA', zorder=1)
    B.plot_exp(pm_avg75, R_WJC2_fsi75, marker='None', linestyle='-', color='orange', label='JVO WJC2 FSI', zorder=1)
    
    #Set axis limits
    B.pl.xlim(0.0, 1.2)
    B.pl.ylim(-0.5, 5.8)  #ORIGINAL

    #log settings
    #B.pl.ylim(0.5, 2.3)
    B.pl.yscale('linear')

    #Set Tick Marks
    ax2.tick_params(axis='both', which='both', direction='in', top=True, bottom=True, left=True, right=True, labelsize=font_size)    

    #Set Axes Labels for subplot 0
    B.pl.xlabel(r'$p_{\mathrm{r}}$ (GeV/c)', fontsize=axes_fontsize,  labelpad=10)                                                                                                                                                                                                   
    B.pl.ylabel('')                                                                                                                                                                        
    B.pl.title('') 
    
    #====================================
    #= SUBPLOTS SIZE / SPACING / LEGEND =
    #====================================
    
    #Remove spacing between subplots
    plt.subplots_adjust(hspace = 0.00000001, bottom=0.12, top=0.98, right=0.95, left=0.18) #, hspace = 0.001, wspace = 0.001)
    B.pl.legend(loc='upper right', fontsize=7.5, frameon=False)

    #B.pl.show()
    B.pl.savefig('./PRL_plot2.pdf', format='pdf')
    



def main():
    print('Entering Main . . .')

    make_prl_plots(1)

if __name__=="__main__":
    main()


