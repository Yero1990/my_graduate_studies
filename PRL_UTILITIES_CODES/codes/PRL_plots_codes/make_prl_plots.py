import LT.box as B
from LT.datafile import dfile
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import least_squares
from scipy.optimize import curve_fit
import numpy.ma as ma
import sys                                     
import os                                                                                                       
from sys import argv  
import matplotlib
from matplotlib import rc
from matplotlib import *
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, InsetPosition, mark_inset)
from prl_utilities import *
#from prl_utilities_test import *


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

    #Read projected relative errors for full deuteron experiment (from simulations, at 120, 700,800,900 MeV/c)
    #  (to compare with PRL data and see how much improvement there is.
    kin35 = dfile('d2_projected_errors_thnq35deg.txt')
    kin45 = dfile('d2_projected_errors_thnq45deg.txt')
    kin75 = dfile('d2_projected_errors_thnq75deg.txt')
    rel_stats_err35 = kin35['rel_stats_err']
    rel_stats_err45 = kin45['rel_stats_err']
    rel_stats_err75 = kin75['rel_stats_err']

    pm_bin_35_proj = kin35['pm_bin']
    pm_bin_45_proj = kin45['pm_bin']
    pm_bin_75_proj = kin75['pm_bin']
    
    rel_stats_err35_m =  np.ma.array(rel_stats_err35, mask=(rel_stats_err35>0.5) | (rel_stats_err35==0))        
    rel_stats_err35_m =  np.ma.filled(rel_stats_err35_m.astype(float), np.nan )        
    rel_stats_err45_m =  np.ma.array(rel_stats_err45, mask=(rel_stats_err45>0.5) | (rel_stats_err45==0))
    rel_stats_err45_m = np.ma.filled(rel_stats_err45_m.astype(float), np.nan)
    rel_stats_err75_m =  np.ma.array(rel_stats_err75, mask=(rel_stats_err75>0.5) | (rel_stats_err75==0))
    rel_stats_err75_m = np.ma.filled(rel_stats_err75_m.astype(float), np.nan)

    
    #Read all Hall C data (without masking >50 % uncertainty, to be able to compare most bins with the projected uncertainties)
    red_dataXsec_avg35, red_dataXsec_avg_stats35 = read_hallc_data(35, verbose=True)
    red_dataXsec_avg45, red_dataXsec_avg_stats45 = read_hallc_data(45, verbose=True)
    red_dataXsec_avg75, red_dataXsec_avg_stats75 = read_hallc_data(75, verbose=True)

    #Calculate statistical relative errors (requires <=50 %)
    rel_dataXsec_err35 = red_dataXsec_avg_stats35 / red_dataXsec_avg35
    rel_dataXsec_err45 = red_dataXsec_avg_stats45 / red_dataXsec_avg45
    rel_dataXsec_err75 = red_dataXsec_avg_stats75 / red_dataXsec_avg75
    
    rel_dataXsec_err35_m = np.ma.array(rel_dataXsec_err35, mask=(rel_dataXsec_err35>0.5))
    rel_dataXsec_err35_m = np.ma.filled(rel_dataXsec_err35_m.astype(float), np.nan)

    rel_dataXsec_err45_m = np.ma.array(rel_dataXsec_err45, mask=(rel_dataXsec_err45>0.5))
    rel_dataXsec_err45_m = np.ma.filled(rel_dataXsec_err45_m.astype(float), np.nan)

    rel_dataXsec_err75_m = np.ma.array(rel_dataXsec_err75, mask=(rel_dataXsec_err75>0.5))
    rel_dataXsec_err75_m = np.ma.filled(rel_dataXsec_err75_m.astype(float), np.nan)
    
    #Read Hall C (E12-10-003) experimental data (requires <=50% statistical uncertainty, units: fm^3, GeV/c)
    pm_bin35, pm_avg35, red_dataXsec_avg_masked35, red_dataXsec_avg_tot_err35 = read_hallc_data(35)
    pm_bin45, pm_avg45, red_dataXsec_avg_masked45, red_dataXsec_avg_tot_err45 = read_hallc_data(45)
    pm_bin75, pm_avg75, red_dataXsec_avg_masked75, red_dataXsec_avg_tot_err75 = read_hallc_data(75)

    #Mask pm_avg bins
    pm_avg35_m = np.ma.array(pm_avg35, mask=(rel_dataXsec_err35_m>0.5) | (rel_dataXsec_err35_m==0.))
    pm_avg35_m = np.ma.filled(pm_avg35_m.astype(float), np.nan)
    pm_avg45_m = np.ma.array(pm_avg45, mask=(rel_dataXsec_err45_m>0.5) | (rel_dataXsec_err45_m==0.))
    pm_avg45_m = np.ma.filled(pm_avg45_m.astype(float), np.nan)
    pm_avg75_m = np.ma.array(pm_avg75, mask=(rel_dataXsec_err75_m>0.5) | (rel_dataXsec_err75_m==0.))
    pm_avg75_m = np.ma.filled(pm_avg75_m.astype(float), np.nan)
    
    #Mask pm_bins
    pm_bin35_m = np.ma.array(pm_bin35, mask=((rel_stats_err35>0.5) | (rel_stats_err35==0)))
    pm_bin35_m = np.ma.filled(pm_bin35_m.astype(float), np.nan)

    pm_bin45_m = np.ma.array(pm_bin45, mask=((rel_stats_err45>0.5) | (rel_stats_err45==0)))
    pm_bin45_m = np.ma.filled(pm_bin45_m.astype(float), np.nan)

    pm_bin75_m = np.ma.array(pm_bin75, mask=((rel_stats_err75>0.5) | (rel_stats_err75==0)))
    pm_bin75_m = np.ma.filled(pm_bin75_m.astype(float), np.nan)
    
    #Calculate the projected absolute stats. uncertainty 
    red_dataXsec_proj_stats35 = red_dataXsec_avg_masked35 * rel_stats_err35_m
    red_dataXsec_proj_stats45 = red_dataXsec_avg_masked45 * rel_stats_err45_m
    red_dataXsec_proj_stats75 = red_dataXsec_avg_masked75 * rel_stats_err75_m
    
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


    #print(pm_v18_45_fsi, f_red_fsiXsec_V18_45(pm_v18_45_fsi))
    #======== Define Ratios of Reduced Cross Sections (to be used in PRL PLOT 2)=====
    
    #Data (and all models) ratio to CD-Bonn PWIA

    #=====================
    #= THETA_NQ = 35 DEG =
    #=====================

    R_data35 = red_dataXsec_avg_masked35 / f_red_pwiaXsec_CD_35(pm_avg35)
    R_data_err35 = red_dataXsec_avg_tot_err35 / f_red_pwiaXsec_CD_35(pm_avg35)
    R_data_proj_err35 = red_dataXsec_proj_stats35 / f_red_pwiaXsec_CD_35(pm_avg35)

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
    R_data_proj_err45 = red_dataXsec_proj_stats45 / f_red_pwiaXsec_CD_45(pm_avg45)

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
    R_data_proj_err75 = red_dataXsec_proj_stats75 / f_red_pwiaXsec_CD_75(pm_avg75)

    R_HAdata75 = red_dataXsec_ha75 / f_red_pwiaXsec_CD_75(pm_ha75)
    R_HAdata_err75 = red_dataXsec_err_ha75 / f_red_pwiaXsec_CD_75(pm_ha75)

    R_CD_fsi75    = f_red_fsiXsec_CD_75(pm_avg75)    / f_red_pwiaXsec_CD_75(pm_avg75)
    R_JML_pwia75  = f_red_pwiaXsec_JML_75(pm_avg75)  / f_red_pwiaXsec_CD_75(pm_avg75)
    R_JML_fsi75   = f_red_fsiXsec_JML_75(pm_avg75)   / f_red_pwiaXsec_CD_75(pm_avg75)
    R_V18_pwia75  = f_red_pwiaXsec_V18_75(pm_avg75)  / f_red_pwiaXsec_CD_75(pm_avg75)
    R_V18_fsi75   = f_red_fsiXsec_V18_75(pm_avg75)   / f_red_pwiaXsec_CD_75(pm_avg75)
    R_WJC2_pwba75 = f_red_pwbaXsec_WJC2_75(pm_avg75) / f_red_pwiaXsec_CD_75(pm_avg75)
    R_WJC2_fsi75  = f_red_fsiXsec_WJC2_75(pm_avg75)  / f_red_pwiaXsec_CD_75(pm_avg75)

    # DEFINE PHENOMENOLOGICAL FIT FUNCTION FOR DATA (For the purpose of projecting the data points of the new deuteron experiment)
    
    #monotonic exponential function with  y = a * exp^ (-m * x) + b
    def monoExp(x, a, m, b):
        
        return a * np.exp(-m * x) + b


    #set initial parameters
    p0_init_35 = (20, 30, 1e-14)

    #mask infs or nans in the redXsec, and pass masked values to pm_avg, and redXsec_err (then compress, to get rid of masked elements for fitting)
    redXsec_fit35 = ma.masked_invalid(red_dataXsec_avg_masked35) 
    redXsec_err_fit35 = ma.masked_array(red_dataXsec_avg_tot_err35, redXsec_fit35.mask)
    pm_avg35_fit = ma.masked_array(pm_avg35, redXsec_fit35.mask)

    #compress masked values (eliminate -- )
    redXsec_fit35 = ma.compressed( redXsec_fit35 )
    redXsec_err_fit35 = ma.compressed(redXsec_err_fit35)
    pm_avg35_fit = ma.compressed( pm_avg35_fit )

    #for optimal phenomenological fit, combine hall a and hall c data points for smooth transition
    pm_avg35_fit_total = np.concatenate((pm_ha35[6:],pm_avg35_fit[6:]))
    redXsec_fit35_total = np.concatenate((red_dataXsec_ha35[6:], redXsec_fit35[6:]))
    redXsec_err_fit35_total = np.concatenate((red_dataXsec_err_ha35[6:],redXsec_err_fit35[6:]))

    print('pm35-----> ', pm_avg35_fit_total)
    #perform fit
    #fit_parms_35_low, cov_35_low = curve_fit(f=monoExp, xdata=pm_avg35_fit[0:7], ydata=redXsec_fit35[0:7], sigma=redXsec_err_fit35[0:7], p0=p0_init_low35)
    #fit_parms_35_hi, cov_35_hi = curve_fit(f=monoExp, xdata=pm_avg35_fit[7:20], ydata=redXsec_fit35[7:20], sigma=redXsec_err_fit35[7:20], p0=p0_init_hi35)

    fit_parms_35, cov_fit_35 = curve_fit(f=monoExp, xdata=pm_avg35_fit_total, ydata=redXsec_fit35_total, sigma=redXsec_err_fit35_total, p0=p0_init_35)

    #polyfit
    #fit_parms_poly = np.poly_fit(f=monoExp, xdata=pm_avg35_fit_total, ydata=redXsec_fit35_total, sigma=redXsec_err_fit35_total, p0=p0_init_low35)
    
    #a35_low, m35_low, b35_low = fit_parms_35_low
    #a35_hi, m35_hi, b35_hi = fit_parms_35_hi

    #a35, m35, b35 = fit_parms_35
    print('fit parms 35 = ', *fit_parms_35)
    #print('a35 = ',a35)
    #print('m35 = ',m35)
    #print('b35 = ',b35)
    #y35_fit_model_low = monoExp(pm_avg35[0:7], a35_low, m35_low, b35_low)
    #y35_fit_model_hi = monoExp(pm_avg35[10:30], a35_hi, m35_hi, b35_hi)

    y35_fit_model = monoExp(pm_bin35[7:], *fit_parms_35)

    '''
    #set initial parameters (for 45 deg)
    p0_init_low45 = (1e-6, 1e-5, 1e-6)
    p0_init_hi45 = (1e-3, 1, 1e-6)

    #mask infs or nans in the redXsec, and pass masked values to pm_avg, and redXsec_err (then compress, to get rid of masked elements for fitting)
    redXsec_fit45 = ma.masked_invalid(red_dataXsec_avg_masked45) 
    redXsec_err_fit45 = ma.masked_array(red_dataXsec_avg_tot_err45, redXsec_fit45.mask)
    pm_avg45_fit = ma.masked_array(pm_avg45, redXsec_fit45.mask)

    #compress masked values (eliminate -- )
    redXsec_fit45 = ma.compressed( redXsec_fit45 )
    redXsec_err_fit45 = ma.compressed(redXsec_err_fit45)
    pm_avg45_fit = ma.compressed( pm_avg45_fit )

    print('pm_avg45_fit ---->',pm_avg45_fit)
    print('redXsec_fit45 ----->',redXsec_fit45)
    print('redXsec_err_fit45 ----->',redXsec_err_fit45)
    
    #perform fit
    fit_parms_45_low, cov_45_low = curve_fit(f=monoExp, xdata=pm_avg45_fit[0:6], ydata=redXsec_fit45[0:6], sigma=redXsec_err_fit45[0:6], p0=p0_init_low45)
    fit_parms_45_hi, cov_45_hi = curve_fit(f=monoExp, xdata=pm_avg45_fit[7:20], ydata=redXsec_fit45[7:20], sigma=redXsec_err_fit45[7:20], p0=p0_init_hi45)
    
    a45_low, m45_low, b45_low = fit_parms_45_low
    a45_hi, m45_hi, b45_hi = fit_parms_45_hi
    
    y45_fit_model_low = monoExp(pm_avg45[0:7], a45_low, m45_low, b45_low)
    y45_fit_model_hi = monoExp(pm_avg45[10:30], a45_hi, m45_hi, b45_hi)
    '''
    
    #------------------------------------------------------------------------------------
    #--------MAKE PRL PLOT 1 (Reduced Cross Sections vs. recoil momenta)-----------------
    #------------------------------------------------------------------------------------
    
    #-----Create Subplots-----
    fig = B.pl.subplots(3, sharex=True, sharey=True, figsize=(13.4, 6))   #figsize (width, height), (13.4,6) is in inches;  PRL requirements: width = 8.6 cm (3.38583 in) --> 1 manuscript column (for wider figs, use 1.5 or 2 columns: (17.2 cm or 6.771654 in)
    gs = gridspec.GridSpec(1, 3) 

    #=====================
    #= THETA_NQ = 35 DEG =
    #====================
    ax0 = B.pl.subplot(gs[0])

    B.pl.text(0.25, 0.5e-6, r'$\theta_{nq}=35\pm5^{o}$', fontsize=19)
    B.pl.text(1.0, 1e0, r'(a)', fontsize=19)

    #Plot Experimental Data (Hall A or Hall C)
    
    #Hall C data (for PRL)
    l1 = B.plot_exp(pm_avg35, red_dataXsec_avg_masked35, red_dataXsec_avg_tot_err35, marker='o', markersize=2, color='lime', capsize=0, markerfacecolor='lime', logy=True, label='This Experiment (Hall C)', zorder=4)

    #lfit_low, = B.plot_exp(pm_avg35[0:7], y35_fit_model_low, linestyle='-', marker='', color='cyan', logy=True, zorder=3)
    #lfit_hi, = B.plot_exp(pm_avg35[10:30], y35_fit_model_hi, linestyle='-', marker='', color='cyan', logy=True,zorder=3)

    #plot fit model curve
    lfit_35, = B.plot_exp(pm_bin35[7:], y35_fit_model, linestyle='-', marker='', color='k', logy=False,zorder=4)

    #Hall C data (for statistical projections of full deuteron experiment)
    l0 = B.plot_exp(pm_bin35[13:29]+0.01, monoExp(pm_bin35[13:29]+0.01,  *fit_parms_35), monoExp(pm_bin35[13:29]+0.01, *fit_parms_35) *rel_stats_err35_m[13:29], marker='o', markersize=3, color='k', ecolor='k', capsize=0, logy=True, zorder=4)

          
    #Hall A data
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
    B.pl.xlim(0.0, 1.2)
    B.pl.ylim(1e-7, 10.)
 
    #Set Tick Marks
    ax0.tick_params(axis='both', which='both', direction='in', top=True, bottom=True, left=True, right=True, labelsize=19)

    #Set Axes Labels for subplot 0
    B.pl.xlabel('')                                                                                                                                                                                                   
    B.pl.ylabel(r'$\sigma_{\mathrm{red}}$ ($\mathrm{fm}^{3}$) ', fontsize=22,  labelpad=10)                                                                                           
    B.pl.title('')

    #Write to File
    print('LENGHTS (data, JML) -----> ', len(pm_bin35), len(pm_avg35), len(red_dataXsec_avg_masked35), len(red_dataXsec_avg_tot_err35), len(pm_jml_35_pwia), len(f_red_pwiaXsec_JML_35(pm_jml_35_pwia)), len(f_red_fsiXsec_JML_35(pm_jml_35_fsi)))
    print('LENGTHS (V18, CD-Bonn)--->', len(pm_v18_35_pwia), len(f_red_pwiaXsec_V18_35(pm_v18_35_pwia)), len(pm_v18_35_fsi), len(f_red_fsiXsec_V18_35(pm_v18_35_fsi)), len(pm_cd_35_pwia), len(f_red_pwiaXsec_CD_35(pm_cd_35_pwia)), len(pm_cd_35_fsi), len(f_red_fsiXsec_CD_35(pm_cd_35_fsi)), len(pm_wjc2_35_pwba), len(f_red_pwbaXsec_WJC2_35(pm_wjc2_35_pwba)), len(pm_wjc2_35_fsi), len(f_red_fsiXsec_WJC2_35(pm_wjc2_35_fsi)))
    print('COLUMN-STACK----> ', np.column_stack((pm_bin35, pm_avg35, pm_jml_35_pwia, pm_jml_35_fsi, pm_v18_35_pwia, pm_v18_35_fsi, pm_cd_35_pwia, pm_cd_35_fsi)))
    print('WJC2 COLUMN----> ', np.column_stack((pm_wjc2_35_pwba, pm_wjc2_35_pwba)))
    '''
    data_35 = np.column_stack([pm_bin35, pm_avg35,      red_dataXsec_avg_masked35, red_dataXsec_avg_tot_err35])
    jml_35 = np.column_stack([pm_bin35, pm_jml_35_pwia, f_red_pwiaXsec_JML_35(pm_jml_35_pwia), f_red_fsiXsec_JML_35(pm_jml_35_fsi) ])

    fout_list = ['data', 'jml_pwia', 'jml_fsi', 'v18_pwia', 'v18_fsi', 'cdbonn_pwia', 'cdbonn_fsi', 'wjc2_pwia', 'wjc2_fsi']
    fout_name = []
    
    for i in enumerate(fout_list):
        idx = i[0]
        prfx = i[1]
        
        fout_name.append('%s_redXsec_35.txt' % (prfx))
        fout = open(fout_name[idx], 'w')
        fout.write('# %s reduced cross sections\n'
            '# thnq             : 35 +/- 5 deg\n'
            '# pm_bin           : central bin value (GeV/c)\n'
            '# pm_avg           : avergae missing momentum (GeV/c)\n'
            '# red_%s_Xsec       : reduced %s cross sections (fm^3)\n'
            '# red_%s_Xsec_err   : reduced %s cross sections error (fm^3)\n'
            '# \n'
            '#! pm_bin[f,0]/ pm_avg[f,1]/ red_%s_Xsec[f,2]/  red_%s_Xsec_err[f,3]/  \n' % (prfx, prfx, prfx, prfx, prfx, prfx, prfx)
            )
        fout.close()

    #np.savetxt(fout_name, data1_35, fmt=['%.3f\t', '%.3f\t', '%.4E\t', '%.4E\t', '%.3f\t', '%.4E\t' ])
    #fout.close()
    '''
        
    #=====================
    #= THETA_NQ = 45 DEG =
    #====================
    ax1 = B.pl.subplot(gs[1], sharey=ax0)

    B.pl.text(0.25, 0.5e-6, r'$\theta_{nq}=45\pm5^{o}$', fontsize=19)
    B.pl.text(1.0, 1e0, r'(b)', fontsize=19)

    #Remove un-necessary Y tick labels from subplot
    B.pl.setp(ax1.get_yticklabels(), visible=False)

    #Hall C data (for statistical projections of full deuteron experiment)
    #l0 = B.plot_exp(pm_avg45_m+0.01, red_dataXsec_avg_masked45, red_dataXsec_proj_stats45, marker='^', ms=2, color='cyan', ecolor='cyan', capsize=0, logy=True, label='Projected Stat. Error', zorder=3)

    #Plot Experimental Data (Hall A or Hall C)
    l1 = B.plot_exp(pm_avg45, red_dataXsec_avg_masked45, red_dataXsec_avg_tot_err45, marker='o', markersize=2, color='k', capsize=0, markerfacecolor='k', logy=True, label='This Experiment (Hall C)', zorder=4)    


    #lfit_low, = B.plot_exp(pm_avg45[0:7], y45_fit_model_low, linestyle='-', marker='', color='cyan', logy=True, zorder=3)
    #lfit_hi, = B.plot_exp(pm_avg45[10:30], y45_fit_model_hi, linestyle='-', marker='', color='cyan', logy=True,zorder=3)

    l2 = B.plot_exp(pm_ha45, red_dataXsec_ha45, red_dataXsec_err_ha45, marker='s',  markersize=5, color='#ff1000', markerfacecolor='white', capsize=0, logy=True,  label='Hall A Data', zorder=3)

    #Hall C data (for statistical projections of full deuteron experiment)
    #l0_lo = B.plot_exp(pm_bin45[0:6]+0.01, monoExp(pm_bin45[0:6]+0.01, a45_low, m45_low, b45_low), monoExp(pm_bin45[0:6]+0.01, a45_low, m45_low, b45_low) *rel_stats_err45_m[0:6], marker='o', markersize=3, color='cyan', mfc='white', mec='cyan', ecolor='cyan', capsize=0, logy=True, label='Projected Stat. Error', zorder=4)
    #l0_hi = B.plot_exp(pm_bin45[14:29]+0.01, monoExp(pm_bin45[14:29]+0.01, a45_hi, m45_hi, b45_hi), monoExp(pm_bin45[14:29]+0.01, a45_hi, m45_hi, b45_hi) *rel_stats_err45_m[14:29], marker='o', markersize=3, color='cyan', mfc='white', mec='cyan', ecolor='cyan', capsize=0, logy=True, label='Projected Stat. Error', zorder=4)

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
    B.pl.xlim(0.0, 1.2)
    B.pl.ylim(1e-7, 10.)

   
    #Set Tick Marks
    ax1.tick_params(axis='both', which='both', direction='in', top=True, bottom=True, left=True, right=True, labelsize=19)

    #Set Axes Labels for subplot 1
    B.pl.xlabel(r'$p_{\mathrm{r}}$ (GeV/c)', fontsize=22,  labelpad=10)                                                                                                                                                                                               
    B.pl.ylabel(r'', fontsize=22,  labelpad=10)                                                                                           
    B.pl.title('')

    #=====================
    #= THETA_NQ = 75 DEG =
    #====================
    ax2 = B.pl.subplot(gs[2], sharey=ax0)

    B.pl.text(0.7, 0.5e-6, r'$\theta_{nq}=75\pm5^{o}$', fontsize=19)
    B.pl.text(1.0, 1e0, r'(c)', fontsize=19)

    #Remove un-necessary Y tick labels from subplot
    B.pl.setp(ax2.get_yticklabels(), visible=False)

    #Hall C data (for statistical projections of full deuteron experiment)
    #l0 = B.plot_exp(pm_avg75_m+0.01, red_dataXsec_avg_masked75, red_dataXsec_proj_stats75, marker='o', ms=2, color='cyan', ecolor='cyan', capsize=0, logy=True, label='Projected Stat. Error', zorder=3)
    
    #Plot Experimental Data (Hall A or Hall C)
    l1 = B.plot_exp(pm_avg75, red_dataXsec_avg_masked75, red_dataXsec_avg_tot_err75, marker='o', markersize=2, color='k', capsize=0, markerfacecolor='k', logy=True, label='This Experiment (Hall C)', zorder=4)
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
    B.pl.xlim(0.0, 1.2)
    B.pl.ylim(1e-7, 10.)

   
    #Set Tick Marks
    ax2.tick_params(axis='both', which='both', direction='in', top=True, bottom=True, left=True, right=True, labelsize=19)

    #Set Axes Labels for subplot 2
    B.pl.xlabel('')                                                                                                                                                                                                   
    B.pl.ylabel(r'', fontsize=22,  labelpad=10)                                                                                           
    B.pl.title('')
    
    #====================================
    #= SUBPLOTS SIZE / SPACING / LEGEND =
    #====================================

    #Remove spacing between subplots
    plt.subplots_adjust(wspace = 0.0000000001, bottom=0.13, top=0.98, left=0.085, right=0.98) #, hspace = 0.001, wspace = 0.001)
    #----PRL PLOT-----
    #line_labels=['This Experiment (Hall C)', 'Hall A Data', 'JML Paris PWIA', 'JML Paris FSI', 'MS AV18 PWBA', 'MS AV18 FSI', 'MS CD-Bonn PWBA', 'MS CD-Bonn FSI', 'JVO WJC2 PWBA', 'JVO WJC2 FSI']
    #ax2.legend([l1, l2, l3, l4, l5, l6, l7, l8, l9, l10], line_labels, loc='upper right', frameon=False, fontsize=12)      #subplot to use for common legend
    #----Projected Errors Plot----
    line_labels=[r'This Experiment (Hall C)', 'Hall A Data', 'JML Paris PWIA', 'JML Paris FSI', 'MS AV18 PWBA', 'MS AV18 FSI', 'MS CD-Bonn PWBA', 'MS CD-Bonn FSI', 'JVO WJC2 PWBA', 'JVO WJC2 FSI', 'Data Fit: f(x) = A$e^{-m x}$+b']
    ax2.legend([l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, lfit_35], line_labels, loc='upper right', frameon=False, fontsize=11)      #subplot to use for common legend


    B.pl.show()
    #B.pl.savefig('./PRL_plot1.pdf')



    
    #------------------------------------------------------------------------------------------
    #--------MAKE PRL PLOT 2 (Reduced Cross Sections Ratio vs. recoil momenta)-----------------
    #------------------------------------------------------------------------------------------
        
    #-----Create Subplots-----
    fig = B.pl.subplots(3, sharex=True, sharey=True, figsize=(6.7, 10))    #ORIGINAL SIZE: (width, height) = (6.7, 10)
    gs = gridspec.GridSpec(3, 1) 

    #=====================
    #= THETA_NQ = 35 DEG =
    #=====================
    
    ax0 = B.pl.subplot(gs[0])

    #B.pl.text(0.1, 4, r'(a)', fontsize=19)  #original txt label
    B.pl.text(1.05, 10, r'(a)', fontsize=19)  #with inset txt label

    #Remove un-necessary X tick labels from subplot
    B.pl.setp(ax0.get_xticklabels(), visible=False)

    #Plot the Data (and all models) to CD-Bonn PWIA model
    B.pl.axhline(y=1.0, xmin = 0.0, xmax=1.2, color='#ff00ff', linestyle='--', label='MS CD-Bonn PWBA',zorder=2)

    #Projected stats. error
    #B.plot_exp(pm_avg35+0.01, R_data35, R_data_proj_err35, marker='^', ms=2, color='cyan', ecolor='cyan', label='Projected Stat. Error', capsize=0, zorder=3)
        
    B.plot_exp(pm_ha35, R_HAdata35, R_HAdata_err35, marker='s', color='#ff0000', markerfacecolor='white',  label='Hall A Data', capsize=0, zorder=3)
    B.plot_exp(pm_avg35, R_data35, R_data_err35, marker='o', ms=2, color='k', label='This Experiment (Hall C)', capsize=0, zorder=4)
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
    ax0.tick_params(axis='both', which='both', direction='in', top=True, bottom=True, left=True, right=True, labelsize=19)    

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

        #Projected stats. error
        #B.plot_exp(pm_avg35+0.01, R_data35, R_data_proj_err35, marker='^', ms=2, color='cyan', ecolor='cyan', label='Projected Stat. Error', capsize=0, zorder=3)
                
        B.plot_exp(pm_ha35, R_HAdata35, R_HAdata_err35, marker='s', markersize=4, color='#ff0000', markerfacecolor='white',  label='Hall A Data', capsize=0, zorder=3)
        B.plot_exp(pm_avg35, R_data35, R_data_err35, marker='o', markersize=2, color='k', label='This Experiment (Hall C)', capsize=0, zorder=4)
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
    B.pl.text(1.05, 18, r'(b)', fontsize=19)  #with inset txt label

    #Remove un-necessary X tick labels from subplot
    B.pl.setp(ax1.get_xticklabels(), visible=False)

    #Plot the Data (and all models) to CD-Bonn PWIA model
    B.pl.axhline(y=1.0, xmin = 0.0, xmax=1.2, color='#ff00ff', linestyle='--', label='MS CD-Bonn PWBA',zorder=2)

    #Projected stats. error
    #B.plot_exp(pm_avg45+0.01, R_data45, R_data_proj_err45, marker='^', color='cyan', ms=2, ecolor='cyan', label='Projected Stat. Error', capsize=0, zorder=3)
        
    B.plot_exp(pm_ha45, R_HAdata45, R_HAdata_err45, marker='s', color='#ff0000', markerfacecolor='white',  label='Hall A Data', capsize=0, zorder=3)
    B.plot_exp(pm_avg45, R_data45, R_data_err45, marker='o', ms=2, color='k', label='This Experiment (Hall C)', capsize=0, zorder=4)
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
    ax1.tick_params(axis='both', which='both', direction='in', top=True, bottom=True, left=True, right=True, labelsize=19)    

    #Set Axes Labels for subplot 0
    B.pl.xlabel('')                                                                                                                                                                                                   
    B.pl.ylabel(r'$R = \sigma_{\mathrm{red}} / \sigma^{\textrm{\large CD-Bonn PWIA}}_{\mathrm{red}}$', fontsize=22,  labelpad=10)                                                                                                                                                                        
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

        B.plot_exp(pm_ha45, R_HAdata45, R_HAdata_err45, marker='s', ms=4, color='#ff0000', markerfacecolor='white',  label='Hall A Data', capsize=0, zorder=3)
        B.plot_exp(pm_avg45, R_data45, R_data_err45, marker='o', ms=2, color='k', label='This Experiment (Hall C)', capsize=0, zorder=4)

        #Projected stats. error
        #B.plot_exp(pm_avg45+0.01, R_data45, R_data_proj_err45, marker='^', ms=2, color='cyan', ecolor='cyan', label='Projected Stat. Error', capsize=0, zorder=3)

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

    B.pl.text(0.1, 4, r'(c)', fontsize=19)  #original txt label
    #B.pl.text(1.1, 15, r'(c)', fontsize=19)  #with inset txt label

    #Remove un-necessary X tick labels from subplot
    B.pl.setp(ax1.get_xticklabels(), visible=False)

    #Plot the Data (and all models) to CD-Bonn PWIA model
    B.pl.axhline(y=1.0, xmin = 0.0, xmax=0.5, color='#ff00ff', linestyle='--', label='MS CD-Bonn PWBA',zorder=2)

    B.plot_exp(pm_ha75, R_HAdata75, R_HAdata_err75, marker='s', color='#ff0000', markerfacecolor='white',  label='Hall A Data', capsize=0, zorder=3)
    B.plot_exp(pm_avg75, R_data75, R_data_err75, marker='o', ms=2, color='k', label='This Experiment (Hall C)', capsize=0, zorder=4)

    #Projected stats. error
    #B.plot_exp(pm_avg75+0.01, R_data75, R_data_proj_err75, marker='^', ms=2, color='cyan', ecolor='cyan', label='Projected Statistical Error', capsize=0, zorder=3)
    
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
    ax2.tick_params(axis='both', which='both', direction='in', top=True, bottom=True, left=True, right=True, labelsize=19)    

    #Set Axes Labels for subplot 0
    B.pl.xlabel(r'$p_{\mathrm{r}}$ (GeV/c)', fontsize=22,  labelpad=10)                                                                                                                                                                                                   
    B.pl.ylabel('')                                                                                                                                                                        
    B.pl.title('') 
    
    #====================================
    #= SUBPLOTS SIZE / SPACING / LEGEND =
    #====================================

    #Remove spacing between subplots
    plt.subplots_adjust(hspace = 0.00000001, bottom=0.09, top=0.98, right=0.95, left=0.18) #, hspace = 0.001, wspace = 0.001)
    B.pl.legend(loc='upper right', fontsize=12, frameon=False)

    #B.pl.show()
    B.pl.savefig('./PRL_plot2.pdf')

    
    '''
    #----------MAKE PLOT OF PROJECTED RELATIVE ERRORS COMPARISON TO COMMISSIONING DATA------------

    
    fig2 = B.pl.subplots(3, sharex=True, sharey=True, figsize=(13.4, 6))
    gs = gridspec.GridSpec(1, 3) 

    #=====================
    #= THETA_NQ = 35 DEG =
    #====================
    ax0 = B.pl.subplot(gs[0])
    l0 = plt.errorbar(pm_bin35_m+0.01, np.repeat(0, len(pm_bin35_m)), rel_stats_err35_m*100., color='cyan', linestyle='none', marker='^', ms=2, label = r'$I_{\textrm{beam}}$=70 $\mu A$' )
    l1 = plt.errorbar(pm_avg35_m, np.repeat(0, len(pm_avg35_m)), rel_dataXsec_err35_m*100., marker='o', color='k', ms=0, linestyle='none',  capsize=2)
    l2 = plt.errorbar(pm_avg35_m, np.repeat(0, len(pm_avg35_m)), (red_dataXsec_avg_tot_err35/red_dataXsec_avg_masked35)*100.,  marker='o', color='k', ms=0,  linestyle='none', capsize=2)
    #Set axis limits
    B.pl.xlim(0.0, 1.2)
    B.pl.ylim(-50, 78.)

    x_coord = [0, 1.2]
    y_coord_1 = [10, 10]
    y_coord_2 = [-10, -10]
    plt.plot(x_coord, y_coord_1, linestyle='dashed', color='k', linewidth=1)
    plt.plot(x_coord, y_coord_2, linestyle='dashed', color='k', linewidth=1)
 
    #Set Tick Marks
    ax0.tick_params(axis='both', which='both', direction='in', top=True, bottom=True, left=True, right=True, labelsize=19)

    #Set Axes Labels for subplot 0
    B.pl.xlabel('')                                                                                                                                                                                                   
    B.pl.ylabel(r'Relative Errors (\%)', fontsize=22,  labelpad=10)                                                                                           
    B.pl.title('')
    
    #=====================
    #= THETA_NQ = 45 DEG =
    #====================
    ax1 = B.pl.subplot(gs[1], sharey=ax0)
    l0 = plt.errorbar(pm_bin45_m+0.01, np.repeat(0, len(pm_bin45_m)), rel_stats_err45_m*100., color='cyan', linestyle='none', marker='^', ms=2 )
    l1 = plt.errorbar(pm_avg45_m, np.repeat(0, len(pm_avg45_m)), rel_dataXsec_err45_m*100., color='k', linestyle='none', marker='', capsize=2)
    l2 = plt.errorbar(pm_avg45_m, np.repeat(0, len(pm_avg45_m)), (red_dataXsec_avg_tot_err45/red_dataXsec_avg_masked45)*100., color='k', linestyle='none', marker='', capsize=2)

    #Remove un-necessary Y tick labels from subplot
    B.pl.setp(ax1.get_yticklabels(), visible=False)
    
    #Set axis limits
    B.pl.xlim(0.0, 1.2)
    B.pl.ylim(-50, 78.)

    x_coord = [0, 1.2]
    y_coord_1 = [10, 10]
    y_coord_2 = [-10, -10]
    plt.plot(x_coord, y_coord_1, linestyle='dashed', color='k', linewidth=1)
    plt.plot(x_coord, y_coord_2, linestyle='dashed', color='k', linewidth=1)
 
    #Set Tick Marks
    ax1.tick_params(axis='both', which='both', direction='in', top=True, bottom=True, left=True, right=True, labelsize=19)

    B.pl.xlabel(r'Missing Momentum, $P_{m}$ (GeV/c)', fontsize=19)
    #B.pl.text(0.25, 0.5e-6, r'Missing Momentum, $P_{m} (GeV/c)$', fontsize=19)

    #=====================
    #= THETA_NQ = 75 DEG =
    #====================
    ax2 = B.pl.subplot(gs[2])
    l0 = plt.errorbar(pm_bin75_m+0.01, np.repeat(0, len(pm_bin75_m)), rel_stats_err75_m*100., color='cyan', linestyle='none', marker='^', ms=2)
    l1 = plt.errorbar(pm_avg75_m, np.repeat(0, len(pm_avg75_m)), rel_dataXsec_err75_m*100., color='k', linestyle='none', marker='', capsize=2)
    l2 = plt.errorbar(pm_avg75_m, np.repeat(0, len(pm_avg75_m)), (red_dataXsec_avg_tot_err75/red_dataXsec_avg_masked75)*100., color='k', linestyle='none', marker='', capsize=2)

    #Remove un-necessary Y tick labels from subplot
    B.pl.setp(ax2.get_yticklabels(), visible=False)
    
    #Set axis limits
    B.pl.xlim(0.0, 1.2)
    B.pl.ylim(-50, 78.)

    x_coord = [0, 1.2]
    y_coord_1 = [10, 10]
    y_coord_2 = [-10, -10]
    plt.plot(x_coord, y_coord_1, linestyle='dashed', color='k', linewidth=1)
    plt.plot(x_coord, y_coord_2, linestyle='dashed', color='k', linewidth=1)
 
    #Set Tick Marks
    ax2.tick_params(axis='both', which='both', direction='in', top=True, bottom=True, left=True, right=True, labelsize=19)

    plt.subplots_adjust(wspace = 0.0000000001, bottom=0.13, top=0.98, left=0.085, right=0.98) #, hspace = 0.001, wspace = 0.001)
    line_labels=[r'Commissioning Experiment (Hall C)' '\n' r'(72 hrs at 45-60 $\mu$A)', 'Projected Statistical Errors ' '\n' r'(352 hrs at 70 $\mu$A, assuming 50\% beam eff.)']
    ax2.legend([l1, l0], line_labels, loc='upper right', frameon=False, fontsize=12)
    
    plt.show()
    '''
    
def main():
    print('Entering Main . . .')

    make_prl_plots(1)

if __name__=="__main__":
    main()


