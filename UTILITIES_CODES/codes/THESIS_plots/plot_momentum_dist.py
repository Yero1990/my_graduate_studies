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


matplotlib.use('Agg')
plt.ioff() # prevent plots froms displaying
#Use latex commands (e.g. \textit ot \textbf)
rc('text', usetex=True)
#Set default font to times new roman
font = {'family' : 'Times New Roman',
        'weight' : 'bold',
        'size'   : 12
}
plt.rc('font', **font)


#Set font
csfont = {'fontname':'Times New Roman'}

def convert2NaN(arr=np.array([]), value=0):
    #method to convert a specified value in a array to nan (not a number)
    
    for i in enumerate(arr):
        if arr[i[0]]==value:
            arr[i[0]] = np.nan
    return arr


#Conversion factor:  1 fm = 1/ (197 MeV),   
#The reduced cross section is in MeV^-3 units
MeV2fm = 197.3**3    #convert MeV^-3 to fm^3

#========User Input (Dir. Name to store output)=======
#usage: ipython plot_momentum_dist.py

dir_name1 = "thesis_redXsec_plots"
dir_name2 = "thesis_relError_plots"
dir_name3 = "thesis_kin_tables"

#check if directory exists, else creates it.
if not os.path.exists(dir_name1):
    os.makedirs(dir_name1)
if not os.path.exists(dir_name2):
    os.makedirs(dir_name2)
if not os.path.exists(dir_name3):
    os.makedirs(dir_name3)     

#Get Reduced Xsec Data File
fname_Q2_3to4 = '../deep_data_files/redXsec_combined_Q2_3to4.txt'   #Read Q2 = 3.5 +/- 0.5 GeV2 data (Hall C E12-10-003)
fname_Q2_4to5 = '../deep_data_files/redXsec_combined_Q2_4to5.txt'   #Read Q2 = 4.5 +/- 0.5 GeV2 data (Hall C E12-10-003)

f3to4 = B.get_file(fname_Q2_3to4)
f4to5 = B.get_file(fname_Q2_4to5)

    
#Get Bin Information (Same info for all files), regardless of kinematic Q2 bin (as long as histogram binning is samee)                                                                                              
i_b = np.array(B.get_data(f4to5, 'i_b'))    #2D bin number                                                                    
i_x = np.array(B.get_data(f4to5, 'i_x'))    #x (th_nq) bin number                                                                                    
i_y = np.array(B.get_data(f4to5, 'i_y'))    #y (pmiss) bin number                                        
thnq = np.array(B.get_data(f4to5, 'xb'))      #th_nq value at bin center
pm =  np.array(B.get_data(f4to5, 'yb'))      #pmiss value at bin center   

#Get averaged missing momentum
pm_avg_3to4 = np.array(B.get_data(f3to4, 'pm_avg'))
pm_avg_4to5 = np.array(B.get_data(f4to5, 'pm_avg'))

#Get Combined Final Red Xsec for all kinematics (Q2 = 3.5 +/- 0.5 GeV2)
red_dataXsec_avg_3to4          = np.array(B.get_data(f3to4,'red_dataXsec_avg'))
red_dataXsec_avg_err_3to4      = np.array(B.get_data(f3to4,'red_dataXsec_avg_err'))  #statistical error
red_dataXsec_avg_syst_err_3to4 = np.array(B.get_data(f3to4,'red_dataXsec_avg_syst_err')) 
red_dataXsec_avg_tot_err_3to4  = np.array(B.get_data(f3to4,'red_dataXsec_avg_tot_err'))
red_pwiaXsec_avg_3to4          = np.array(B.get_data(f3to4,'red_pwiaXsec_avg'))
red_fsiXsec_avg_3to4           = np.array(B.get_data(f3to4,'red_fsiXsec_avg'))

#Get Combined Final Red Xsec for all kinematics (Q2 = 4.5 +/- 0.5 GeV2)
red_dataXsec_avg_4to5          = np.array(B.get_data(f4to5,'red_dataXsec_avg'))
red_dataXsec_avg_err_4to5      = np.array(B.get_data(f4to5,'red_dataXsec_avg_err'))
red_dataXsec_avg_syst_err_4to5 = np.array(B.get_data(f4to5,'red_dataXsec_avg_syst_err'))
red_dataXsec_avg_tot_err_4to5  = np.array(B.get_data(f4to5,'red_dataXsec_avg_tot_err'))
red_pwiaXsec_avg_4to5          = np.array(B.get_data(f4to5,'red_pwiaXsec_avg'))
red_fsiXsec_avg_4to5           = np.array(B.get_data(f4to5,'red_fsiXsec_avg'))

#Get total relative errors / (Pm, thnq) bin for plotting (write to table)
kin_syst_tot_3to4 = np.array(dfile(fname_Q2_3to4)['kin_syst_tot'])    #kinematic systematics
norm_syst_tot_3to4 = np.array(dfile(fname_Q2_3to4)['norm_syst_tot'])  #normalization systematics
tot_syst_err_3to4 = np.array(dfile(fname_Q2_3to4)['tot_syst_err'])    #total systematics
tot_stats_err_3to4 = np.array(dfile(fname_Q2_3to4)['tot_stats_err'])  #total statistical
tot_err_3to4 = np.array(dfile(fname_Q2_3to4)['tot_err'])              #overall error

kin_syst_tot_4to5 = np.array(dfile(fname_Q2_4to5)['kin_syst_tot'])    #kinematic systematics
norm_syst_tot_4to5 = np.array(dfile(fname_Q2_4to5)['norm_syst_tot'])  #normalization systematics
tot_syst_err_4to5 = np.array(dfile(fname_Q2_4to5)['tot_syst_err'])    #total systematics
tot_stats_err_4to5 = np.array(dfile(fname_Q2_4to5)['tot_stats_err'])  #total statistical
tot_err_4to5 = np.array(dfile(fname_Q2_4to5)['tot_err'])              #overall error


def read_halla_data(thnq=0):

    #The code read data files containing momentum distributions from Hall A

    #Conversion factor: 1 fm = 1 / (0.1973 GeV)  or  1 fm^-1 = 0.1973 GeV
    #fminv2GeV = 1./(1./0.1973)
    ifm2GeV = 0.1973  #inverse fermi to GeV

    thnq_f = '%.1f' % (thnq)

    #fname = '../HallA_data/theta_%s_th_corr.data' % (thnq_f)  #older version Hall A data
    fname = './HallA_data/pm_distribution_results_q3_1_%s.data' % (thnq_f)  #Hall A data

    kin = dfile(fname)
    p_miss = np.array(kin['p_miss_av']) * ifm2GeV    #Pm here is not the central bin value, but an average over that bin
    red_dataXsec = kin['rho']             #units fm^3
    red_dataXsec_err = kin['delta_rho']

    return p_miss, red_dataXsec, red_dataXsec_err

def read_misak_pdist(model=""):

    #The code read data files containing misak's calculation of theory momentum distirbutions using Hall C kinematics

    #Conversion factor: 1 fm = 1 / (0.1973 GeV)  or  1 fm^-1 = 0.1973 GeV
    #fminv2GeV = 1./(1./0.1973)
    GeV2fm = 0.1973**3  #inverse fermi^3 to GeV^3

    if(model=="Paris"):
        fname = './HallA_data/mom_dist_paris_mevc.data'  
    elif(model=="AV18"):
        fname = './HallA_data/mom_dist_v18_mevc.data'  
    elif(model=="CD-Bonn"):
        fname = './HallA_data/mom_dist_cdbonn_mevc.data'  


    kin = dfile(fname)
    p_miss = np.array(kin['pm']) / 1000.            #convert momentum from  MeV to GeV
    red_dataXsec = np.array(kin['rho']) / GeV2fm           #Convert GeV^-3 to fm^3 units 

    #iNTERPOLATE
    f_red_misak = interp1d(p_miss, red_dataXsec, fill_value='extrapolate', kind='linear')  

    return p_miss, f_red_misak

def read_theoretical_models(theory="", model="", thnq=0):

    #This code read the averaged red.Xsec and returns arrays in pm and reduced Xsec
    #theory: V18, CD-Bonn    model: PWIA, FSI
    thnq_f = "%.2f" %(thnq)

    fname = './theoretical_models/updated_%s_models/theoryXsec_%s%s_thnq%s_combined.data' % (model, theory, model, thnq_f)
    kin = dfile(fname)
    
    pm_bin = np.array(kin['pm_avg'])
    

    #print(pm_bin)
    if(model=="PWIA" and theory=="V18"):
        red_pwiaXsec_V18 = kin['red_pwiaXsec_theory']
        return pm_bin, red_pwiaXsec_V18 
    if(model=="PWIA" and theory=="CD-Bonn"):                                                                                                                          
        red_pwiaXsec_CD_Bonn = kin['red_pwiaXsec_theory']                                                                          
        return pm_bin, red_pwiaXsec_CD_Bonn
    if(model=="FSI" and theory=="V18"):                                                 
        red_fsiXsec_V18 = kin['red_fsiXsec_theory']                                               
        return pm_bin, red_fsiXsec_V18 
    if(model=="FSI" and theory=="CD-Bonn"):                                                                                            
        red_fsiXsec_CD_Bonn = kin['red_fsiXsec_theory'] 
        return pm_bin, red_fsiXsec_CD_Bonn


def plot_momentum_dist():
    #Make THESIS paper plots

    #Read This Experiment (Hall C) Data, and require better than 50% statistics for both Q2 = 3.5 and Q2 = 4.5
    red_dataXsec_avg_masked_3to4 = np.ma.array(red_dataXsec_avg_3to4, mask=(red_dataXsec_avg_err_3to4>0.5*red_dataXsec_avg_3to4))    
    red_dataXsec_avg_masked_4to5 = np.ma.array(red_dataXsec_avg_4to5, mask=(red_dataXsec_avg_err_4to5>0.5*red_dataXsec_avg_4to5))

    #Apply masks to relative errors
    kin_syst_tot_3to4_m  = np.ma.masked_where(np.ma.getmask(red_dataXsec_avg_masked_3to4), kin_syst_tot_3to4) # applies the mask of m on x
    norm_syst_tot_3to4_m  = np.ma.masked_where(np.ma.getmask(red_dataXsec_avg_masked_3to4), norm_syst_tot_3to4) # applies the mask of m on x
    tot_syst_err_3to4_m  = np.ma.masked_where(np.ma.getmask(red_dataXsec_avg_masked_3to4), tot_syst_err_3to4)
    tot_stats_err_3to4_m  = np.ma.masked_where(np.ma.getmask(red_dataXsec_avg_masked_3to4), tot_stats_err_3to4)
    tot_err_3to4_m =  np.ma.masked_where(np.ma.getmask(red_dataXsec_avg_masked_3to4), tot_err_3to4)

    kin_syst_tot_4to5_m  = np.ma.masked_where(np.ma.getmask(red_dataXsec_avg_masked_4to5), kin_syst_tot_4to5) # applies the mask of m on x
    norm_syst_tot_4to5_m  = np.ma.masked_where(np.ma.getmask(red_dataXsec_avg_masked_4to5), norm_syst_tot_4to5) # applies the mask of m on x
    tot_syst_err_4to5_m  = np.ma.masked_where(np.ma.getmask(red_dataXsec_avg_masked_4to5), tot_syst_err_4to5)
    tot_stats_err_4to5_m  = np.ma.masked_where(np.ma.getmask(red_dataXsec_avg_masked_4to5), tot_stats_err_4to5)
    tot_err_4to5_m =  np.ma.masked_where(np.ma.getmask(red_dataXsec_avg_masked_4to5), tot_err_4to5)

    #convert masks to nan
    red_dataXsec_avg_masked_3to4 = np.ma.filled(red_dataXsec_avg_masked_3to4.astype(float), np.nan)
    red_dataXsec_avg_masked_4to5 = np.ma.filled(red_dataXsec_avg_masked_4to5.astype(float), np.nan)

    #convert rel. err masks to nan
    kin_syst_tot_3to4_m = np.ma.filled(kin_syst_tot_3to4_m.astype(float), np.nan)
    norm_syst_tot_3to4_m = np.ma.filled(norm_syst_tot_3to4_m.astype(float), np.nan)
    tot_syst_err_3to4_m = np.ma.filled(tot_syst_err_3to4_m.astype(float), np.nan)
    tot_stats_err_3to4_m = np.ma.filled(tot_stats_err_3to4_m.astype(float), np.nan)
    tot_err_3to4_m = np.ma.filled(tot_err_3to4_m.astype(float), np.nan)

    kin_syst_tot_4to5_m = np.ma.filled(kin_syst_tot_4to5_m.astype(float), np.nan)
    norm_syst_tot_4to5_m = np.ma.filled(norm_syst_tot_4to5_m.astype(float), np.nan)
    tot_syst_err_4to5_m = np.ma.filled(tot_syst_err_4to5_m.astype(float), np.nan)
    tot_stats_err_4to5_m = np.ma.filled(tot_stats_err_4to5_m.astype(float), np.nan)
    tot_err_4to5_m = np.ma.filled(tot_err_4to5_m.astype(float), np.nan)

    #Read Hall A data (reduced Xsec already in fm^3 units, and Precoil in GeV)
    pm_ha35, red_dataXsec_ha35, red_dataXsec_err_ha35 = read_halla_data(35)
    pm_ha45, red_dataXsec_ha45, red_dataXsec_err_ha45 = read_halla_data(45)  
    pm_ha75, red_dataXsec_ha75, red_dataXsec_err_ha75 = read_halla_data(75)  


    #TO-DO: Mask certain values of theoretical cross section for which there is no data
    #(so they don't get plotted at thnq=5, 15, 25 deg, etc. extreme angles)

    #Plot momentum distribution vs. Pmiss for different theta_nq
    thnq_arr = [5, 15, 25, 35, 45, 55, 65, 75, 85, 95, 105]

    for i, ithnq in enumerate(thnq_arr):

        th_nq_min = ithnq - 5
        th_nq_max = ithnq + 5

        
        #Read Other Theoretical Models (V18, CD-BONN) (Only at Q2 = 4.5 +/- 0.5)
        pm_avg1, red_pwiaXsec_V18 = read_theoretical_models("V18", "PWIA", ithnq)
        pm_avg2, red_fsiXsec_V18 = read_theoretical_models("V18", "FSI", ithnq) 
        pm_avg3, red_pwiaXsec_CD = read_theoretical_models("CD-Bonn", "PWIA", ithnq)                                   
        pm_avg4, red_fsiXsec_CD = read_theoretical_models("CD-Bonn", "FSI", ithnq)    


        #Read Momentum Dsitributions from Misak (Alternative to calculateing reduced Cross section)
        pm_paris, redXsec_Paris_MS = read_misak_pdist(model="Paris")
        pm_av18, redXsec_AV18_MS = read_misak_pdist(model="AV18")
        pm_cd, redXsec_CD_MS = read_misak_pdist(model="CD-Bonn")

        
        #------Convert MeV^-3 to fm^3 for theory--------

        #LAGET THEORY
        #red_pwiaXsec_avg_3to4[thnq==ithnq] = red_pwiaXsec_avg_3to4[thnq==ithnq]*MeV2fm    #Laget PWIA
        #red_fsiXsec_avg_3to4[thnq==ithnq] = red_fsiXsec_avg_3to4[thnq==ithnq]*MeV2fm    #Laget FSI

        red_pwiaXsec_avg_4to5[thnq==ithnq] = red_pwiaXsec_avg_4to5[thnq==ithnq]*MeV2fm    #Laget PWIA
        red_fsiXsec_avg_4to5[thnq==ithnq] = red_fsiXsec_avg_4to5[thnq==ithnq]*MeV2fm    #Laget FSI
              

        red_pwiaXsec_V18 = np.array(red_pwiaXsec_V18)*MeV2fm
        red_fsiXsec_V18 = np.array(red_fsiXsec_V18)*MeV2fm
        red_pwiaXsec_CD = np.array(red_pwiaXsec_CD)*MeV2fm
        red_fsiXsec_CD = np.array(red_fsiXsec_CD)*MeV2fm

        #Paris LAGET (Due to limited data, adjustements on interpolation have to be made,depending on the data)
        #thnq: 5, 15, 95 deg (plot with pm_avg_4to5[thnq])   configuration1, c1 (interpolation)
        f_red_pwiaXsec_avg_4to5_c1 = interp1d(pm_avg_4to5[thnq==ithnq], red_pwiaXsec_avg_4to5[thnq==ithnq],fill_value='extrapolate', kind='linear') 
        f_red_fsiXsec_avg_4to5_c1 = interp1d(pm_avg_4to5[thnq==ithnq], red_fsiXsec_avg_4to5[thnq==ithnq],fill_value='extrapolate', kind='linear') 

        #thnq: 35, 45 deg (plot with pm[thnq])
        f_red_pwiaXsec_avg_4to5_c2 = interp1d(pm[thnq==ithnq], red_pwiaXsec_avg_4to5[thnq==ithnq],fill_value='extrapolate', kind='cubic') 
        f_red_fsiXsec_avg_4to5_c2 = interp1d(pm[thnq==ithnq], red_fsiXsec_avg_4to5[thnq==ithnq],fill_value='extrapolate', kind='cubic') 

        #thnq:  25, 55, 65, 75, 85 (plot with pm[thnq])
        f_red_pwiaXsec_avg_4to5_c3 = interp1d(pm[thnq==ithnq], red_pwiaXsec_avg_4to5[thnq==ithnq],fill_value='extrapolate', kind='linear') 
        f_red_fsiXsec_avg_4to5_c3 = interp1d(pm[thnq==ithnq], red_fsiXsec_avg_4to5[thnq==ithnq],fill_value='extrapolate', kind='linear') 

        
        
        f_red_pwiaXsec_V18 = interp1d(pm_avg1, red_pwiaXsec_V18,fill_value='extrapolate', kind='cubic')                        #AV18 (M. Sargsian calculation)
        f_red_fsiXsec_V18 = interp1d(pm_avg2, red_fsiXsec_V18,fill_value='extrapolate', kind='cubic')
       
        f_red_pwiaXsec_CD = interp1d(pm_avg3, red_pwiaXsec_CD,fill_value='extrapolate', kind='cubic')                          #CD-Bonn (M. Sargsian calculation)
        f_red_fsiXsec_CD = interp1d(pm_avg4, red_fsiXsec_CD,fill_value='extrapolate', kind='cubic')

        
        print(ithnq)

        '''
        #----------PLOT FINAL REDUCED XSEC VS. THETA_nq------------

        B.pl.clf()
        B.pl.figure(i)


        B.plot_exp(pm_avg[thnq==ithnq], red_dataXsec_avg_masked[thnq==ithnq]*MeV2fm, red_dataXsec_avg_tot_err[thnq==ithnq]*MeV2fm, marker='o', markersize=6, color='k', markerfacecolor='k', label='This Experiment (Hall C)', zorder=3)
                
        #--------------HALL A DATA---------------
        if(ithnq==35):
            #Plot Hall A Data for thnq==35 or 45 deg
            B.plot_exp(pm_ha35, red_dataXsec_ha35, red_dataXsec_err_ha35, marker='s', markersize=6, mfc='white', color='r', label='Hall A Data', zorder=2)
        if(ithnq==45):                                                                                                             
            #Plot Hall A Data for thnq==35 or 45 deg                                                                                                                                                    
            B.plot_exp(pm_ha45, red_dataXsec_ha45, red_dataXsec_err_ha45, marker='s', markersize=6, mfc='white', color='r', label='Hall A Data', zorder=2)
        if(ithnq==75):                                                                                                             
            #Plot Hall A Data for thnq==35 or 45 deg                                                                                                                                                    
            B.plot_exp(pm_ha75, red_dataXsec_ha75, red_dataXsec_err_ha75, marker='s', markersize=6, mfc='white', color='r', label='Hall A Data', zorder=2)
        #----------------------------------------

        
        #Plot theoretical curves
        B.plot_exp(pm_avg[thnq==ithnq], f_red_pwiaXsec_avg(pm_avg[thnq==ithnq]), linestyle='--', marker='None', color='blue', logy=True, label='Paris PWIA')
        B.plot_exp(pm_avg[thnq==ithnq], f_red_fsiXsec_avg(pm_avg[thnq==ithnq]), linestyle='-', marker='None', color='blue', logy=True, label='Paris FSI')
       
        #Q2 = (3,4) GeV2
        #B.plot_exp(pm_avg0[thnq==ithnq], f_red_pwiaXsec_avg0(pm_avg0[thnq==ithnq]), linestyle='--', marker='None', color='red', logy=True, label=r'Paris PWIA ($Q^{2}=3.5\pm0.5$ GeV$^{2}$)')
        #B.plot_exp(pm_avg0[thnq==ithnq], f_red_fsiXsec_avg0(pm_avg0[thnq==ithnq]), linestyle='-', marker='None', color='red', logy=True, label=r'Paris FSI ($Q^{2}=3.5\pm0.5$ GeV$^{2}$)')


        B.plot_exp(pm_avg1, f_red_pwiaXsec_V18(pm_avg1), linestyle='--', marker='None', color='green', logy=True, label='AV18 PWIA')   
        B.plot_exp(pm_avg2, f_red_fsiXsec_V18(pm_avg2), linestyle='-', marker='None', color='green', logy=True, label='AV18 FSI') 

        B.plot_exp(pm_avg3, f_red_pwiaXsec_CD(pm_avg3), linestyle='--', marker='None', color='magenta', logy=True, label='CD-Bonn PWIA')     
        B.plot_exp(pm_avg4, f_red_fsiXsec_CD(pm_avg4), linestyle='-', marker='None', color='magenta', logy=True, label='CD-Bonn FSI') 


        B.pl.xlabel('')
        B.pl.ylabel(r'$\sigma_{red} (fm^{3})$')
        B.pl.xlim(0., 1.2)
        B.pl.title(r'Reduced Cross Section, $\theta_{nq} = %i \pm 5$ deg '%(ithnq))
        B.pl.legend(frameon=False, fontsize=12, loc='upper right')            

        #ax2 = B.pl.subplot(212, sharex = ax1)
        B.pl.xlabel('$p_{r}$ (GeV/c)')

        B.pl.show()
        
        #--------------------------------Calculate and PLOT RATIO-----------------------------
 
        #--------Calculate Ratio of redXsec relative to CD-Bonn
        #Calculate ref and data ratios
        R_ref = f_red_pwiaXsec_CD(pm_avg[thnq==ithnq]) / f_red_pwiaXsec_CD(pm_avg[thnq==ithnq])
        R_data = red_dataXsec_avg_masked[thnq==ithnq]*MeV2fm / f_red_pwiaXsec_CD(pm_avg[thnq==ithnq])
        R_data_err = red_dataXsec_avg_tot_err[thnq==ithnq]*MeV2fm / f_red_pwiaXsec_CD(pm_avg[thnq==ithnq])

        #Theory curves ratio to CD-BONN
        R_CDBonn_fsi = f_red_fsiXsec_CD(pm_avg[thnq==ithnq]) / f_red_pwiaXsec_CD(pm_avg[thnq==ithnq])
        R_Paris_pwia = f_red_pwiaXsec_avg(pm_avg[thnq==ithnq]) / f_red_pwiaXsec_CD(pm_avg[thnq==ithnq])
        R_Paris_fsi = f_red_fsiXsec_avg(pm_avg[thnq==ithnq]) / f_red_pwiaXsec_CD(pm_avg[thnq==ithnq])
        R_AV18_pwia = f_red_pwiaXsec_V18(pm_avg[thnq==ithnq]) / f_red_pwiaXsec_CD(pm_avg[thnq==ithnq])
        R_AV18_fsi = f_red_fsiXsec_V18(pm_avg[thnq==ithnq]) / f_red_pwiaXsec_CD(pm_avg[thnq==ithnq])

        #Plot the Data (and all models) to CD-Bonn PWIA model
        B.plot_exp(pm_avg[thnq==ithnq], R_data, R_data_err, marker='o', color='k', label='Hall C Data', capsize=0)
        B.plot_exp(pm_avg[thnq==ithnq], R_ref,  marker='None', linestyle='--', color='#ff00ff', label='CD-Bonn PWIA')
        B.plot_exp(pm_avg[thnq==ithnq], R_CDBonn_fsi, marker='None', linestyle='-', color='#ff00ff', label='CD-Bonn FSI')
        B.plot_exp(pm_avg[thnq==ithnq], R_Paris_pwia, marker='None', linestyle='--', color='b', label='Paris PWIA')
        B.plot_exp(pm_avg[thnq==ithnq], R_Paris_fsi, marker='None', linestyle='-', color='b', label='Paris FSI')
        B.plot_exp(pm_avg[thnq==ithnq], R_AV18_pwia, marker='None', linestyle='--', color='g', label='AV18 PWIA')
        B.plot_exp(pm_avg[thnq==ithnq], R_AV18_fsi, marker='None', linestyle='-', color='g', label='AV18 FSI')
        print(R_CDBonn_fsi)
        
        if ithnq==35:
            pm_ha, red_dataXsec_ha, red_dataXsec_err_ha = read_halla_data(ithnq)
            R_HAdata = red_dataXsec_ha / f_red_pwiaXsec_CD(pm_ha)                #Hall A data (Q2=3.5) ratio to CD-Bonn PWIA (Q2=4.5)
            R_HAdata_err =  red_dataXsec_err_ha / f_red_pwiaXsec_CD(pm_ha)
            B.plot_exp(pm_ha, R_HAdata, R_HAdata_err, marker='s', color='r', mfc='white',  label='Hall A Data', capsize=0)
        if ithnq==45:
            pm_ha, red_dataXsec_ha, red_dataXsec_err_ha = read_halla_data(ithnq)
            R_HAdata = red_dataXsec_ha / f_red_pwiaXsec_CD(pm_ha)                #Hall A data (Q2=3.5) ratio to CD-Bonn PWIA (Q2=4.5)
            R_HAdata_err =  red_dataXsec_err_ha / f_red_pwiaXsec_CD(pm_ha) 
            B.plot_exp(pm_ha, R_HAdata, R_HAdata_err, marker='s', color='r', mfc='white',  label='Hall A Data', capsize=0)
        if ithnq==75:
            pm_ha, red_dataXsec_ha, red_dataXsec_err_ha = read_halla_data(ithnq)
            R_HAdata = red_dataXsec_ha / f_red_pwiaXsec_CD(pm_ha)                #Hall A data (Q2=3.5) ratio to CD-Bonn PWIA (Q2=4.5)
            R_HAdata_err =  red_dataXsec_err_ha / f_red_pwiaXsec_CD(pm_ha)
            B.plot_exp(pm_ha, R_HAdata, R_HAdata_err, marker='s', color='r', mfc='white',  label='Hall A Data', capsize=0)
        
        
        #Set axis limits
        B.pl.xlim(0.0, 1.2)
        #B.pl.ylim(-0.05, 22.8)

        #Set Axes Labels for subplot 0
        B.pl.xlabel(r'$p_{r} [GeV/c]$')                                                                                                                                                                                                   
        B.pl.ylabel(r'$R=\sigma_{\mathrm{red}}/\sigma^{\mathrm{CD-Bonn PWIA}}_{\mathrm{red}}$')                                                                                                                                                                        
        B.pl.title(r'Reduced Cross Section Ratio: $\theta_{nq} = %i \pm 5$ deg '%(ithnq))
        B.pl.legend()
        B.pl.show()
        '''

        #-----------------------------ALTERNATIVE: MAKE SUBPLOTS----------------------------------
        B.pl.clf()
        #-----Create Subplots-----
        fig = B.pl.subplots(2, sharex=True, sharey=False, figsize=(7, 10)) 
        gs = gridspec.GridSpec(2, 1) 
        ax0 = B.pl.subplot(gs[0])
        ax0.tick_params(axis='both', which='both', direction='in', top=True, bottom=True, left=True, right=True, labelsize=19)    

        #Remove un-necessary X tick labels from subplot
        B.pl.setp(ax0.get_xticklabels(), visible=False)
        
        #----------PLOT HALL C DATA: FINAL REDUCED XSEC VS. THETA_nq------------
        B.plot_exp(pm_avg_4to5[(thnq==ithnq)], red_dataXsec_avg_masked_4to5[thnq==ithnq]*MeV2fm, red_dataXsec_avg_tot_err_4to5[thnq==ithnq]*MeV2fm, marker='o', markersize=6, c='k', mec='k', mfc='k', label=r'$Q^{2}=4.5\pm0.5$ GeV$^{2}$ (Hall C)', zorder=4,capsize=0)
        B.plot_exp(pm_avg_3to4[(thnq==ithnq)], red_dataXsec_avg_masked_3to4[(thnq==ithnq)]*MeV2fm, red_dataXsec_avg_tot_err_3to4[(thnq==ithnq)]*MeV2fm, marker='o', markersize=6, c='c', mec='c', mfc='white', label=r'$Q^{2}=3.5\pm0.5$ GeV$^{2}$ (Hall C)', zorder=3,capsize=0)

        #---------PLOT HALL A DATA---------------
        if(ithnq==35):
            #Plot Hall A Data for thnq==35 or 45 deg
            B.plot_exp(pm_ha35, red_dataXsec_ha35, red_dataXsec_err_ha35, marker='s', markersize=6, mfc='white', color='r', label='Hall A Data', zorder=2, capsize=0)
        if(ithnq==45):                                                                                                             
            #Plot Hall A Data for thnq==35 or 45 deg                                                                                                                                                    
            B.plot_exp(pm_ha45, red_dataXsec_ha45, red_dataXsec_err_ha45, marker='s', markersize=6, mfc='white', color='r', label='Hall A Data', zorder=2, capsize=0)
        if(ithnq==75):                                                                                                             
            #Plot Hall A Data for thnq==35 or 45 deg                                                                                                                                                    
            B.plot_exp(pm_ha75, red_dataXsec_ha75, red_dataXsec_err_ha75, marker='s', markersize=6, mfc='white', color='r', label='Hall A Data', zorder=2, capsize=0)
        #----------------------------------------

        
        #----PLOT THEORETICAL CURVES-----

        #---Paris (Laget calculation)  Q2 = 4.5 +/- 0.5 (multiple interpolations configurations depending on angle setting)
        if(ithnq==5 or ithnq==15 or ithnq==95):
            B.plot_exp(pm_avg_4to5[thnq==ithnq], f_red_pwiaXsec_avg_4to5_c1(pm_avg_4to5[thnq==ithnq]), linestyle='--', marker='None', color='blue', logy=True, label='Paris PWIA')
            B.plot_exp(pm_avg_4to5[thnq==ithnq], f_red_fsiXsec_avg_4to5_c1(pm_avg_4to5[thnq==ithnq]), linestyle='-', marker='None', color='blue', logy=True, label='Paris FSI')
        elif(ithnq==35 or ithnq==45):
            B.plot_exp(pm[thnq==ithnq], f_red_pwiaXsec_avg_4to5_c2(pm[thnq==ithnq]), linestyle='--', marker='None', color='blue', logy=True, label='Paris PWIA')
            B.plot_exp(pm[thnq==ithnq], f_red_fsiXsec_avg_4to5_c2(pm[thnq==ithnq]), linestyle='-', marker='None', color='blue', logy=True, label='Paris FSI')
        else:
            B.plot_exp(pm[thnq==ithnq], f_red_pwiaXsec_avg_4to5_c3(pm[thnq==ithnq]), linestyle='--', marker='None', color='blue', logy=True, label='Paris PWIA')
            B.plot_exp(pm[thnq==ithnq], f_red_fsiXsec_avg_4to5_c3(pm[thnq==ithnq]), linestyle='-', marker='None', color='blue', logy=True, label='Paris FSI')
        
        #AV18
        B.plot_exp(pm_avg1, f_red_pwiaXsec_V18(pm_avg1), linestyle='--', marker='None', color='green', logy=True, label='AV18 PWIA')   
        B.plot_exp(pm_avg2, f_red_fsiXsec_V18(pm_avg2), linestyle='-', marker='None', color='green', logy=True, label='AV18 FSI') 

        #CD-Bonn
        B.plot_exp(pm_avg3, f_red_pwiaXsec_CD(pm_avg3), linestyle='--', marker='None', color='magenta', logy=True, label='CD-Bonn PWIA')     
        B.plot_exp(pm_avg4, f_red_fsiXsec_CD(pm_avg4), linestyle='-', marker='None', color='magenta', logy=True, label='CD-Bonn FSI') 

        B.pl.xlabel('')
        B.pl.ylabel(r'$\sigma_{\mathrm{red}} [\mathrm{fm}^{3}]$', fontsize=18)
        B.pl.xlim(0., 1.2)
        B.pl.title(r'Reduced Cross Section, $\theta_{nq} = %i \pm 5^{\circ}$ '%(ithnq), fontsize=20)
        B.pl.legend(frameon=False, fontsize=12, loc='upper right')            

        
        #--------------------------------Calculate and PLOT RATIO-----------------------------
        #2nd subplot (ratio)

        ax1 = plt.subplot(gs[1], sharex=ax0)
        ax1.tick_params(axis='both', which='both', direction='in', top=True, bottom=True, left=True, right=True, labelsize=19)    

        
        #--------Calculate Ratio of redXsec relative to CD-Bonn-----
        #Reference
        R_ref = f_red_pwiaXsec_CD(pm[thnq==ithnq]) / f_red_pwiaXsec_CD(pm[thnq==ithnq])

        #DATA / CD-BONN PWIA
        R_data_4to5 = red_dataXsec_avg_masked_4to5[thnq==ithnq]*MeV2fm / f_red_pwiaXsec_CD(pm_avg_4to5[thnq==ithnq])
        R_data_err_4to5 = red_dataXsec_avg_tot_err_4to5[thnq==ithnq]*MeV2fm / f_red_pwiaXsec_CD(pm_avg_4to5[thnq==ithnq])
        
        R_data_3to4 = red_dataXsec_avg_masked_3to4[thnq==ithnq]*MeV2fm / f_red_pwiaXsec_CD(pm_avg_3to4[thnq==ithnq])
        R_data_err_3to4 = red_dataXsec_avg_tot_err_3to4[thnq==ithnq]*MeV2fm / f_red_pwiaXsec_CD(pm_avg_3to4[thnq==ithnq])
        
        #Theory curves ratio to CD-BONN

        #Paris LAGET            
        #thnq: 5, 15, 95 deg
        R_Paris_pwia_4to5_c1 = f_red_pwiaXsec_avg_4to5_c1(pm_avg_4to5[thnq==ithnq]) / f_red_pwiaXsec_CD(pm_avg_4to5[thnq==ithnq])
        R_Paris_fsi_4to5_c1 = f_red_fsiXsec_avg_4to5_c1(pm_avg_4to5[thnq==ithnq]) / f_red_pwiaXsec_CD(pm_avg_4to5[thnq==ithnq])
        #thnq: 35, 45 deg
        R_Paris_pwia_4to5_c2 = f_red_pwiaXsec_avg_4to5_c2(pm[thnq==ithnq]) / f_red_pwiaXsec_CD(pm[thnq==ithnq])
        R_Paris_fsi_4to5_c2 = f_red_fsiXsec_avg_4to5_c2(pm[thnq==ithnq]) / f_red_pwiaXsec_CD(pm[thnq==ithnq])
        #thnq: 25, 55, 65, 75, 85
        R_Paris_pwia_4to5_c3 = f_red_pwiaXsec_avg_4to5_c3(pm[thnq==ithnq]) / f_red_pwiaXsec_CD(pm[thnq==ithnq])
        R_Paris_fsi_4to5_c3 = f_red_fsiXsec_avg_4to5_c3(pm[thnq==ithnq]) / f_red_pwiaXsec_CD(pm[thnq==ithnq])

        #AV 18
        R_AV18_pwia = f_red_pwiaXsec_V18(pm_avg1) / f_red_pwiaXsec_CD(pm_avg1)
        R_AV18_fsi = f_red_fsiXsec_V18(pm_avg2) / f_red_pwiaXsec_CD(pm_avg2)
        #CD-BONN
        R_CDBonn_fsi = f_red_fsiXsec_CD(pm_avg4) / f_red_pwiaXsec_CD(pm_avg4)

        
        #----Plot the RATIO OF Data to CD-Bonn PWIA model
        B.plot_exp(pm_avg_4to5[thnq==ithnq], R_data_4to5, R_data_err_4to5, marker='o', markersize=6, c='k', mec='k', mfc='k',     label=r'$Q^{2}=4.5\pm0.5$ GeV$^{2}$ (Hall C)', capsize=0, zorder=4)
        B.plot_exp(pm_avg_3to4[thnq==ithnq], R_data_3to4, R_data_err_3to4, marker='o', markersize=6, c='c', mec='c', mfc='white', label=r'$Q^{2}=3.5\pm0.5$ GeV$^{2}$ (Hall C)', capsize=0, zorder=3)

        #-----PLOT THEORETICAL RATIO-----
        #Paris LAGET
        if(ithnq==5 or ithnq==15 or ithnq==95):
            B.plot_exp(pm_avg_4to5[thnq==ithnq], R_Paris_pwia_4to5_c1, marker='None', linestyle='--', color='b', label='Paris PWIA')
            B.plot_exp(pm_avg_4to5[thnq==ithnq], R_Paris_fsi_4to5_c1, marker='None', linestyle='-', color='b', label='Paris FSI')
        if(ithnq==35 or ithnq==45):
            B.plot_exp(pm[thnq==ithnq], R_Paris_pwia_4to5_c2, marker='None', linestyle='--', color='b', label='Paris PWIA')
            B.plot_exp(pm[thnq==ithnq], R_Paris_fsi_4to5_c2, marker='None', linestyle='-', color='b', label='Paris FSI')
        else:
            B.plot_exp(pm[thnq==ithnq], R_Paris_pwia_4to5_c3, marker='None', linestyle='--', color='b', label='Paris PWIA')
            B.plot_exp(pm[thnq==ithnq], R_Paris_fsi_4to5_c3, marker='None', linestyle='-', color='b', label='Paris FSI')

        #Plot CD-BONN
        B.plot_exp(pm[thnq==ithnq], R_ref,  marker='None', linestyle='--', color='#ff00ff', label='CD-Bonn PWIA')
        B.plot_exp(pm_avg4, R_CDBonn_fsi, marker='None', linestyle='-', color='#ff00ff', label='CD-Bonn FSI')
        #Plot AV 18
        B.plot_exp(pm_avg1, R_AV18_pwia, marker='None', linestyle='--', color='g', label='AV18 PWIA')
        B.plot_exp(pm_avg2, R_AV18_fsi, marker='None', linestyle='-', color='g', label='AV18 FSI')
        
        if ithnq==35:
            pm_ha, red_dataXsec_ha, red_dataXsec_err_ha = read_halla_data(ithnq)
            R_HAdata = red_dataXsec_ha / f_red_pwiaXsec_CD(pm_ha)                #Hall A data (Q2=3.5) ratio to CD-Bonn PWIA (Q2=4.5)
            R_HAdata_err =  red_dataXsec_err_ha / f_red_pwiaXsec_CD(pm_ha)
            B.plot_exp(pm_ha, R_HAdata, R_HAdata_err, marker='s', color='r', mfc='white',  label='$Q^{2}=3.5\pm0.25$ GeV$^{2}$ (Hall A)', capsize=0, zorder=2)
        if ithnq==45:
            pm_ha, red_dataXsec_ha, red_dataXsec_err_ha = read_halla_data(ithnq)
            R_HAdata = red_dataXsec_ha / f_red_pwiaXsec_CD(pm_ha)                #Hall A data (Q2=3.5) ratio to CD-Bonn PWIA (Q2=4.5)
            R_HAdata_err =  red_dataXsec_err_ha / f_red_pwiaXsec_CD(pm_ha) 
            B.plot_exp(pm_ha, R_HAdata, R_HAdata_err, marker='s', color='r', mfc='white',  label='$Q^{2}=3.5\pm0.25$ GeV$^{2}$ (Hall A)', capsize=0, zorder=2)
        if ithnq==75:
            pm_ha, red_dataXsec_ha, red_dataXsec_err_ha = read_halla_data(ithnq)
            R_HAdata = red_dataXsec_ha / f_red_pwiaXsec_CD(pm_ha)                #Hall A data (Q2=3.5) ratio to CD-Bonn PWIA (Q2=4.5)
            R_HAdata_err =  red_dataXsec_err_ha / f_red_pwiaXsec_CD(pm_ha)
            B.plot_exp(pm_ha, R_HAdata, R_HAdata_err, marker='s', color='r', mfc='white',  label='$Q^{2}=3.5\pm0.25$ GeV$^{2}$ (Hall A)', capsize=0, zorder=2)
        
        
        #Set axis limits
        B.pl.xlim(0.0, 1.2)
        if(ithnq==5 or ithnq==15):
            B.pl.ylim(-1., 10.5)
        if(ithnq==15):
            B.pl.ylim(-0.5, 12.3)
        elif(ithnq==25):
            B.pl.ylim(-1., 20.5)
        elif(ithnq==35):
            B.pl.ylim(-1., 29.9)            
        elif(ithnq==45 or ithnq==55):
            B.pl.ylim(-5., 119.9)
        elif(ithnq==65):
            B.pl.ylim(-1.5, 66.5)
        elif(ithnq==75):
            B.pl.ylim(-1., 27.5)
        elif(ithnq==85):
            B.pl.ylim(0.2, 13.5)
        elif(ithnq==95):
            B.pl.ylim(0.2, 10.5)
        else:
            print('OK')

        #Remove spacing between subplots
        plt.subplots_adjust(hspace = 0.01, bottom=0.08, top=0.95, right=0.95, left=0.2) #, hspace = 0.001, wspace = 0.001)

        #Set Axes Labels for subplot 0
        B.pl.xlabel(r'$p_{\mathrm{r}} $ [GeV/c]', fontsize=18)                                                                                                                                                                                                   
        B.pl.ylabel(r'$R=\sigma_{\mathrm{red}}/\sigma^{\mathrm{CD-Bonn PWIA}}_{\mathrm{red}}$', fontsize=18)                                                                                                                                                                        
        B.pl.title('')

        
        #CREATE INSET PLOT
        if(ithnq==5 or ithnq==25 or ithnq==95 or ithnq==75 or ithnq==85):
            axins = inset_axes(ax1, width='60%', height='60%', loc='upper right')
            axins.yaxis.tick_left()
        elif(ithnq==45 or ithnq==55 or ithnq==35 or ithnq==25):
            axins = inset_axes(ax1, width='70%', height='70%', loc='upper left')
            axins.yaxis.tick_right()            
        else:
            axins = inset_axes(ax1, width='60%', height='60%', loc='upper left')
            axins.yaxis.tick_right()
       
        #----Plot the RATIO OF Data to CD-Bonn PWIA model
        B.plot_exp(pm_avg_4to5[thnq==ithnq], R_data_4to5, R_data_err_4to5, marker='o', markersize=6, c='k', mec='k', mfc='k',     label=r'', capsize=0, zorder=4)
        B.plot_exp(pm_avg_3to4[thnq==ithnq], R_data_3to4, R_data_err_3to4, marker='o', markersize=6, c='c', mec='c', mfc='white', label=r'', capsize=0, zorder=3)

        
         #-----PLOT THEORETICAL RATIO-----
        #Paris LAGET
        if(ithnq==5 or ithnq==15 or ithnq==95):
            B.plot_exp(pm_avg_4to5[thnq==ithnq], R_Paris_pwia_4to5_c1, marker='None', linestyle='--', color='b', label='')
            B.plot_exp(pm_avg_4to5[thnq==ithnq], R_Paris_fsi_4to5_c1, marker='None', linestyle='-', color='b', label='')
        if(ithnq==35 or ithnq==45):
            B.plot_exp(pm[thnq==ithnq], R_Paris_pwia_4to5_c2, marker='None', linestyle='--', color='b', label='')
            B.plot_exp(pm[thnq==ithnq], R_Paris_fsi_4to5_c2, marker='None', linestyle='-', color='b', label='Paris FSI')
        else:
            B.plot_exp(pm[thnq==ithnq], R_Paris_pwia_4to5_c3, marker='None', linestyle='--', color='b', label='Paris PWIA')
            B.plot_exp(pm[thnq==ithnq], R_Paris_fsi_4to5_c3, marker='None', linestyle='-', color='b', label='Paris FSI')

        
        #Plot CD-BONN
        B.plot_exp(pm[thnq==ithnq], R_ref,  marker='None', linestyle='--', color='#ff00ff', label='CD-Bonn PWIA')
        B.plot_exp(pm_avg4, R_CDBonn_fsi, marker='None', linestyle='-', color='#ff00ff', label='CD-Bonn FSI')
        #Plot AV 18
        B.plot_exp(pm_avg1, R_AV18_pwia, marker='None', linestyle='--', color='g', label='AV18 PWIA')
        B.plot_exp(pm_avg2, R_AV18_fsi, marker='None', linestyle='-', color='g', label='AV18 FSI')
        
        if ithnq==35:
            pm_ha, red_dataXsec_ha, red_dataXsec_err_ha = read_halla_data(ithnq)
            R_HAdata = red_dataXsec_ha / f_red_pwiaXsec_CD(pm_ha)                #Hall A data (Q2=3.5) ratio to CD-Bonn PWIA (Q2=4.5)
            R_HAdata_err =  red_dataXsec_err_ha / f_red_pwiaXsec_CD(pm_ha)
            B.plot_exp(pm_ha, R_HAdata, R_HAdata_err, marker='s', color='r', mfc='white',  label='$Q^{2}=3.5\pm0.25$ GeV$^{2}$ (Hall A)', capsize=0, zorder=2)
        if ithnq==45:
            pm_ha, red_dataXsec_ha, red_dataXsec_err_ha = read_halla_data(ithnq)
            R_HAdata = red_dataXsec_ha / f_red_pwiaXsec_CD(pm_ha)                #Hall A data (Q2=3.5) ratio to CD-Bonn PWIA (Q2=4.5)
            R_HAdata_err =  red_dataXsec_err_ha / f_red_pwiaXsec_CD(pm_ha) 
            B.plot_exp(pm_ha, R_HAdata, R_HAdata_err, marker='s', color='r', mfc='white',  label='$Q^{2}=3.5\pm0.25$ GeV$^{2}$ (Hall A)', capsize=0, zorder=2)
        if ithnq==75:
            pm_ha, red_dataXsec_ha, red_dataXsec_err_ha = read_halla_data(ithnq)
            R_HAdata = red_dataXsec_ha / f_red_pwiaXsec_CD(pm_ha)                #Hall A data (Q2=3.5) ratio to CD-Bonn PWIA (Q2=4.5)
            R_HAdata_err =  red_dataXsec_err_ha / f_red_pwiaXsec_CD(pm_ha)
            B.plot_exp(pm_ha, R_HAdata, R_HAdata_err, marker='s', color='r', mfc='white',  label='$Q^{2}=3.5\pm0.25$ GeV$^{2}$ (Hall A)', capsize=0, zorder=2)



        if(ithnq==5):
            B.pl.xlim(-0.01, 0.195); B.pl.ylim(0.5, 1.45)
        elif(ithnq==15):
            B.pl.xlim(-0.01, 0.23); B.pl.ylim(0.7, 1.63)
        elif(ithnq==25):
            B.pl.xlim(-0.01, 0.823); B.pl.ylim(0.2, 2.7)
        elif(ithnq==35):
            B.pl.xlim(-0.01, 1.2); B.pl.ylim(0.3, 2.7) 
        elif(ithnq==45):
            B.pl.xlim(-0.01, 1.2); B.pl.ylim(0.5, 4.7)
        elif(ithnq==55):
            B.pl.xlim(-0.01, 0.823); B.pl.ylim(0.5, 7.7)            
        elif(ithnq==65):
            B.pl.xlim(-0.01, 0.73); B.pl.ylim(0.2, 10.7)
        elif(ithnq==75):
            B.pl.xlim(-0.01, 1.2); B.pl.ylim(0.3, 6.7)
        elif(ithnq==85):
            B.pl.xlim(-0.01, 0.47); B.pl.ylim(0.5, 3.7)
        elif(ithnq==95):
            B.pl.xlim(-0.01, 0.33); B.pl.ylim(0.2, 2.)
        else:
            print('OK')

            
        B.pl.xlabel('')
        B.pl.ylabel('')
        B.pl.title('')
        B.pl.tight_layout()

        B.pl.savefig(dir_name1+'/redXsec_thnq%i_deg.pdf'%(ithnq))

        #B.pl.show()

        #----------------------------PLOT RELATIVE ERRORS------------------------
        B.pl.figure(figsize=(10,7))

        yref = np.array((red_dataXsec_avg_masked_3to4[thnq==ithnq])) * 0.0
        B.plot_exp(pm_avg_3to4[thnq==ithnq], yref, tot_stats_err_3to4_m[thnq==ithnq]*100, marker='o', c='b', label='Statistical Error')
        B.plot_exp(pm_avg_3to4[thnq==ithnq], yref, tot_syst_err_3to4_m[thnq==ithnq]*100, marker='o', c='r', label='Systematics Error')
        B.plot_exp(pm_avg_3to4[thnq==ithnq], yref, tot_err_3to4_m[thnq==ithnq]*100, marker='o', c='k', label='Total Error')

        B.pl.title(r'Relative Errors: $\theta_{nq}=%i\pm5^{\circ}$, $Q^{2}=3.5\pm0.5$ GeV$^{2}$'%(ithnq), fontsize=20)
        B.pl.xlabel(r'$p_{\mathrm{r}} $ [GeV/c]', fontsize=18)
        B.pl.ylabel(r'Relative Errors $[\%]$', fontsize=18)
        plt.xticks(fontsize=19)
        plt.yticks(fontsize=19)

        B.pl.hlines(10., 0, 1.2, colors='r', linestyles='dashed', label=r'$\pm$ 10 $\%$')
        B.pl.hlines(-10., 0, 1.2, colors='r', linestyles='dashed')
        B.pl.hlines(25., 0, 1.2, colors='b', linestyles='dashed', label=r'$\pm$ 25 $\%$')
        B.pl.hlines(-25., 0, 1.2, colors='b', linestyles='dashed')
        
        B.pl.ylim(-60,80)
        B.pl.legend(frameon=False, fontsize=12, loc='upper left')
        B.pl.savefig(dir_name2+'/relXsec_err_Q2_3to4_thnq%i_deg.pdf'%(ithnq))
        #B.pl.show()
        B.pl.tight_layout()

        B.pl.figure(figsize=(10,7))

        yref = np.array((red_dataXsec_avg_masked_4to5[thnq==ithnq])) * 0.0
        B.plot_exp(pm_avg_4to5[thnq==ithnq], yref, tot_stats_err_4to5_m[thnq==ithnq]*100, marker='o', c='b', label='Statistical Error')
        B.plot_exp(pm_avg_4to5[thnq==ithnq], yref, tot_syst_err_4to5_m[thnq==ithnq]*100, marker='o', c='r', label='Systematics Error')
        B.plot_exp(pm_avg_4to5[thnq==ithnq], yref, tot_err_4to5_m[thnq==ithnq]*100, marker='o', c='k', label='Total Error')

        B.pl.title(r'Relative Errors: $\theta_{nq}=%i\pm5^{\circ}$, $Q^{2}=4.5\pm0.5$ GeV$^{2}$'%(ithnq), fontsize=20)
        B.pl.xlabel(r'$p_{\mathrm{r}} $ [GeV/c]', fontsize=18)
        B.pl.ylabel(r'Relative Errors $[\%]$', fontsize=18)
        plt.xticks(fontsize=19)
        plt.yticks(fontsize=19)
        
        B.pl.hlines(10., 0, 1.2, colors='r', linestyles='dashed', label=r'$\pm$ 10 $\%$')
        B.pl.hlines(-10., 0, 1.2, colors='r', linestyles='dashed')
        B.pl.hlines(25., 0, 1.2, colors='b', linestyles='dashed', label=r'$\pm$ 25 $\%$')
        B.pl.hlines(-25., 0, 1.2, colors='b', linestyles='dashed')
        
        B.pl.ylim(-60,80)
        B.pl.legend(frameon=False, fontsize=12, loc='upper left')
        B.pl.savefig(dir_name2+'/relXsec_err_Q2_4to5_thnq%i_deg.pdf'%(ithnq))
        #B.pl.show()
        B.pl.tight_layout()


        
        
        #-----------------------------WRITE DATA XSEC (And other relevant quatities) to FILE----------------------------------
        #Create output files to write relative uncertainties for (Pm) bins for each th_nq setting for kinematic bin Q2 = 3.5 +/- 0.5
        fout_name3to4 = dir_name3+'/relative_errors_thnq%i_Q2_3to4.txt' % (ithnq)
        fout_3to4 = open(fout_name3to4, 'w')
        fout_3to4.write('#theta_nq = %i +/- 5 deg :: Q2 = 3.5 +/- 0.5 GeV2 :: All (*_syst) errors are relative, dsig / sig [%%],  all(*_err) are absolute. \n'  %(ithnq))
        fout_3to4.write('#pm_bin: central bin with +/- 0.02 GeV. The pwia/fsi Xsec are from Laget calculation. red_dataXsec_avg with nan values have > 50%% stats uncertainty should be ignored.\n')
        fout_3to4.write('#Each of the quantities here have been averaged over overlapping (Pm, thnq) bins for Pm=80, 580(sets 1,2) and 750(sets 1,2,3). Energy and angle units are in [Gev] and [deg]  \n')
        fout_3to4.write('#!pm_bin[i,0]/  pm_avg[f,1]/  red_dataXsec_avg[f,2]/   red_dataXsec_avg_stat_err[f,3]/   red_dataXsec_avg_syst_err[f,4]/   red_dataXsec_avg_tot_err[f,5]/   rel_stats_err[f,6]/   rel_norm_syst[f,7]/   rel_kin_syst[f,8]/   rel_syst_tot[f,9]/   rel_tot_err[f,10]/\n')

        for i in range(len(pm[thnq==ithnq])):
            line = "{:<16.5f}{:<14.5f}{:<25.5E}{:<34.5E}{:<34.5E}{:<33.5E}{:<22.5f}{:<22.5f}{:<21.5f}{:<21.5f}{:<21.5f}\n".format(
                                                                                                                                  (pm[thnq==ithnq])[i], (pm_avg_3to4[thnq==ithnq])[i],
                                                                                                                                  (red_dataXsec_avg_masked_3to4[thnq==ithnq])[i]*MeV2fm,
                                                                                                                                  (red_dataXsec_avg_err_3to4[thnq==ithnq])[i]*MeV2fm,
                                                                                                                                  (red_dataXsec_avg_syst_err_3to4[thnq==ithnq])[i]*MeV2fm,
                                                                                                                                  (red_dataXsec_avg_tot_err_3to4[thnq==ithnq])[i]*MeV2fm,
                                                                                                                                  (tot_stats_err_3to4[thnq==ithnq])[i], (norm_syst_tot_3to4[thnq==ithnq])[i],
                                                                                                                                  (kin_syst_tot_3to4[thnq==ithnq])[i], (tot_syst_err_3to4[thnq==ithnq])[i],
                                                                                                                                  (tot_err_3to4[thnq==ithnq])[i])

            fout_3to4.write(line)
        fout_3to4.close()

      #Create output files to write relative uncertainties for (Pm) bins for each th_nq setting for kinematic bin Q2 = 4.5 +/- 0.5
        fout_name4to5 = dir_name3+'/relative_errors_thnq%i_Q2_4to5.txt' % (ithnq)
        fout_4to5 = open(fout_name4to5, 'w')
        fout_4to5.write('#theta_nq = %i +/- 5 deg :: Q2 = 4.5 +/- 0.5 GeV2 :: All (*_syst) errors are relative, dsig / sig [%%],  all(*_err) are absolute. \n'  %(ithnq))
        fout_4to5.write('#pm_bin: central bin with +/- 0.02 GeV. The pwia/fsi Xsec are from Laget calculation. red_dataXsec_avg with nan values have > 50%% stats uncertainty should be ignored.\n')
        fout_4to5.write('#Each of the quantities here have been averaged over overlapping (Pm, thnq) bins for Pm=80, 580(sets 1,2) and 750(sets 1,2,3). Energy and angle units are in [Gev] and [deg]  \n')
        fout_4to5.write('#!pm_bin[i,0]/  pm_avg[f,1]/  red_dataXsec_avg[f,2]/   red_dataXsec_avg_stat_err[f,3]/   red_dataXsec_avg_syst_err[f,4]/   red_dataXsec_avg_tot_err[f,5]/   rel_stats_err[f,6]/   rel_norm_syst[f,7]/   rel_kin_syst[f,8]/   rel_syst_tot[f,9]/   rel_tot_err[f,10]/\n')

        for i in range(len(pm[thnq==ithnq])):
            line = "{:<16.5f}{:<14.5f}{:<25.5E}{:<34.5E}{:<34.5E}{:<33.5E}{:<22.5f}{:<22.5f}{:<21.5f}{:<21.5f}{:<21.5f}\n".format(
                                                                                                                                  (pm[thnq==ithnq])[i], (pm_avg_4to5[thnq==ithnq])[i],
                                                                                                                                  (red_dataXsec_avg_masked_4to5[thnq==ithnq])[i]*MeV2fm,
                                                                                                                                  (red_dataXsec_avg_err_4to5[thnq==ithnq])[i]*MeV2fm,
                                                                                                                                  (red_dataXsec_avg_syst_err_4to5[thnq==ithnq])[i]*MeV2fm,
                                                                                                                                  (red_dataXsec_avg_tot_err_4to5[thnq==ithnq])[i]*MeV2fm,
                                                                                                                                  (tot_stats_err_4to5[thnq==ithnq])[i], (norm_syst_tot_4to5[thnq==ithnq])[i],
                                                                                                                                  (kin_syst_tot_4to5[thnq==ithnq])[i], (tot_syst_err_4to5[thnq==ithnq])[i],
                                                                                                                                  (tot_err_4to5[thnq==ithnq])[i])

            fout_4to5.write(line)
        fout_4to5.close()

        
    
    '''
    #Read Models at all other angles
    #Read All Other Theoretical Models At all angles (V18, CD-BONN)

    pm_avg1, red_pwiaXsec_V18_35 = read_theoretical_models("V18", "PWIA", 35)   
    pm_avg2, red_fsiXsec_V18_35 = read_theoretical_models("V18", "FSI", 35)      
    pm_avg3, red_pwiaXsec_CD_35 = read_theoretical_models("CD-Bonn", "PWIA", 35)                                 
    pm_avg4, red_fsiXsec_CD_35 = read_theoretical_models("CD-Bonn", "FSI", 35)   

    pm_avg5, red_pwiaXsec_V18_45 = read_theoretical_models("V18", "PWIA", 45)   
    pm_avg6, red_fsiXsec_V18_45 = read_theoretical_models("V18", "FSI", 45)      
    pm_avg7, red_pwiaXsec_CD_45 = read_theoretical_models("CD-Bonn", "PWIA", 45)                                 
    pm_avg8, red_fsiXsec_CD_45 = read_theoretical_models("CD-Bonn", "FSI", 45)

    pm_avg9, red_pwiaXsec_V18_75 = read_theoretical_models("V18", "PWIA", 75)   
    pm_avg10, red_fsiXsec_V18_75 = read_theoretical_models("V18", "FSI", 75)      
    pm_avg11, red_pwiaXsec_CD_75 = read_theoretical_models("CD-Bonn", "PWIA", 75)                                 
    pm_avg12, red_fsiXsec_CD_75 = read_theoretical_models("CD-Bonn", "FSI", 75)

    #Convert to sig_reduced from MeV^-3 to fm^3
    red_pwiaXsec_avg_35 = np.array(red_pwiaXsec_avg[thnq==35])*MeV2fm    #Laget PWIA
    red_fsiXsec_avg_35 = np.array(red_fsiXsec_avg[thnq==35])*MeV2fm      #Laget FSI
    red_pwiaXsec_V18_35 = np.array(red_pwiaXsec_V18_35)*MeV2fm
    red_fsiXsec_V18_35 = np.array(red_fsiXsec_V18_35)*MeV2fm
    red_pwiaXsec_CD_35 = np.array(red_pwiaXsec_CD_35)*MeV2fm
    red_fsiXsec_CD_35 = np.array(red_fsiXsec_CD_35)*MeV2fm

    
    red_pwiaXsec_avg_45 = red_pwiaXsec_avg[thnq==45]*MeV2fm    #Laget PWIA
    red_fsiXsec_avg_45 = red_fsiXsec_avg[thnq==45]*MeV2fm      #Laget FSI
    red_pwiaXsec_V18_45 = red_pwiaXsec_V18_45*MeV2fm
    red_fsiXsec_V18_45 = red_fsiXsec_V18_45*MeV2fm
    red_pwiaXsec_CD_45 = red_pwiaXsec_CD_45*MeV2fm
    red_fsiXsec_CD_45 = red_fsiXsec_CD_45*MeV2fm

    red_pwiaXsec_avg_75 = red_pwiaXsec_avg[thnq==75]*MeV2fm    #Laget PWIA
    red_fsiXsec_avg_75 = red_fsiXsec_avg[thnq==75]*MeV2fm      #Laget FSI
    red_pwiaXsec_V18_75 = red_pwiaXsec_V18_75*MeV2fm
    red_fsiXsec_V18_75 = red_fsiXsec_V18_75*MeV2fm
    red_pwiaXsec_CD_75 = red_pwiaXsec_CD_75*MeV2fm
    red_fsiXsec_CD_75 = red_fsiXsec_CD_75*MeV2fm
    
    #Do Interpolation to make theory curves smooth
    #35 deg
    #f_red_pwiaXsec_avg_35 = interp1d(pm_avg[thnq==35], red_pwiaXsec_avg_35,fill_value='extrapolate', kind='cubic')   #Paris (J.M. Laget Calculation)
    #f_red_fsiXsec_avg_35 = interp1d(pm_avg[thnq==35], red_fsiXsec_avg_35,fill_value='extrapolate', kind='cubic')    
 
    f_red_pwiaXsec_avg_35 = interp1d(pm_avg[thnq==35], red_pwiaXsec_avg_35,fill_value='extrapolate', kind='cubic')   #Paris (J.M. Laget Calculation)
    f_red_fsiXsec_avg_35 = interp1d(pm_avg[thnq==35], red_fsiXsec_avg_35,fill_value='extrapolate', kind='cubic')    
    
    f_red_pwiaXsec_V18_35 = interp1d(pm_avg1, red_pwiaXsec_V18_35,fill_value='extrapolate', kind='cubic')                        #AV18 (M. Sargsian calculation)
    f_red_fsiXsec_V18_35 = interp1d(pm_avg2, red_fsiXsec_V18_35,fill_value='extrapolate', kind='cubic')
    
    f_red_pwiaXsec_CD_35 = interp1d(pm_avg3, red_pwiaXsec_CD_35,fill_value='extrapolate', kind='cubic')                          #CD-Bonn (M. Sargsian calculation)
    f_red_fsiXsec_CD_35 = interp1d(pm_avg4, red_fsiXsec_CD_35,fill_value='extrapolate', kind='cubic')
   
    #45 eg
    f_red_pwiaXsec_avg_45 = interp1d(pm_avg[thnq==45], red_pwiaXsec_avg_45,fill_value='extrapolate', kind='cubic')   #Paris (J.M. Laget Calculation)
    f_red_fsiXsec_avg_45 = interp1d(pm_avg[thnq==45], red_fsiXsec_avg_45,fill_value='extrapolate', kind='cubic')    
    
    f_red_pwiaXsec_V18_45 = interp1d(pm_avg5, red_pwiaXsec_V18_45,fill_value='extrapolate', kind='cubic')                        #AV18 (M. Sargsian calculation)
    f_red_fsiXsec_V18_45 = interp1d(pm_avg6, red_fsiXsec_V18_45,fill_value='extrapolate', kind='cubic')
    
    f_red_pwiaXsec_CD_45 = interp1d(pm_avg7, red_pwiaXsec_CD_45,fill_value='extrapolate', kind='cubic')                          #CD-Bonn (M. Sargsian calculation)
    f_red_fsiXsec_CD_45 = interp1d(pm_avg8, red_fsiXsec_CD_45,fill_value='extrapolate', kind='cubic')
    #75 eg
    f_red_pwiaXsec_avg_75 = interp1d(pm_avg[thnq==75], red_pwiaXsec_avg_75,fill_value='extrapolate', kind='linear')   #Paris (J.M. Laget Calculation)
    f_red_fsiXsec_avg_75 = interp1d(pm_avg[thnq==75], red_fsiXsec_avg_75,fill_value='extrapolate', kind='linear')    
    
    f_red_pwiaXsec_V18_75 = interp1d(pm_avg9, red_pwiaXsec_V18_75,fill_value='extrapolate', kind='cubic')                        #AV18 (M. Sargsian calculation)
    f_red_fsiXsec_V18_75 = interp1d(pm_avg10, red_fsiXsec_V18_75,fill_value='extrapolate', kind='cubic')
    
    f_red_pwiaXsec_CD_75 = interp1d(pm_avg11, red_pwiaXsec_CD_75,fill_value='extrapolate', kind='cubic')                          #CD-Bonn (M. Sargsian calculation)
    f_red_fsiXsec_CD_75 = interp1d(pm_avg12, red_fsiXsec_CD_75,fill_value='extrapolate', kind='cubic')
    

    #------Get Averaged Missing Momentum-----
    pmiss_avg_35 = pm_avg[thnq==35]
    pmiss_avg_45 = pm_avg[thnq==45]
    pmiss_avg_75 = pm_avg[thnq==75]
            
    

    #====================================================PRL PLOT 1=============================================================
    
    #-----Create Subplots-----
    fig = B.pl.subplots(3, sharex=True, sharey=True, figsize=(13.4, 6)) 
    gs = gridspec.GridSpec(1, 3) 

    
    #THETA_NQ = 35 DEG
    ax0 = B.pl.subplot(gs[0])

    B.pl.text(0.25, 1e-7, r'$\theta_{nq}=35\pm5^{o}$', fontsize=19)
    B.pl.text(1.0, 1e0, r'(a)', fontsize=19)
    
    
    l2 = B.plot_exp(pm_ha35, red_dataXsec_ha35, red_dataXsec_err_ha35, marker='s',  markersize=5, color='#ff1000', markerfacecolor='white', capsize=0, logy=True,  label='Hall A Data')
    l1 = B.plot_exp(pm_avg[thnq==35], red_dataXsec_avg_masked[thnq==35]*MeV2fm, red_dataXsec_avg_tot_err[thnq==35]*MeV2fm, marker='o', markersize=5, color='k', capsize=0, markerfacecolor='k', logy=True, label='This Experiment (Hall C)' )
    B.pl.hlines(3.56356E-06, 0.920, 0.960, colors='k', linestyles='-')
    #B.pl.hlines(3.62178722, 0.920, 0.960, colors='k', linestyles='-')
    #Plot theoretical curves
    print('f_red_pwiaXsec_avg_35(pm_avg[thnq==35])=',f_red_pwiaXsec_avg_35(pm_avg[thnq==35]))
    print('f_red_fsiXsec_avg_35(pm_avg[thnq==35])=',f_red_fsiXsec_avg_35(pm_avg[thnq==35]))

    l3, = B.plot_exp(pm_avg1, f_red_pwiaXsec_avg_35(pm_avg1), linestyle='--', marker='None', color='#0000ff', logy=True, label='Paris PWIA')
    l4, = B.plot_exp(pm_avg1, f_red_fsiXsec_avg_35(pm_avg1), linestyle='-', marker='None', color='#0000ff', logy=True, label='Paris FSI')
    
    l5, = B.plot_exp(pm_avg1, f_red_pwiaXsec_V18_35(pm_avg1), linestyle='--', marker='None', color='#009000', logy=True, label='AV18 PWIA')   
    l6, = B.plot_exp(pm_avg2, f_red_fsiXsec_V18_35(pm_avg2), linestyle='-', marker='None', color='#009000', logy=True, label='AV18 FSI') 
    
    l7, = B.plot_exp(pm_avg3, f_red_pwiaXsec_CD_35(pm_avg3), linestyle='--', marker='None', color='#ff00ff', logy=True, label='CD-Bonn PWIA')     
    l8, = B.plot_exp(pm_avg4, f_red_fsiXsec_CD_35(pm_avg4), linestyle='-', marker='None', color='#ff00ff', logy=True, label='CD-Bonn FSI') 


    #Set axis limits
    B.pl.xlim(0.0, 1.2)

    #Set Tick Marks
    ax0.tick_params(axis='both', which='both', direction='in', top=True, bottom=True, left=True, right=True, labelsize=19)

    #Set Axes Labels for subplot 0
    B.pl.xlabel('')                                                                                                                                                                                                   
    B.pl.ylabel(r'$\sigma_{\mathrm{red}}$ ($\mathrm{fm}^{3}$) ', fontsize=22,  labelpad=10)                                                                                                                                                                        
    B.pl.title('') 

    #THETA_NQ = 45 DEG
    ax1 = plt.subplot(gs[1], sharey=ax0)
    
    B.pl.text(0.25, 1e-7, r'$\theta_{nq}=45\pm5^{o}$', fontsize=19)
    B.pl.text(1.0, 1e0, r'(b)', fontsize=19)

    #Remove un-necessary Y tick labels from subplot
    B.pl.setp(ax1.get_yticklabels(), visible=False)

    B.plot_exp(pm_ha45, red_dataXsec_ha45, red_dataXsec_err_ha45, marker='s', markersize=5, capsize=0, color='#ff0000', markerfacecolor='white',  logy=True, label='Hall A Data')
    B.plot_exp(pm_avg[thnq==45], red_dataXsec_avg_masked[thnq==45]*MeV2fm, red_dataXsec_avg_tot_err[thnq==45]*MeV2fm, marker='o', markersize=5, capsize=0, color='k', markerfacecolor='k', logy=True, label='This Experiment (Hall C)' )
    B.pl.hlines(6.96669E-06, 0.980, 1.02, colors='k', linestyles='-') 

    #Plot theoretical curves
    B.plot_exp(pm_avg[thnq==45], f_red_pwiaXsec_avg_45(pm_avg[thnq==45]), linestyle='--', marker='None', color='#0000ff', logy=True, label='Paris PWIA')
    B.plot_exp(pm_avg[thnq==45], f_red_fsiXsec_avg_45(pm_avg[thnq==45]), linestyle='-', marker='None', color='#0000ff', logy=True, label='Paris FSI')
    
    B.plot_exp(pm_avg5, f_red_pwiaXsec_V18_45(pm_avg5), linestyle='--', marker='None', color='#009000', logy=True, label='AV18 PWIA')   
    B.plot_exp(pm_avg6, f_red_fsiXsec_V18_45(pm_avg6), linestyle='-', marker='None', color='#009000', logy=True, label='AV18 FSI') 
    
    B.plot_exp(pm_avg7, f_red_pwiaXsec_CD_45(pm_avg7), linestyle='--', marker='None', color='#ff00ff', logy=True, label='CD-Bonn PWIA')     
    B.plot_exp(pm_avg8, f_red_fsiXsec_CD_45(pm_avg8), linestyle='-', marker='None', color='#ff00ff', logy=True, label='CD-Bonn FSI') 

    #Set axis limits
    B.pl.xlim(0.0, 1.2)

    #Set Tick Marks
    ax1.tick_params(axis='both', which='both', direction='in', top=True, bottom=True, left=True, right=True, labelsize=19)

    #Set Axes Labels for subplot 1
    B.pl.xlabel(r'$p_{\mathrm{r}}$ (GeV/c)', fontsize=22,  labelpad=10)                                                                                                                                                                                                   
    B.pl.ylabel('')                                                                                                                                                                        
    B.pl.title('')

    #THETA_NQ = 75
    ax2 = plt.subplot(gs[2], sharey=ax0)

    B.pl.text(0.25, 1e-7, r'$\theta_{nq}=75\pm5^{o}$', fontsize=19)
    B.pl.text(1.0, 1e0, r'(c)', fontsize=19)

    #Remove un-necessary Y tick labels from subplot
    B.pl.setp(ax2.get_yticklabels(), visible=False)

    B.plot_exp(pm_ha75, red_dataXsec_ha75, red_dataXsec_err_ha75, marker='s', markersize=5, capsize=0, color='#ff0000',  markerfacecolor='white', logy=True, label='Hall A Data')
    B.plot_exp(pm_avg[thnq==75], red_dataXsec_avg_masked[thnq==75]*MeV2fm, red_dataXsec_avg_tot_err[thnq==75]*MeV2fm, marker='o', markersize=5, capsize=0, color='k', markerfacecolor='k', logy=True, label='This Experiment (Hall C)' )
    B.pl.hlines(1.41619E-03, 0.360, 0.400, colors='k', linestyles='-')
    #Plot theoretical curves
    B.plot_exp(pm_avg[thnq==75], f_red_pwiaXsec_avg_75(pm_avg[thnq==75]), linestyle='--', marker='None', color='#0000ff', logy=True, label='Paris PWIA')
    B.plot_exp(pm_avg[thnq==75], f_red_fsiXsec_avg_75(pm_avg[thnq==75]), linestyle='-', marker='None', color='#0000ff', logy=True, label='Paris FSI')
    
    B.plot_exp(pm_avg9, f_red_pwiaXsec_V18_75(pm_avg9), linestyle='--', marker='None', color='#009000', logy=True, label='AV18 PWIA')   
    B.plot_exp(pm_avg10, f_red_fsiXsec_V18_75(pm_avg10), linestyle='-', marker='None', color='#009000', logy=True, label='AV18 FSI') 
    
    B.plot_exp(pm_avg11, f_red_pwiaXsec_CD_75(pm_avg11), linestyle='--', marker='None', color='#ff00ff', logy=True, label='CD-Bonn PWIA')     
    B.plot_exp(pm_avg12, f_red_fsiXsec_CD_75(pm_avg12), linestyle='-', marker='None', color='#ff00ff', logy=True, label='CD-Bonn FSI') 

    #Set axis limits
    B.pl.xlim(0.0, 1.2)

    #Set Tick Marks
    ax2.tick_params(axis='both', which='both', direction='in', top=True, bottom=True, left=True, right=True, labelsize=19)

    #Set Axes Labels for subplot 2
    B.pl.xlabel('')                                                                                                                                                                                                   
    B.pl.ylabel('')                                                                                                                                                                        
    B.pl.title('')
  
    #Remove spacing between subplots
    plt.subplots_adjust(wspace = 0.0000000001, bottom=0.13, top=0.98, left=0.085, right=0.98) #, hspace = 0.001, wspace = 0.001)
    
    #Create labels for specified plots in given order
    line_labels=['This Experiment (Hall C)', 'Hall A Data', 'Paris PWIA', 'Paris FSI', 'AV18 PWIA', 'AV18 FSI', 'CD-Bonn PWIA', 'CD-Bonn FSI']   #define legend labels
    ax2.legend([l1, l2, l3, l4, l5, l6, l7, l8], line_labels, loc='upper right', frameon=False, fontsize=14)      #subplot to use for common legend
      

    #B.pl.show()
    B.pl.savefig('./PRL_plot1.pdf')

    
    B.pl.figure()
    #Take Ratio of red. Xsec to Misak's redXsec to determine K_sig_cc1 ratio
    pm_paris, redXsec_Paris_MS = read_misak_pdist(model="Paris")
    pm_av18, redXsec_AV18_MS = read_misak_pdist(model="AV18")
    pm_cd, redXsec_CD_MS = read_misak_pdist(model="CD-Bonn")

    R_cd_35 =  f_red_pwiaXsec_CD_35(pmiss_avg_35) / redXsec_CD_MS(pmiss_avg_35) #f_red_pwiaXsec_CD_35(pm_cd) / redXsec_CD_MS(pm_cd)
    R_cd_45 =  f_red_pwiaXsec_CD_45(pmiss_avg_45) / redXsec_CD_MS(pmiss_avg_45)
    R_cd_75 =  f_red_pwiaXsec_CD_75(pmiss_avg_75) / redXsec_CD_MS(pmiss_avg_75)
    
    R_paris_35 =  f_red_pwiaXsec_avg_35(pmiss_avg_35) / redXsec_Paris_MS(pmiss_avg_35) #f_red_pwiaXsec_CD_35(pm_cd) / redXsec_CD_MS(pm_cd)
    R_paris_45 =  f_red_pwiaXsec_avg_45(pmiss_avg_45) / redXsec_Paris_MS(pmiss_avg_45)
    R_paris_75 =  f_red_pwiaXsec_avg_75(pmiss_avg_75) / redXsec_Paris_MS(pmiss_avg_75)
    
    R_AV_35 =  f_red_pwiaXsec_V18_35(pmiss_avg_35) / redXsec_AV18_MS(pmiss_avg_35) #f_red_pwiaXsec_CD_35(pm_cd) / redXsec_CD_MS(pm_cd)
    R_AV_45 =  f_red_pwiaXsec_V18_45(pmiss_avg_45) / redXsec_AV18_MS(pmiss_avg_45)
    R_AV_75 =  f_red_pwiaXsec_V18_75(pmiss_avg_75) / redXsec_AV18_MS(pmiss_avg_75)

    R_ha =  red_dataXsec_ha45 / redXsec_CD_MS(pm_ha45)
    R_ha_err =  red_dataXsec_err_ha45 / redXsec_CD_MS(pm_ha45)

    R_hc = red_dataXsec_avg_masked[thnq==45]*MeV2fm / redXsec_CD_MS(pm_avg[thnq==45])
    R_hc_err = red_dataXsec_avg_tot_err[thnq==75]*MeV2fm / redXsec_CD_MS(pm_avg[thnq==45])

    R_ref =  redXsec_CD_MS(pm_ha45) /  redXsec_CD_MS(pm_ha45)

    #B.plot_exp(pm_ha45, R_ha, R_ha_err, marker='o', markersize=5, capsize=0, color='#ff0000',   label='Hall A Data / CD-Bonn (Misak)')
    #B.plot_exp(pm_avg[thnq==45], R_hc, R_hc_err, marker='o', markersize=5, capsize=0, color='k', label='Hall C Data / CD-Bonn (Misak)')
    #B.plot_exp(pm_ha45, R_ref,  marker='None', linestyle='--', color='#ff00ff',  label='Ref.')

    B.pl.yscale('log')
    
    
    B.plot_exp(pmiss_avg_35, R_cd_35, marker='o', linestyle='--', color='r', label=r'$CD-Bonn: 35 deg$')
    B.plot_exp(pmiss_avg_45, R_cd_45, marker='^', linestyle='--', color='r', label=r'$CD-Bonn: 45 deg$')
    
    B.plot_exp(pmiss_avg_35, R_paris_35, marker='o', linestyle='--', color='b', label=r'$Paris: 35 deg$')
    B.plot_exp(pmiss_avg_45, R_paris_45, marker='^', linestyle='--', color='b', label=r'$Paris: 45 deg$')
        
    B.plot_exp(pmiss_avg_35, R_AV_35, marker='o', linestyle='--', color='g', label=r'$AV18: 35 deg$')
    B.plot_exp(pmiss_avg_45, R_AV_45, marker='^', linestyle='--', color='g', label=r'$AV18: 45 deg$')
    
    #B.plot_exp(pm_cd, redXsec_CD_MS(pm_cd), linestyle='--', marker='None', color='#ff00ff', logy=False, label='CD-Bonn PWIA MS')     
    #B.plot_exp(pm_paris, redXsec_Paris_MS(pm_paris), linestyle='--', marker='None', color='#0000ff', logy=False, label='Paris PWIA MS')     
    #B.plot_exp(pm_av18, redXsec_AV18_MS(pm_av18), linestyle='--', marker='None', color='g', logy=False, label='AV18 PWIA MS')     
    
    #B.plot_exp(pm_cd, f_red_pwiaXsec_CD_35(pm_cd), linestyle='-', marker='None', color='#ff00ff', logy=False, label='CD-Bonn PWIA')   
    #B.plot_exp(pm_paris, f_red_pwiaXsec_avg_35(pm_paris), linestyle='-', marker='None', color='#0000ff', logy=False, label='Paris PWIA')   
    #B.plot_exp(pm_av18, f_red_pwiaXsec_V18_35(pm_av18), linestyle='-', marker='None', color='g', logy=False, label='AV18 PWIA')   

    B.pl.xlabel(r'$P_{m}$ [GeV/c]')
    B.pl.ylabel(r'R = $\sigma_{red} / \sigma^{Misak}_{red}$')
    B.pl.legend()
    #B.pl.show()
    
    

    #Write PRL PLOT 1 data to file (avg. Pm, data Xsec and errors, and theory Xsec

    fout_name = 'redXsec_HallC_thnq35_deg.txt'
    fout = open(fout_name, 'w')
    comment1='#This datafile contains redXsec from Hall C Deuteron Experiment: E12-10-003\n'
    comment2='#Units: pm_avg [GeV/c] :: redXsec [fm^3],   theta_nq = 35 +\- 5 deg \n'
    header='#!pm_avg[f,0]/   data_redXsec[f,1]/   data_redXsec_tot_err[f,2]/   stats_rel_err[f,3]/  kin_rel_err[f,4]/   norm_rel_err[f,5]   tot_rel_syst_err[f,6]/   tot_rel_err[f,7]/ \n'
    fout.write(comment1)
    fout.write(comment2)
    fout.write(header)

    for i in range(len(pm_avg[thnq==35])):
        fout.write('%.5f %.5E  %.5E  %.5f  %.5f  %.5f  %.5f  %.5f\n' % ((pm_avg[thnq==35])[i], (red_dataXsec_avg_masked[thnq==35]*MeV2fm)[i], (red_dataXsec_avg_tot_err[thnq==35]*MeV2fm)[i], (tot_stats_err[thnq==35])[i], (kin_syst_tot[thnq==35])[i], (norm_syst_tot[thnq==35])[i], (tot_syst_err[thnq==35])[i], (tot_err[thnq==35])[i]))

    fout.close()

    fout_name = 'redXsec_HallC_thnq45_deg.txt'
    fout = open(fout_name, 'w')
    comment1='#This datafile contains redXsec from Hall C Deuteron Experiment: E12-10-003\n'
    comment2='#Units: pm_avg [GeV/c] :: redXsec [fm^3],   theta_nq = 45 +\- 5 deg \n'
    header='#!pm_avg[f,0]/   data_redXsec[f,1]/   data_redXsec_tot_err[f,2]/   stats_rel_err[f,3]/  kin_rel_err[f,4]/   norm_rel_err[f,5]   tot_rel_syst_err[f,6]/   tot_rel_err[f,7]/ \n'
    fout.write(comment1)
    fout.write(comment2)
    fout.write(header)

    for i in range(len(pm_avg[thnq==45])):
        fout.write('%.5f %.5E  %.5E  %.5f  %.5f  %.5f  %.5f  %.5f\n' % ((pm_avg[thnq==45])[i], (red_dataXsec_avg_masked[thnq==45]*MeV2fm)[i], (red_dataXsec_avg_tot_err[thnq==45]*MeV2fm)[i], (tot_stats_err[thnq==45])[i], (kin_syst_tot[thnq==45])[i], (norm_syst_tot[thnq==45])[i], (tot_syst_err[thnq==45])[i], (tot_err[thnq==45])[i]))

    fout.close()


    fout_name = 'redXsec_HallC_thnq75_deg.txt'
    fout = open(fout_name, 'w')
    comment1='#This datafile contains redXsec from Hall C Deuteron Experiment: E12-10-003\n'
    comment2='#Units: pm_avg [GeV/c] :: redXsec [fm^3],   theta_nq = 75 +\- 5 deg \n'
    header='#!pm_avg[f,0]/   data_redXsec[f,1]/   data_redXsec_tot_err[f,2]/   stats_rel_err[f,3]/  kin_rel_err[f,4]/   norm_rel_err[f,5]   tot_rel_syst_err[f,6]/   tot_rel_err[f,7]/ \n'
    fout.write(comment1)
    fout.write(comment2)
    fout.write(header)

    for i in range(len(pm_avg[thnq==75])):
        fout.write('%.5f %.5E  %.5E  %.5f  %.5f  %.5f  %.5f  %.5f\n' % ((pm_avg[thnq==75])[i], (red_dataXsec_avg_masked[thnq==75]*MeV2fm)[i], (red_dataXsec_avg_tot_err[thnq==75]*MeV2fm)[i], (tot_stats_err[thnq==75])[i], (kin_syst_tot[thnq==75])[i], (norm_syst_tot[thnq==75])[i], (tot_syst_err[thnq==75])[i], (tot_err[thnq==75])[i]))

    fout.close()
    
    #====================================================PRL PLOT 2=============================================================
    
    #Define the Data (and all models) to CD-Bonn PWIA Ratio
    #theta_nq = 35 deg
    R_ref35 = f_red_pwiaXsec_CD_35(pmiss_avg_35) / f_red_pwiaXsec_CD_35(pmiss_avg_35)
    R_data35 = red_dataXsec_avg_masked[thnq==35]*MeV2fm / f_red_pwiaXsec_CD_35(pmiss_avg_35)
    R_data_err35 = red_dataXsec_avg_tot_err[thnq==35]*MeV2fm / f_red_pwiaXsec_CD_35(pmiss_avg_35)

    R_HAdata35 = red_dataXsec_ha35 / f_red_pwiaXsec_CD_35(pm_ha35)   #Hall A data (Q2=3.5) ratio to CD-Bonn PWIA (Q2=4.5)
    R_HAdata_err35 =  red_dataXsec_err_ha35 / f_red_pwiaXsec_CD_35(pm_ha35)

    R_CDBonn_fsi35 = f_red_fsiXsec_CD_35(pmiss_avg_35) / f_red_pwiaXsec_CD_35(pmiss_avg_35)
    R_Paris_pwia35 = f_red_pwiaXsec_avg_35(pmiss_avg_35) / f_red_pwiaXsec_CD_35(pmiss_avg_35)
    R_Paris_fsi35 = f_red_fsiXsec_avg_35(pmiss_avg_35) / f_red_pwiaXsec_CD_35(pmiss_avg_35)
    R_AV18_pwia35 = f_red_pwiaXsec_V18_35(pmiss_avg_35) / f_red_pwiaXsec_CD_35(pmiss_avg_35)
    R_AV18_fsi35 = f_red_fsiXsec_V18_35(pmiss_avg_35) / f_red_pwiaXsec_CD_35(pmiss_avg_35)

    #theta_nq = 45 deg    
    R_ref45 = f_red_pwiaXsec_CD_45(pmiss_avg_45) / f_red_pwiaXsec_CD_45(pmiss_avg_45)
    R_data45 = red_dataXsec_avg_masked[thnq==45]*MeV2fm / f_red_pwiaXsec_CD_45(pmiss_avg_45)
    R_data_err45 = red_dataXsec_avg_tot_err[thnq==45]*MeV2fm / f_red_pwiaXsec_CD_45(pmiss_avg_45)

    R_HAdata45 = red_dataXsec_ha45 / f_red_pwiaXsec_CD_45(pm_ha45)   #Hall A data (Q2=3.5) ratio to CD-Bonn PWIA (Q2=4.5)
    R_HAdata_err45 =  red_dataXsec_err_ha45 / f_red_pwiaXsec_CD_45(pm_ha45)

    R_CDBonn_fsi45 = f_red_fsiXsec_CD_45(pmiss_avg_45) / f_red_pwiaXsec_CD_45(pmiss_avg_45)
    R_Paris_pwia45 = f_red_pwiaXsec_avg_45(pm_avg[thnq==45]) / f_red_pwiaXsec_CD_45(pmiss_avg_45)
    R_Paris_fsi45 = f_red_fsiXsec_avg_45(pm_avg[thnq==45]) / f_red_pwiaXsec_CD_45(pmiss_avg_45)
    R_AV18_pwia45 = f_red_pwiaXsec_V18_45(pmiss_avg_45) / f_red_pwiaXsec_CD_45(pmiss_avg_45)
    R_AV18_fsi45 = f_red_fsiXsec_V18_45(pmiss_avg_45) / f_red_pwiaXsec_CD_45(pmiss_avg_45)

    #theta_nq = 75 deg    
    R_ref75 = f_red_pwiaXsec_CD_75(pmiss_avg_75) / f_red_pwiaXsec_CD_75(pmiss_avg_75)
    R_data75 = red_dataXsec_avg_masked[thnq==75]*MeV2fm / f_red_pwiaXsec_CD_75(pmiss_avg_75)
    R_data_err75 = red_dataXsec_avg_tot_err[thnq==75]*MeV2fm / f_red_pwiaXsec_CD_75(pmiss_avg_75)

    R_HAdata75 = red_dataXsec_ha75 / f_red_pwiaXsec_CD_75(pm_ha75)   #Hall A data (Q2=3.5) ratio to CD-Bonn PWIA (Q2=4.5)
    R_HAdata_err75 =  red_dataXsec_err_ha75 / f_red_pwiaXsec_CD_75(pm_ha75)

    R_CDBonn_fsi75 = f_red_fsiXsec_CD_75(pmiss_avg_75) / f_red_pwiaXsec_CD_75(pmiss_avg_75)
    R_Paris_pwia75 = f_red_pwiaXsec_avg_75(pmiss_avg_75) / f_red_pwiaXsec_CD_75(pmiss_avg_75)
    R_Paris_fsi75 = f_red_fsiXsec_avg_75(pmiss_avg_75) / f_red_pwiaXsec_CD_75(pmiss_avg_75)
    R_AV18_pwia75 = f_red_pwiaXsec_V18_75(pmiss_avg_75) / f_red_pwiaXsec_CD_75(pmiss_avg_75)
    R_AV18_fsi75 = f_red_fsiXsec_V18_75(pmiss_avg_75) / f_red_pwiaXsec_CD_75(pmiss_avg_75)


    #B.plot_exp(pm_ha35, R_HAdata35, R_HAdata_err35, marker='o', color='r', label='Hall A Data, 35 deg', capsize=0)
    #B.plot_exp(pm_ha45, R_HAdata45, R_HAdata_err45, marker='^', color='g', label='Hall A Data, 45 deg', capsize=0)
    #B.plot_exp(pm_ha75, R_HAdata75, R_HAdata_err75, marker='s', color='b', label='Hall A Data, 75 deg', capsize=0)
    #B.pl.ylabel('Ratio')
    #B.pl.xlabel('GeV/c')
    #B.pl.title('Hall A / CD-Bonn PWIA red. Xsec (Updated Hall A data)')
    #B.pl.legend()

    #-----Create Subplots-----
    fig = B.pl.subplots(3, sharex=True, sharey=True, figsize=(6.7, 10)) 
    gs = gridspec.GridSpec(3, 1) 

    #THETA_NQ = 35 DEG
    ax0 = B.pl.subplot(gs[0])

    B.pl.text(0.1, 4, r'(a)', fontsize=19)  #original txt label
    #B.pl.text(0.1, 2, r'(a)', fontsize=19)  #ZOOM plot txt label
    #B.pl.text(1.1, 8, r'(a)', fontsize=19)  #with inset txt label

    #Remove un-necessary X tick labels from subplot
    B.pl.setp(ax0.get_xticklabels(), visible=False)

    #Plot the Data (and all models) to CD-Bonn PWIA model
    B.plot_exp(pm_avg[thnq==35], R_ref35,  marker='None', linestyle='--', color='#ff00ff', label='CD-Bonn PWIA')
       
    B.plot_exp(pm_ha35, R_HAdata35, R_HAdata_err35, marker='s', color='#ff0000', markerfacecolor='white',  label='Hall A Data', capsize=0)
    B.plot_exp(pm_avg[thnq==35], R_data35, R_data_err35, marker='o', color='k', label='Hall C Data', capsize=0)
    #print('pm, Rdata_35 = ',pm_avg[thnq==35],'::',R_data35)
    B.pl.hlines(3.62178722, 0.920, 0.960, colors='k', linestyles='-') 
    B.plot_exp(pm_avg[thnq==35], R_CDBonn_fsi35, marker='None', linestyle='-', color='#ff00ff', label='CD-Bonn FSI')
    B.plot_exp(pm_avg[thnq==35], R_Paris_pwia35, marker='None', linestyle='--', color='#0000ff', label='Paris PWIA')
    B.plot_exp(pm_avg[thnq==35], R_Paris_fsi35, marker='None', linestyle='-', color='#0000ff', label='Paris FSI')
    
    B.plot_exp(pm_avg[thnq==35], R_AV18_pwia35, marker='None', linestyle='--', color='#009000', label='AV18 PWIA')
    B.plot_exp(pm_avg[thnq==35], R_AV18_fsi35, marker='None', linestyle='-', color='#009000', label='AV18 FSI')
 
    #Set axis limits
    B.pl.xlim(0.0, 1.2)
    #B.pl.ylim(0.5, 2.3) ZOOM 
    B.pl.ylim(-0.5, 5.8)  #ORIGINAL

    #log settings
    #B.pl.ylim(0.5, 2.3)
    B.pl.yscale('linear')

    #Set Tick Marks
    ax0.tick_params(axis='both', which='both', direction='in', top=True, bottom=True, left=True, right=True, labelsize=19)    

    #Set Axes Labels for subplot 0
    B.pl.xlabel('')                                                                                                                                                                                                   
    B.pl.ylabel('')                                                                                                                                                                        
    B.pl.title('') 

    
    #----CRTE INSET PLOT(35 DEG)
    axins35 = inset_axes(ax0, width='50%', height='50%', loc='upper left')
    B.pl.xlim(0,0.8)
    B.pl.ylim(0.3,2)
    B.plot_exp(pm_avg[thnq==35], R_ref35,  marker='None', linestyle='--', color='#ff00ff')
    B.plot_exp(pm_ha35, R_HAdata35, R_HAdata_err35, marker='o', color='#ff0000',   capsize=0)
    B.plot_exp(pm_avg[thnq==35], R_data35, R_data_err35, marker='o', color='k', capsize=0)

    B.plot_exp(pm_avg[thnq==35], R_CDBonn_fsi35, marker='None', linestyle='-', color='#ff00ff')
    B.plot_exp(pm_avg[thnq==35], R_Paris_pwia35, marker='None', linestyle='--', color='#0000ff')
    B.plot_exp(pm_avg[thnq==35], R_Paris_fsi35, marker='None', linestyle='-', color='#0000ff')
    
    B.plot_exp(pm_avg[thnq==35], R_AV18_pwia35, marker='None', linestyle='--', color='#009000')
    B.plot_exp(pm_avg[thnq==35], R_AV18_fsi35, marker='None', linestyle='-', color='#009000')
    #axins35.tick_params(labelleft=False, labelbottom=False)
    axins35.yaxis.tick_right()
    

    B.pl.xlabel('')
    B.pl.ylabel('')
    B.pl.title('')

    #THETA_NQ = 45 DEG
    ax1 = B.pl.subplot(gs[1])

    B.pl.text(0.1, 16, r'(b)', fontsize=19)  #original  txt
    #B.pl.text(0.1, 1.8, r'(b)', fontsize=19)   #new txt
    #B.pl.text(1.1, 16, r'(b)', fontsize=19)  #with inset txt label

    #Remove un-necessary X tick labels from subplot
    B.pl.setp(ax1.get_xticklabels(), visible=False)
    
    #Plot the Data (and all models) to CD-Bonn PWIA model
    B.plot_exp(pm_avg[thnq==45], R_ref45,  marker='None', linestyle='--', color='#ff00ff', label='CD-Bonn PWIA')
    
    B.plot_exp(pm_ha45, R_HAdata45, R_HAdata_err45, marker='s', color='#ff0000', markerfacecolor='white',  label='Hall A Data', capsize=0)
    B.plot_exp(pm_avg[thnq==45], R_data45, R_data_err45, marker='o', color='k', label='Hall C Data)', capsize=0)
    #print('pm, Rdata_45 = ',pm_avg[thnq==45],'::',R_data45)                                                                   
    B.pl.hlines(16.65557, 0.980, 1.02, colors='k', linestyles='-')  

    B.plot_exp(pm_avg[thnq==45], R_CDBonn_fsi45, marker='None', linestyle='-', color='#ff00ff', label='CD-Bonn FSI')

    B.plot_exp(pm_avg[thnq==45], R_Paris_pwia45, marker='None', linestyle='--', color='#0000ff', label='Paris PWIA')
    B.plot_exp(pm_avg[thnq==45], R_Paris_fsi45, marker='None', linestyle='-', color='#0000ff', label='Paris FSI')
    
    B.plot_exp(pm_avg[thnq==45], R_AV18_pwia45, marker='None', linestyle='--', color='#009000', label='AV18 PWIA')
    B.plot_exp(pm_avg[thnq==45], R_AV18_fsi45, marker='None', linestyle='-', color='#009000', label='AV18 FSI')
 

    #Set axis limits
    B.pl.xlim(0.0, 1.2)
    #B.pl.ylim(0.5, 2.3) ZOOM 
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

    
    #----CREATE INSET PLOT(45 DEG)
    axins45 = inset_axes(ax1, width='50%', height='50%', loc='upper left')
    B.pl.xlim(0,0.8)
    B.pl.ylim(0.3,2)
    B.plot_exp(pm_avg[thnq==45], R_ref45,  marker='None', linestyle='--', color='#ff00ff')
    B.plot_exp(pm_ha45, R_HAdata45, R_HAdata_err45, marker='o', color='#ff0000',   capsize=0)
    B.plot_exp(pm_avg[thnq==45], R_data45, R_data_err45, marker='o', color='k', capsize=0)

    B.plot_exp(pm_avg[thnq==45], R_CDBonn_fsi45, marker='None', linestyle='-', color='#ff00ff')
    B.plot_exp(pm_avg[thnq==45], R_Paris_pwia45, marker='None', linestyle='--', color='#0000ff')
    B.plot_exp(pm_avg[thnq==45], R_Paris_fsi45, marker='None', linestyle='-', color='#0000ff')
    
    B.plot_exp(pm_avg[thnq==45], R_AV18_pwia45, marker='None', linestyle='--', color='#009000')
    B.plot_exp(pm_avg[thnq==45], R_AV18_fsi45, marker='None', linestyle='-', color='#009000')
    #axins45.tick_params(labelleft=False, labelbottom=False)
    axins45.yaxis.tick_right()
    

    B.pl.xlabel('')
    B.pl.ylabel('')
    B.pl.title('')
    

    #THETA_NQ = 75 DEG
    ax2 = B.pl.subplot(gs[2])
    
    B.pl.text(0.1, 4, r'(c)', fontsize=19)   # original 
    #B.pl.text(0.1, 4, r'(c)', fontsize=19)   # new txt

    #B.pl.text(1.1, 15, r'(c)', fontsize=19)  #with inset txt label


    #Plot the Data (and all models) to CD-Bonn PWIA model
    B.plot_exp(pm_avg[thnq==35], R_ref35,  marker='None', linestyle='--', color='#ff00ff', label='CD-Bonn PWIA')
    
    B.plot_exp(pm_ha75, R_HAdata75, R_HAdata_err75, marker='s', color='#ff0000', markerfacecolor='white',  label='Hall A Data', capsize=0)
    B.plot_exp(pm_avg[thnq==75], R_data75, R_data_err75, marker='o', color='k', label='Hall C Data', capsize=0)

    #print('pm, Rdata_75 = ',pm_avg[thnq==75],'::',R_data75)                                                                     
    B.pl.hlines(3.86100677, 0.360, 0.400, colors='k', linestyles='-')   

    B.plot_exp(pm_avg[thnq==75], R_CDBonn_fsi75, marker='None', linestyle='-', color='#ff00ff', label='CD-Bonn FSI')

    B.plot_exp(pm_avg[thnq==75], R_Paris_pwia75, marker='None', linestyle='--', color='#0000ff', label='Paris PWIA')
    B.plot_exp(pm_avg[thnq==75], R_Paris_fsi75, marker='None', linestyle='-', color='#0000ff', label='Paris FSI')
    
    B.plot_exp(pm_avg[thnq==75], R_AV18_pwia75, marker='None', linestyle='--', color='#009000', label='AV18 PWIA')
    B.plot_exp(pm_avg[thnq==75], R_AV18_fsi75, marker='None', linestyle='-', color='#009000', label='AV18 FSI')
 
    #Set axis limits
    B.pl.xlim(0.0, 1.2)
    B.pl.ylim(-0.5, 5.8)

    #log settings
    #B.pl.ylim(0.1, 5.3)
    B.pl.yscale('linear')

    #Set Tick Marks
    ax2.tick_params(axis='both', which='both', direction='in', top=True, bottom=True, left=True, right=True, labelsize=19)

    #Set Axes Labels for subplot 0
    B.pl.xlabel(r'$p_{\mathrm{r}}$ (GeV/c)', fontsize=22,  labelpad=10)                                                                               
    B.pl.ylabel('')                                                                                 
    B.pl.title('') 
    
    #Remove spacing between subplots
    plt.subplots_adjust(hspace = 0.00000001, bottom=0.09, top=0.98, right=0.95, left=0.18) #, hspace = 0.001, wspace = 0.001)
    


    B.pl.legend(loc='upper right', fontsize=12, frameon=False)
    
    #----CREATE INSET PLOT(75 DEG)
    axins75 = inset_axes(ax2, width='50%', height='50%', loc='upper left')
    B.pl.xlim(0,0.6)
    B.pl.ylim(0.3,4)
    B.plot_exp(pm_avg[thnq==75], R_ref75,  marker='None', linestyle='--', color='#ff00ff')
    B.plot_exp(pm_ha75, R_HAdata75, R_HAdata_err75, marker='o', color='#ff0000',   capsize=0)
    B.plot_exp(pm_avg[thnq==75], R_data75, R_data_err75, marker='o', color='k', capsize=0)

    B.plot_exp(pm_avg[thnq==75], R_CDBonn_fsi75, marker='None', linestyle='-', color='#ff00ff')
    B.plot_exp(pm_avg[thnq==75], R_Paris_pwia75, marker='None', linestyle='--', color='#0000ff')
    B.plot_exp(pm_avg[thnq==75], R_Paris_fsi75, marker='None', linestyle='-', color='#0000ff')
    
    B.plot_exp(pm_avg[thnq==75], R_AV18_pwia75, marker='None', linestyle='--', color='#009000')
    B.plot_exp(pm_avg[thnq==75], R_AV18_fsi75, marker='None', linestyle='-', color='#009000')
    #axins75.tick_params(labelleft=False, labelbottom=False)
    axins75.yaxis.tick_right()
    

    #B.pl.xlabel('')
    #B.pl.ylabel('')
    B.pl.title('')
    #B.pl.show()
    B.pl.savefig('./PRL_plot2.pdf')
    
   
    #Write to file
    fout_name = 'redXsec_file.txt'
    fout = open(fout_name, 'w')
    comment='pm_avg in GeV/c :: redXsec in fm^3 :: R = redXsec_HallA / redXsec_CD-Bonn_PWIA\n'
    comment2='Conv. factor: 1 fm^-1 = 0.1973 GeV\n'
    header='#!pm_avg[f,0]/   redXsec_HallA[f,1]/   redXsec_CD-Bonn_PWIA[f,2]/   R[f,3]/\n'
    fout.write(comment)
    fout.write(comment2)
    fout.write(header)

    for i in range(len(pm_ha45)):
        print('ith=',i,' :: pmiss_avg=',pm_ha45[i])
        fout.write('%.5f  %.10f   %.10f   %.10f\n' % (pm_ha45[i], red_dataXsec_ha45[i], f_red_pwiaXsec_CD_45(pm_ha45)[i], R_HAdata45[i]))

    fout.close()
    '''

    #===========================================================================================================================

    
def main():
    print('Entering Main . . .')

    #plot_data_sets()

    #plot_theory_sets('pwia')
    #plot_theory_sets('fsi')
    #plot_Xsec_vs_thnq()
    #plot_final()
    plot_momentum_dist()

if __name__=="__main__":
    main()


