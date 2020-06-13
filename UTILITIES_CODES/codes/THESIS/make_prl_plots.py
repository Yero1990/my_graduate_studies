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


#matplotlib.use('Agg')

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

def convert2NaN(arr=np.array([]), value=0):
    #method to convert a specified value in a array to nan (not a number)
    
    for i in enumerate(arr):
        if arr[i[0]]==value:
            arr[i[0]] = np.nan
    return arr


#Conversion factor:  1 fm = 1/ (197 MeV),   
#The reduced cross section is in MeV^-3 units
MeV2fm = 197.3**3    #convert MeV^-3 to fm^3

#User Input (Dir. Name to store output)
#sys_ext = sys.argv[1]   
    
#dir_name = sys_ext+"_plots"
#dir_name_misc = dir_name + "/misc_plots"        #directory to store miscellaneous plots (trigger rates, live time, etc. . .) #obtained from the report files
#print(dir_name)
#check if directory exists, else creates it.
#if not os.path.exists(dir_name):
#    os.makedirs(dir_name)
#if not os.path.exists(dir_name_misc):
#    os.makedirs(dir_name_misc)
     

#Get Reduced Xsec Data File
fname = '../deep_data_files/Q2_4to5GeV/redXsec_combined.txt'
f = B.get_file(fname)

#Get Bin Information (Same info for all files)                                                                                                  
i_b = B.get_data(f, 'i_b')    #2D bin number                                                                    
i_x = B.get_data(f, 'i_x')    #x (th_nq) bin number                                                                                    
i_y = B.get_data(f, 'i_y')    #y (pmiss) bin number                                        
thnq = B.get_data(f, 'xb')      #th_nq value at bin center                                                          
pm =  B.get_data(f, 'yb')      #pmiss value at bin center   
pm_avg = B.get_data(f, 'pm_avg')

#Get Combined Final Red Xsec for all kinematics
red_dataXsec_avg     = np.array(B.get_data(f,'red_dataXsec_avg'))
red_dataXsec_avg_err = np.array(B.get_data(f,'red_dataXsec_avg_err'))
red_dataXsec_avg_syst_err = np.array(B.get_data(f,'red_dataXsec_avg_syst_err'))
red_dataXsec_avg_tot_err = np.array(B.get_data(f,'red_dataXsec_avg_tot_err'))
red_pwiaXsec_avg     = np.array(B.get_data(f,'red_pwiaXsec_avg'))
red_fsiXsec_avg      = np.array(B.get_data(f,'red_fsiXsec_avg'))


#Get total relative errors / (Pm, thnq) bin for plotting (probably also write to table)
kin_syst_tot = np.array(dfile(fname)['kin_syst_tot'])    #kinematic systematics
norm_syst_tot = np.array(dfile(fname)['norm_syst_tot'])  #normalization systematics
tot_syst_err = np.array(dfile(fname)['tot_syst_err'])    #total systematics
tot_stats_err = np.array(dfile(fname)['tot_stats_err'])  #total statistical
tot_err = np.array(dfile(fname)['tot_err'])              #overall error

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
    red_dataXsec = np.array(kin['rho'])             #units fm^3
    red_dataXsec_err = np.array(kin['delta_rho'])

    return p_miss, red_dataXsec, red_dataXsec_err

def read_misak_pdist(model=""):

    #The code read data files containing misak's calculation of theory momentum distirbutions using Hall C kinematics

    #Conversion factor: 1 fm = 1 / (0.1973 GeV)  or  1 fm^-1 = 0.1973 GeV
    #fminv2GeV = 1./(1./0.1973)
    GeV2fm = 0.1973**3  #inverse fermi^3 to GeV^3

    if(model=="Paris"):
        fname = '../HallA_data/mom_dist_paris_mevc.data'  
    elif(model=="AV18"):
        fname = '../HallA_data/mom_dist_v18_mevc.data'  
    elif(model=="CD-Bonn"):
        fname = '../HallA_data/mom_dist_cdbonn_mevc.data'  


    kin = dfile(fname)
    p_miss = kin['pm'] / 1000.            #convert momentum from  MeV to GeV
    red_dataXsec = kin['rho'] / GeV2fm           #Convert GeV^-3 to fm^3 units 

    #iNTERPOLATE
    f_red_misak = interp1d(p_miss, red_dataXsec, fill_value='extrapolate')  

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
        red_pwiaXsec_V18 = np.array(kin['red_pwiaXsec_theory'])
        return pm_bin, red_pwiaXsec_V18 
    if(model=="PWIA" and theory=="CD-Bonn"):                                                                                                                          
        red_pwiaXsec_CD_Bonn = np.array(kin['red_pwiaXsec_theory'])                                                                          
        return pm_bin, red_pwiaXsec_CD_Bonn
    if(model=="FSI" and theory=="V18"):                                                 
        red_fsiXsec_V18 = np.array(kin['red_fsiXsec_theory'])                                               
        return pm_bin, red_fsiXsec_V18 
    if(model=="FSI" and theory=="CD-Bonn"):                                                                                            
        red_fsiXsec_CD_Bonn = np.array(kin['red_fsiXsec_theory']) 
        return pm_bin, red_fsiXsec_CD_Bonn



def read_digitized_data():
    
    #Read Digitized Data from W.V.Orden (2014) article
    #This data shows a band with lower/upper limits of calculations for PWIA and FSI
    fname_digit_pwia = './digitized_data/PWIA_digitized_final.txt'
    fname_digit_fsi = './digitized_data/FSI_digitized_final.txt'

    pwia_dgt = dfile(fname_digit_pwia)
    fsi_dgt = dfile(fname_digit_fsi)
    
    pm_pwia_dgt = np.array(pwia_dgt['Pm'])                       #digitized PWIA missing momentum [GeV/c]
    redXsec_pwia_dgt = np.array(pwia_dgt['red_theoryXsec'])  #reduced PWIA Xsec [fm^3]
    
    pm_fsi_dgt = np.array(fsi_dgt['Pm'])    #digitized FSI missing momentum [GeV/c]
    redXsec_fsi_dgt = np.array(fsi_dgt['red_theoryXsec']) #reduced FSI Xsec [fm^3]
    
    #Define the lower/upper bounds for each array
    pm_pwia_dgt_lower =  pm_pwia_dgt[0:19]               
    redXsec_pwia_dgt_lower = redXsec_pwia_dgt[0:19]
    pm_pwia_dgt_upper =  pm_pwia_dgt[20:40]              
    redXsec_pwia_dgt_upper = redXsec_pwia_dgt[20:40]
    
    pm_fsi_dgt_lower =  pm_fsi_dgt[0:22]               
    redXsec_fsi_dgt_lower = redXsec_fsi_dgt[0:22]
    pm_fsi_dgt_upper =  pm_fsi_dgt[23:46]              
    redXsec_fsi_dgt_upper = redXsec_fsi_dgt[23:46]
    
    #interpolate
    f_pwia_lower = interp1d(pm_pwia_dgt_lower,redXsec_pwia_dgt_lower , fill_value='extrapolate', kind='linear')  
    f_pwia_upper = interp1d(pm_pwia_dgt_upper,redXsec_pwia_dgt_upper , fill_value='extrapolate', kind='linear')  
    
    f_fsi_lower = interp1d(pm_fsi_dgt_lower,redXsec_fsi_dgt_lower , fill_value='extrapolate', kind='linear')  
    f_fsi_upper = interp1d(pm_fsi_dgt_upper,redXsec_fsi_dgt_upper , fill_value='extrapolate', kind='linear')  

    return[pm_pwia_dgt_lower, pm_fsi_dgt_lower, f_pwia_lower(pm_pwia_dgt_lower), f_pwia_upper(pm_pwia_dgt_lower), f_fsi_lower(pm_fsi_dgt_lower), f_fsi_upper(pm_fsi_dgt_lower)]

    

def read_JWVO_data(ithnq=0, theory='', fofa='', model='', pmavg=np.array([])):
    
    #Read theoretical cross sections from W.V.Orden
    #The cross sections used the WJC2, AV18 and CD-Bonn models

    fname = '../deep_data_files/JWVOrden_calculations/JWV_Orden_redXsec_%i_deg_AVERAGE.txt' % (ithnq)

    f = dfile(fname)
    
    pm_avg = np.array(f['pm_avg'])                       #average missing momentum (GeV/c)
    WJC2_GKex05_PWBA_red = np.array(f['WJC2_GKex05_PWBA_red'])
    WJC2_GKex05_DWBA_red = np.array(f['WJC2_GKex05_DWBA_red'])
    WJC2_AMT_PWBA_red = np.array(f['WJC2_AMT_PWBA_red'])
    WJC2_AMT_DWBA_red = np.array(f['WJC2_AMT_DWBA_red'])

    AV18_GKex05_PWBA_red = np.array(f['AV18_GKex05_PWBA_red'])
    AV18_GKex05_DWBA_red = np.array(f['AV18_GKex05_DWBA_red'])
    AV18_AMT_PWBA_red = np.array(f['AV18_AMT_PWBA_red'])
    AV18_AMT_DWBA_red = np.array(f['AV18_AMT_DWBA_red'])

    CD_GKex05_PWBA_red = np.array(f['CD_GKex05_PWBA_red'])
    CD_GKex05_DWBA_red = np.array(f['CD_GKex05_DWBA_red'])
    CD_AMT_PWBA_red = np.array(f['CD_AMT_PWBA_red'])
    CD_AMT_DWBA_red = np.array(f['CD_AMT_DWBA_red'])
    
    #interpolate
    f_WJC2_GKex05_PWBA_red = interp1d(pm_avg, WJC2_GKex05_PWBA_red, fill_value='extrapolate', kind='linear')  
    f_WJC2_GKex05_DWBA_red = interp1d(pm_avg, WJC2_GKex05_DWBA_red, fill_value='extrapolate', kind='linear')  
    f_WJC2_AMT_PWBA_red = interp1d(pm_avg, WJC2_AMT_PWBA_red, fill_value='extrapolate', kind='linear')  
    f_WJC2_AMT_DWBA_red = interp1d(pm_avg, WJC2_AMT_DWBA_red, fill_value='extrapolate', kind='linear') 

    f_AV18_GKex05_PWBA_red = interp1d(pm_avg, AV18_GKex05_PWBA_red, fill_value='extrapolate', kind='linear')  
    f_AV18_GKex05_DWBA_red = interp1d(pm_avg, AV18_GKex05_DWBA_red, fill_value='extrapolate', kind='linear')  
    f_AV18_AMT_PWBA_red = interp1d(pm_avg, AV18_AMT_PWBA_red, fill_value='extrapolate', kind='linear')  
    f_AV18_AMT_DWBA_red = interp1d(pm_avg, AV18_AMT_DWBA_red, fill_value='extrapolate', kind='linear') 

    f_CD_GKex05_PWBA_red = interp1d(pm_avg, CD_GKex05_PWBA_red, fill_value='extrapolate', kind='linear')  
    f_CD_GKex05_DWBA_red = interp1d(pm_avg, CD_GKex05_DWBA_red, fill_value='extrapolate', kind='linear')  
    f_CD_AMT_PWBA_red = interp1d(pm_avg, CD_AMT_PWBA_red, fill_value='extrapolate', kind='linear')  
    f_CD_AMT_DWBA_red = interp1d(pm_avg, CD_AMT_DWBA_red, fill_value='extrapolate', kind='linear')

    if(theory=='WJC2' and fofa=='GKex05' and model=='PWBA'):
        return[pmavg, f_WJC2_GKex05_PWBA_red(pmavg)]
    elif(theory=='WJC2' and fofa=='GKex05' and model=='DWBA'):
        return[pmavg, f_WJC2_GKex05_DWBA_red(pmavg)]
    elif(theory=='WJC2' and fofa=='AMT' and model=='PWBA'):
        return[pmavg, f_WJC2_AMT_PWBA_red(pmavg)]
    elif(theory=='WJC2' and fofa=='AMT' and model=='DWBA'):
        return[pmavg, f_WJC2_AMT_DWBA_red(pmavg)]

    elif(theory=='AV18' and fofa=='GKex05' and model=='PWBA'):
        return[pmavg, f_AV18_GKex05_PWBA_red(pmavg)]
    elif(theory=='AV18' and fofa=='GKex05' and model=='DWBA'):
        return[pmavg, f_AV18_GKex05_DWBA_red(pmavg)]
    elif(theory=='AV18' and fofa=='AMT' and model=='PWBA'):
        return[pmavg, f_AV18_AMT_PWBA_red(pmavg)]
    elif(theory=='AV18' and fofa=='AMT' and model=='DWBA'):
        return[pmavg, f_AV18_AMT_DWBA_red(pmavg)]

    elif(theory=='CD' and fofa=='GKex05' and model=='PWBA'):
        return[pmavg, f_CD_GKex05_PWBA_red(pmavg)]
    elif(theory=='CD' and fofa=='GKex05' and model=='DWBA'):
        return[pmavg, f_CD_GKex05_DWBA_red(pmavg)]
    elif(theory=='CD' and fofa=='AMT' and model=='PWBA'):
        return[pmavg, f_CD_AMT_PWBA_red(pmavg)]
    elif(theory=='CD' and fofa=='AMT' and model=='DWBA'):
        return[pmavg, f_CD_AMT_DWBA_red(pmavg)] 

def make_prl_plots():
    #Make PRL paper plots

    

    #Read This Experiment (Hall C) Data, and require better than 50% statistics
    red_dataXsec_avg_masked = np.ma.array(red_dataXsec_avg, mask=(red_dataXsec_avg_err>0.5*red_dataXsec_avg))
    red_dataXsec_avg_masked = np.ma.filled(red_dataXsec_avg_masked.astype(float), np.nan)

    #Read Hall A data (reduced Xsec already in fm^3 units, and Precoil in GeV)
    pm_ha35, red_dataXsec_ha35, red_dataXsec_err_ha35 = read_halla_data(35)
    pm_ha45, red_dataXsec_ha45, red_dataXsec_err_ha45 = read_halla_data(45)  
    pm_ha75, red_dataXsec_ha75, red_dataXsec_err_ha75 = read_halla_data(75)  

    #Read Models 35, 45 75 (for PRL plot)
    #Read Other Theoretical Models (V18, CD-BONN)

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
    red_pwiaXsec_avg_35 = red_pwiaXsec_avg[thnq==35]*MeV2fm    #Laget PWIA
    red_fsiXsec_avg_35 = red_fsiXsec_avg[thnq==35]*MeV2fm      #Laget FSI
    red_pwiaXsec_V18_35 = red_pwiaXsec_V18_35*MeV2fm
    red_fsiXsec_V18_35 = red_fsiXsec_V18_35*MeV2fm
    red_pwiaXsec_CD_35 = red_pwiaXsec_CD_35*MeV2fm
    red_fsiXsec_CD_35 = red_fsiXsec_CD_35*MeV2fm

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

    B.pl.text(0.25, 0.5e-6, r'$\theta_{nq}=35\pm5^{o}$', fontsize=19)
    B.pl.text(1.0, 1e0, r'(a)', fontsize=19)

    #Plot Digitized Data
    #l0_pwia_band = B.pl.fill_between(read_digitized_data()[0], read_digitized_data()[2], read_digitized_data()[3], color='gray', alpha=0.3, label='W.V.Orden PWIA calculations')
    #B.pl.fill_between(pm_fsi_dgt_lower, f_fsi_lower(pm_fsi_dgt_lower),f_fsi_upper(pm_fsi_dgt_lower), color='gray', alpha=0.2, label='W.V.Orden FSI calculations')

    #Plot J.W.V. Orden Calculations using WJC2, AV18 and CD-Bonn
    WJC2_GK_PWBA, = B.plot_exp(read_JWVO_data(35,'WJC2','GKex05','PWBA',pm_avg[thnq==35])[0],  read_JWVO_data(35,'WJC2','GKex05','PWBA',pm_avg[thnq==35])[1], linestyle='--', marker='None', color='orange', logy=True, label='WJC2 (GKex05) PWBA', zorder=2 )
    WJC2_GK_DWBA, = B.plot_exp(read_JWVO_data(35,'WJC2','GKex05','DWBA',pm_avg[thnq==35])[0],  read_JWVO_data(35,'WJC2','GKex05','DWBA',pm_avg[thnq==35])[1], linestyle='-', marker='None', color='orange', logy=True, label='WJC2 (GKex05) DWBA', zorder=2 )

    WJC2_AMT_PWBA, = B.plot_exp(read_JWVO_data(35,'WJC2','AMT','PWBA',pm_avg[thnq==35])[0],  read_JWVO_data(35,'WJC2','AMT','PWBA',pm_avg[thnq==35])[1], linestyle='--', marker='None', color='darkgoldenrod', logy=True, label='WJC2 (AMT) PWBA', zorder=2 )
    WJC2_AMT_DWBA, = B.plot_exp(read_JWVO_data(35,'WJC2','AMT','DWBA',pm_avg[thnq==35])[0],  read_JWVO_data(35,'WJC2','AMT','DWBA',pm_avg[thnq==35])[1], linestyle='-', marker='None', color='darkgoldenrod', logy=True, label='WJC2 (AMT) DWBA', zorder=2 )

    AV18_GK_PWBA, = B.plot_exp(read_JWVO_data(35,'AV18','GKex05','PWBA',pm_avg[thnq==35])[0],  read_JWVO_data(35,'AV18','GKex05','PWBA',pm_avg[thnq==35])[1], linestyle='--', marker='None', color='cyan', logy=True, label='AV18 (GKex05) PWBA', zorder=2 )
    AV18_GK_DWBA, = B.plot_exp(read_JWVO_data(35,'AV18','GKex05','DWBA',pm_avg[thnq==35])[0],  read_JWVO_data(35,'AV18','GKex05','DWBA',pm_avg[thnq==35])[1], linestyle='-', marker='None', color='cyan', logy=True, label='AV18 (GKex05) DWBA', zorder=2 )

    AV18_AMT_PWBA, = B.plot_exp(read_JWVO_data(35,'AV18','AMT','PWBA',pm_avg[thnq==35])[0],  read_JWVO_data(35,'AV18','AMT','PWBA',pm_avg[thnq==35])[1], linestyle='--', marker='None', color='darkcyan', logy=True, label='AV18 (AMT) PWBA', zorder=2 )
    AV18_AMT_DWBA, = B.plot_exp(read_JWVO_data(35,'AV18','AMT','DWBA',pm_avg[thnq==35])[0],  read_JWVO_data(35,'AV18','AMT','DWBA',pm_avg[thnq==35])[1], linestyle='-', marker='None', color='darkcyan', logy=True, label='AV18 (AMT) DWBA', zorder=2 )

    CD_GK_PWBA, = B.plot_exp(read_JWVO_data(35,'CD','GKex05','PWBA',pm_avg[thnq==35])[0],  read_JWVO_data(35,'CD','GKex05','PWBA',pm_avg[thnq==35])[1], linestyle='--', marker='None', color='purple', logy=True, label='CD-Bonn (GKex05) PWBA', zorder=2 )
    CD_GK_DWBA, = B.plot_exp(read_JWVO_data(35,'CD','GKex05','DWBA',pm_avg[thnq==35])[0],  read_JWVO_data(35,'CD','GKex05','DWBA',pm_avg[thnq==35])[1], linestyle='-', marker='None', color='purple', logy=True, label='CD-Bonn (GKex05) DWBA', zorder=2 )

    CD_AMT_PWBA, = B.plot_exp(read_JWVO_data(35,'CD','AMT','PWBA',pm_avg[thnq==35])[0],  read_JWVO_data(35,'CD','AMT','PWBA',pm_avg[thnq==35])[1], linestyle='--', marker='None', color='darkviolet', logy=True, label='CD-Bonn (AMT) PWBA', zorder=2 )
    CD_AMT_DWBA, = B.plot_exp(read_JWVO_data(35,'CD','AMT','DWBA',pm_avg[thnq==35])[0],  read_JWVO_data(35,'CD','AMT','DWBA',pm_avg[thnq==35])[1], linestyle='-', marker='None', color='darkviolet', logy=True, label='CD-Bonn (AMT) DWBA', zorder=2 )

    
    #Plot Experimental Data (Hall A or Hall C)
    l2 = B.plot_exp(pm_ha35, red_dataXsec_ha35, red_dataXsec_err_ha35, marker='s',  markersize=5, color='#ff1000', markerfacecolor='white', capsize=0, logy=True,  label='Hall A Data', zorder=3)
    l1 = B.plot_exp(pm_avg[thnq==35], red_dataXsec_avg_masked[thnq==35]*MeV2fm, red_dataXsec_avg_tot_err[thnq==35]*MeV2fm, marker='o', markersize=5, color='k', capsize=0, markerfacecolor='k', logy=True, label='This Experiment (Hall C)', zorder=4)
    #B.pl.hlines(3.56356E-06, 0.920, 0.960, colors='k', linestyles='-')
    #B.pl.hlines(3.62178722, 0.920, 0.960, colors='k', linestyles='-')
    #Plot theoretical curves
    print('f_red_pwiaXsec_avg_35(pm_avg[thnq==35])=',f_red_pwiaXsec_avg_35(pm_avg[thnq==35]))
    print('f_red_fsiXsec_avg_35(pm_avg[thnq==35])=',f_red_fsiXsec_avg_35(pm_avg[thnq==35]))

    l3, = B.plot_exp(pm_avg1, f_red_pwiaXsec_avg_35(pm_avg1), linestyle='--', marker='None', color='#0000ff', logy=True, label='Paris PWIA', zorder=2)
    l4, = B.plot_exp(pm_avg1, f_red_fsiXsec_avg_35(pm_avg1), linestyle='-', marker='None', color='#0000ff', logy=True, label='Paris FSI', zorder=2)
    
    l5, = B.plot_exp(pm_avg1, f_red_pwiaXsec_V18_35(pm_avg1), linestyle='--', marker='None', color='#009000', logy=True, label='AV18 PWIA', zorder=2)   
    l6, = B.plot_exp(pm_avg2, f_red_fsiXsec_V18_35(pm_avg2), linestyle='-', marker='None', color='#009000', logy=True, label='AV18 FSI', zorder=2) 
    
    l7, = B.plot_exp(pm_avg3, f_red_pwiaXsec_CD_35(pm_avg3), linestyle='--', marker='None', color='#ff00ff', logy=True, label='CD-Bonn PWIA', zorder=2)     
    l8, = B.plot_exp(pm_avg4, f_red_fsiXsec_CD_35(pm_avg4), linestyle='-', marker='None', color='#ff00ff', logy=True, label='CD-Bonn FSI', zorder=2) 


    #Set axis limits
    B.pl.xlim(0.0, 1.2)
    B.pl.ylim(1e-7, 10.)
    
    #Set Tick Marks
    ax0.tick_params(axis='both', which='both', direction='in', top=True, bottom=True, left=True, right=True, labelsize=19)

    #Set Axes Labels for subplot 0
    B.pl.xlabel('')                                                                                                                                                                                                   
    B.pl.ylabel(r'$\sigma_{\mathrm{red}}$ [$\mathrm{fm}^{3}$] ', fontsize=22,  labelpad=10)                                                                                                                                                                        
    B.pl.title('') 

    #THETA_NQ = 45 DEG
    ax1 = plt.subplot(gs[1], sharey=ax0)
    
    B.pl.text(0.25, 0.5e-6, r'$\theta_{nq}=45\pm5^{o}$', fontsize=19)
    B.pl.text(1.0, 1e0, r'(b)', fontsize=19)

    #Remove un-necessary Y tick labels from subplot
    B.pl.setp(ax1.get_yticklabels(), visible=False)

    #Plot Digitized Data
    #B.pl.fill_between(read_digitized_data()[0], read_digitized_data()[2], read_digitized_data()[3], color='gray', alpha=0.3, label='W.V.Orden PWIA calculations')
    #B.pl.fill_between(pm_fsi_dgt_lower, f_fsi_lower(pm_fsi_dgt_lower),f_fsi_upper(pm_fsi_dgt_lower), color='gray', alpha=0.2, label='W.V.Orden FSI calculations')

    #Plot J.W.V. Orden Calculations using WJC2, AV18 and CD-Bonn
    WJC2_GK_PWBA, = B.plot_exp(read_JWVO_data(45,'WJC2','GKex05','PWBA',pm_avg[thnq==45])[0],  read_JWVO_data(45,'WJC2','GKex05','PWBA',pm_avg[thnq==45])[1], linestyle='--', marker='None', color='orange', logy=True, label='WJC2 (GKex05) PWBA', zorder=2 )
    WJC2_GK_DWBA, = B.plot_exp(read_JWVO_data(45,'WJC2','GKex05','DWBA', pm_avg[thnq==45])[0],  read_JWVO_data(45,'WJC2','GKex05','DWBA',pm_avg[thnq==45])[1], linestyle='-', marker='None', color='orange', logy=True, label='WJC2 (GKex05) DWBA', zorder=2 )

    WJC2_AMT_PWBA, = B.plot_exp(read_JWVO_data(45,'WJC2','AMT','PWBA',pm_avg[thnq==45])[0],  read_JWVO_data(45,'WJC2','AMT','PWBA',pm_avg[thnq==45])[1], linestyle='--', marker='None', color='darkgoldenrod', logy=True, label='WJC2 (AMT) PWBA', zorder=2 )
    WJC2_AMT_DWBA, = B.plot_exp(read_JWVO_data(45,'WJC2','AMT','DWBA',pm_avg[thnq==45])[0],  read_JWVO_data(45,'WJC2','AMT','DWBA',pm_avg[thnq==45])[1], linestyle='-', marker='None', color='darkgoldenrod', logy=True, label='WJC2 (AMT) DWBA', zorder=2 )

    AV18_GK_PWBA, = B.plot_exp(read_JWVO_data(45,'AV18','GKex05','PWBA',pm_avg[thnq==45])[0],  read_JWVO_data(45,'AV18','GKex05','PWBA',pm_avg[thnq==45])[1], linestyle='--', marker='None', color='cyan', logy=True, label='AV18 (GKex05) PWBA', zorder=2 )
    AV18_GK_DWBA, = B.plot_exp(read_JWVO_data(45,'AV18','GKex05','DWBA',pm_avg[thnq==45])[0],  read_JWVO_data(45,'AV18','GKex05','DWBA',pm_avg[thnq==45])[1], linestyle='-', marker='None', color='cyan', logy=True, label='AV18 (GKex05) DWBA', zorder=2 )

    AV18_AMT_PWBA, = B.plot_exp(read_JWVO_data(45,'AV18','AMT','PWBA',pm_avg[thnq==45])[0],  read_JWVO_data(45,'AV18','AMT','PWBA',pm_avg[thnq==45])[1], linestyle='--', marker='None', color='darkcyan', logy=True, label='AV18 (AMT) PWBA', zorder=2 )
    AV18_AMT_DWBA, = B.plot_exp(read_JWVO_data(45,'AV18','AMT','DWBA',pm_avg[thnq==45])[0],  read_JWVO_data(45,'AV18','AMT','DWBA',pm_avg[thnq==45])[1], linestyle='-', marker='None', color='darkcyan', logy=True, label='AV18 (AMT) DWBA', zorder=2 )

    CD_GK_PWBA, = B.plot_exp(read_JWVO_data(45,'CD','GKex05','PWBA',pm_avg[thnq==45])[0],  read_JWVO_data(45,'CD','GKex05','PWBA',pm_avg[thnq==45])[1], linestyle='--', marker='None', color='purple', logy=True, label='CD-Bonn (GKex05) PWBA', zorder=2 )
    CD_GK_DWBA, = B.plot_exp(read_JWVO_data(45,'CD','GKex05','DWBA',pm_avg[thnq==45])[0],  read_JWVO_data(45,'CD','GKex05','DWBA',pm_avg[thnq==45])[1], linestyle='-', marker='None', color='purple', logy=True, label='CD-Bonn (GKex05) DWBA', zorder=2 )

    CD_AMT_PWBA, = B.plot_exp(read_JWVO_data(45,'CD','AMT','PWBA',pm_avg[thnq==45])[0],  read_JWVO_data(45,'CD','AMT','PWBA',pm_avg[thnq==45])[1], linestyle='--', marker='None', color='darkviolet', logy=True, label='CD-Bonn (AMT) PWBA', zorder=2 )
    CD_AMT_DWBA, = B.plot_exp(read_JWVO_data(45,'CD','AMT','DWBA',pm_avg[thnq==45])[0],  read_JWVO_data(45,'CD','AMT','DWBA',pm_avg[thnq==45])[1], linestyle='-', marker='None', color='darkviolet', logy=True, label='CD-Bonn (AMT) DWBA', zorder=2 )

    
    B.plot_exp(pm_ha45, red_dataXsec_ha45, red_dataXsec_err_ha45, marker='s', markersize=5, capsize=0, color='#ff0000', markerfacecolor='white',  logy=True, label='Hall A Data', zorder=3)
    B.plot_exp(pm_avg[thnq==45], red_dataXsec_avg_masked[thnq==45]*MeV2fm, red_dataXsec_avg_tot_err[thnq==45]*MeV2fm, marker='o', markersize=5, capsize=0, color='k', markerfacecolor='k', logy=True, label='This Experiment (Hall C)', zorder=4)
    #B.pl.hlines(6.96669E-06, 0.980, 1.02, colors='k', linestyles='-') 

    #Plot theoretical curves
    B.plot_exp(pm_avg[thnq==45], f_red_pwiaXsec_avg_45(pm_avg[thnq==45]), linestyle='--', marker='None', color='#0000ff', logy=True, label='Paris PWIA', zorder=2)
    B.plot_exp(pm_avg[thnq==45], f_red_fsiXsec_avg_45(pm_avg[thnq==45]), linestyle='-', marker='None', color='#0000ff', logy=True, label='Paris FSI', zorder=2)
    
    B.plot_exp(pm_avg5, f_red_pwiaXsec_V18_45(pm_avg5), linestyle='--', marker='None', color='#009000', logy=True, label='AV18 PWIA', zorder=2)   
    B.plot_exp(pm_avg6, f_red_fsiXsec_V18_45(pm_avg6), linestyle='-', marker='None', color='#009000', logy=True, label='AV18 FSI', zorder=2) 
    
    B.plot_exp(pm_avg7, f_red_pwiaXsec_CD_45(pm_avg7), linestyle='--', marker='None', color='#ff00ff', logy=True, label='CD-Bonn PWIA', zorder=2)     
    B.plot_exp(pm_avg8, f_red_fsiXsec_CD_45(pm_avg8), linestyle='-', marker='None', color='#ff00ff', logy=True, label='CD-Bonn FSI', zorder=2) 

    #Set axis limits
    B.pl.xlim(0.0, 1.2)
    B.pl.ylim(1e-7, 10.)

    #Set Tick Marks
    ax1.tick_params(axis='both', which='both', direction='in', top=True, bottom=True, left=True, right=True, labelsize=19)

    #Set Axes Labels for subplot 1
    B.pl.xlabel(r'$p_{\mathrm{r}}$ (GeV/c)', fontsize=22,  labelpad=10)                                                                                                                                                                                                   
    B.pl.ylabel('')                                                                                                                                                                        
    B.pl.title('')

    #THETA_NQ = 75
    ax2 = plt.subplot(gs[2], sharey=ax0)

    B.pl.text(0.7, 0.5e-6, r'$\theta_{nq}=75\pm5^{o}$', fontsize=19)
    B.pl.text(2.3, 1e0, r'(c)', fontsize=19)

    #Remove un-necessary Y tick labels from subplot
    B.pl.setp(ax2.get_yticklabels(), visible=False)

    #Plot J.W.V. Orden Calculations using WJC2, AV18 and CD-Bonn
    WJC2_GK_PWBA, = B.plot_exp(read_JWVO_data(75,'WJC2','GKex05','PWBA',pm_avg[thnq==75])[0],  read_JWVO_data(75,'WJC2','GKex05','PWBA',pm_avg[thnq==75])[1], linestyle='--', marker='None', color='orange', logy=True, label='WJC2 (GKex05) PWBA', zorder=2 )
    WJC2_GK_DWBA, = B.plot_exp(read_JWVO_data(75,'WJC2','GKex05','DWBA',pm_avg[thnq==75])[0],  read_JWVO_data(75,'WJC2','GKex05','DWBA',pm_avg[thnq==75])[1], linestyle='-', marker='None', color='orange', logy=True, label='WJC2 (GKex05) DWBA', zorder=2 )

    WJC2_AMT_PWBA, = B.plot_exp(read_JWVO_data(75,'WJC2','AMT','PWBA',pm_avg[thnq==75])[0],  read_JWVO_data(75,'WJC2','AMT','PWBA',pm_avg[thnq==75])[1], linestyle='--', marker='None', color='darkgoldenrod', logy=True, label='WJC2 (AMT) PWBA', zorder=2 )
    WJC2_AMT_DWBA, = B.plot_exp(read_JWVO_data(75,'WJC2','AMT','DWBA',pm_avg[thnq==75])[0],  read_JWVO_data(75,'WJC2','AMT','DWBA',pm_avg[thnq==75])[1], linestyle='-', marker='None', color='darkgoldenrod', logy=True, label='WJC2 (AMT) DWBA', zorder=2 )

    AV18_GK_PWBA, = B.plot_exp(read_JWVO_data(75,'AV18','GKex05','PWBA',pm_avg[thnq==75])[0],  read_JWVO_data(75,'AV18','GKex05','PWBA',pm_avg[thnq==75])[1], linestyle='--', marker='None', color='cyan', logy=True, label='AV18 (GKex05) PWBA', zorder=2 )
    AV18_GK_DWBA, = B.plot_exp(read_JWVO_data(75,'AV18','GKex05','DWBA',pm_avg[thnq==75])[0],  read_JWVO_data(75,'AV18','GKex05','DWBA',pm_avg[thnq==75])[1], linestyle='-', marker='None', color='cyan', logy=True, label='AV18 (GKex05) DWBA', zorder=2 )

    AV18_AMT_PWBA, = B.plot_exp(read_JWVO_data(75,'AV18','AMT','PWBA',pm_avg[thnq==75])[0],  read_JWVO_data(75,'AV18','AMT','PWBA',pm_avg[thnq==75])[1], linestyle='--', marker='None', color='darkcyan', logy=True, label='AV18 (AMT) PWBA', zorder=2 )
    AV18_AMT_DWBA, = B.plot_exp(read_JWVO_data(75,'AV18','AMT','DWBA',pm_avg[thnq==75])[0],  read_JWVO_data(75,'AV18','AMT','DWBA',pm_avg[thnq==75])[1], linestyle='-', marker='None', color='darkcyan', logy=True, label='AV18 (AMT) DWBA', zorder=2 )

    CD_GK_PWBA, = B.plot_exp(read_JWVO_data(75,'CD','GKex05','PWBA',pm_avg[thnq==75])[0],  read_JWVO_data(75,'CD','GKex05','PWBA',pm_avg[thnq==75])[1], linestyle='--', marker='None', color='purple', logy=True, label='CD-Bonn (GKex05) PWBA', zorder=2 )
    CD_GK_DWBA, = B.plot_exp(read_JWVO_data(75,'CD','GKex05','DWBA',pm_avg[thnq==75])[0],  read_JWVO_data(75,'CD','GKex05','DWBA',pm_avg[thnq==75])[1], linestyle='-', marker='None', color='purple', logy=True, label='CD-Bonn (GKex05) DWBA', zorder=2 )

    CD_AMT_PWBA, = B.plot_exp(read_JWVO_data(75,'CD','AMT','PWBA',pm_avg[thnq==75])[0],  read_JWVO_data(75,'CD','AMT','PWBA',pm_avg[thnq==75])[1], linestyle='--', marker='None', color='darkviolet', logy=True, label='CD-Bonn (AMT) PWBA', zorder=2 )
    CD_AMT_DWBA, = B.plot_exp(read_JWVO_data(75,'CD','AMT','DWBA',pm_avg[thnq==75])[0],  read_JWVO_data(75,'CD','AMT','DWBA',pm_avg[thnq==75])[1], linestyle='-', marker='None', color='darkviolet', logy=True, label='CD-Bonn (AMT) DWBA', zorder=2 )

    
    B.plot_exp(pm_ha75, red_dataXsec_ha75, red_dataXsec_err_ha75, marker='s', markersize=5, capsize=0, color='#ff0000',  markerfacecolor='white', logy=True, label='Hall A Data',zorder=3)
    B.plot_exp(pm_avg[thnq==75], red_dataXsec_avg_masked[thnq==75]*MeV2fm, red_dataXsec_avg_tot_err[thnq==75]*MeV2fm, marker='o', markersize=5, capsize=0, color='k', markerfacecolor='k', logy=True, label='This Experiment (Hall C)', zorder=4 )
    #B.pl.hlines(1.41619E-03, 0.360, 0.400, colors='k', linestyles='-')
    #Plot theoretical curves
    B.plot_exp(pm_avg[thnq==75], f_red_pwiaXsec_avg_75(pm_avg[thnq==75]), linestyle='--', marker='None', color='#0000ff', logy=True, label='Paris PWIA', zorder=2)
    B.plot_exp(pm_avg[thnq==75], f_red_fsiXsec_avg_75(pm_avg[thnq==75]), linestyle='-', marker='None', color='#0000ff', logy=True, label='Paris FSI', zorder=2)
    
    B.plot_exp(pm_avg9, f_red_pwiaXsec_V18_75(pm_avg9), linestyle='--', marker='None', color='#009000', logy=True, label='AV18 PWIA', zorder=2)   
    B.plot_exp(pm_avg10, f_red_fsiXsec_V18_75(pm_avg10), linestyle='-', marker='None', color='#009000', logy=True, label='AV18 FSI', zorder=2) 
    
    B.plot_exp(pm_avg11, f_red_pwiaXsec_CD_75(pm_avg11), linestyle='--', marker='None', color='#ff00ff', logy=True, label='CD-Bonn PWIA', zorder=2)     
    B.plot_exp(pm_avg12, f_red_fsiXsec_CD_75(pm_avg12), linestyle='-', marker='None', color='#ff00ff', logy=True, label='CD-Bonn FSI', zorder=2) 

    #Set axis limits
    B.pl.xlim(0.0, 1.4)
    B.pl.ylim(1e-7, 10.)

    #Set Tick Marks
    ax2.tick_params(axis='both', which='both', direction='in', top=True, bottom=True, left=True, right=True, labelsize=19)

    #Set Axes Labels for subplot 2
    B.pl.xlabel('')                                                                                                                                                                                                   
    B.pl.ylabel('')                                                                                                                                                                        
    B.pl.title('')
  
    #Remove spacing between subplots
    plt.subplots_adjust(wspace = 0.0000000001, bottom=0.13, top=0.98, left=0.085, right=0.98) #, hspace = 0.001, wspace = 0.001)
    
    #Create labels for specified plots in given order
    line_labels=['This Experiment (Hall C)', 'Hall A Data', 'Paris PWIA', 'Paris FSI', 'AV18 PWIA', 'AV18 FSI', 'CD-Bonn PWIA', 'CD-Bonn FSI', 'WJC2 (GKex05) PWBA', 'WJC2 (GKex05) DWBA', 'WJC2 (AMT) PWBA', 'WJC2 (AMT) DWBA', 'AV18 (GKex05) PWBA', 'AV18 (GKex05) DWBA', 'AV18 (AMT) PWBA', 'AV18 (AMT) DWBA', 'CD-Bonn (GKex05) PWBA', 'CD-Bonn (GKex05) DWBA', 'CD-Bonn (AMT) PWBA', 'CD-Bonn (AMT) DWBA']   #define legend labels
    ax2.legend([l1, l2, l3, l4, l5, l6, l7, l8, WJC2_GK_PWBA, WJC2_GK_DWBA, WJC2_AMT_PWBA, WJC2_AMT_DWBA, AV18_GK_PWBA, AV18_GK_DWBA, AV18_AMT_PWBA, AV18_AMT_DWBA, CD_GK_PWBA, CD_GK_DWBA, CD_AMT_PWBA, CD_AMT_DWBA], line_labels, loc='upper right', frameon=False, fontsize=9)      #subplot to use for common legend
      

    B.pl.show()
    #B.pl.savefig('./PRL_plot1.pdf')


    #Write PRL PLOT 1 data to file (avg. Pm, data Xsec and errors, and theory Xsec
    '''
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
    '''
    
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

    #RATIO using the JWV Orden calculations
    #pmavg_35 = read_JWVO_data(ithnq=35, theory='CD', fofa='GKex05', model='PWBA',pmiss_avg_35)[0]
    #Ref_jvo_35 = read_JWVO_data(ithnq=35, theory='CD', fofa='GKex05', model='PWBA',pmiss_avg_35)[1]  #other reference

    #R_Ref_jvo_35 = Ref_jvo_35 / Ref_jvo_35

    R_WJC2_GK_pwba_35 =  read_JWVO_data(ithnq=35, theory='WJC2', fofa='GKex05', model='PWBA',pmavg=pmiss_avg_35)[1] /  f_red_pwiaXsec_CD_35(pmiss_avg_35)
    R_WJC2_GK_dwba_35 =  read_JWVO_data(ithnq=35, theory='WJC2', fofa='GKex05', model='DWBA',pmavg=pmiss_avg_35)[1] /  f_red_pwiaXsec_CD_35(pmiss_avg_35)
    R_WJC2_AMT_pwba_35 =  read_JWVO_data(ithnq=35, theory='WJC2', fofa='AMT', model='PWBA',pmavg=pmiss_avg_35)[1] /  f_red_pwiaXsec_CD_35(pmiss_avg_35)
    R_WJC2_AMT_dwba_35 =  read_JWVO_data(ithnq=35, theory='WJC2', fofa='AMT', model='DWBA',pmavg=pmiss_avg_35)[1] /  f_red_pwiaXsec_CD_35(pmiss_avg_35)

    R_AV18_GK_pwba_35 =  read_JWVO_data(ithnq=35, theory='AV18', fofa='GKex05', model='PWBA',pmavg=pmiss_avg_35)[1] /  f_red_pwiaXsec_CD_35(pmiss_avg_35)
    R_AV18_GK_dwba_35 =  read_JWVO_data(ithnq=35, theory='AV18', fofa='GKex05', model='DWBA',pmavg=pmiss_avg_35)[1] /  f_red_pwiaXsec_CD_35(pmiss_avg_35)
    R_AV18_AMT_pwba_35 =  read_JWVO_data(ithnq=35, theory='AV18', fofa='AMT', model='PWBA',pmavg=pmiss_avg_35)[1] /  f_red_pwiaXsec_CD_35(pmiss_avg_35)
    R_AV18_AMT_dwba_35 =  read_JWVO_data(ithnq=35, theory='AV18', fofa='AMT', model='DWBA',pmavg=pmiss_avg_35)[1] /  f_red_pwiaXsec_CD_35(pmiss_avg_35)

    R_CD_GK_pwba_35 =  read_JWVO_data(ithnq=35, theory='CD', fofa='GKex05', model='PWBA',pmavg=pmiss_avg_35)[1] /  f_red_pwiaXsec_CD_35(pmiss_avg_35)
    R_CD_GK_dwba_35 =  read_JWVO_data(ithnq=35, theory='CD', fofa='GKex05', model='DWBA',pmavg=pmiss_avg_35)[1] /  f_red_pwiaXsec_CD_35(pmiss_avg_35)
    R_CD_AMT_pwba_35 =  read_JWVO_data(ithnq=35, theory='CD', fofa='AMT', model='PWBA',pmavg=pmiss_avg_35)[1] /  f_red_pwiaXsec_CD_35(pmiss_avg_35)
    R_CD_AMT_dwba_35 =  read_JWVO_data(ithnq=35, theory='CD', fofa='AMT', model='DWBA',pmavg=pmiss_avg_35)[1] /  f_red_pwiaXsec_CD_35(pmiss_avg_35)
    
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

    #RATIO using the JWV Orden calculations
    #pmavg_45 = read_JWVO_data(ithnq=45, theory='CD', fofa='GKex05', model='PWBA')[0]
    #Ref_jvo_45 = read_JWVO_data(ithnq=45, theory='CD', fofa='GKex05', model='PWBA')[1]  #reference

    #R_Ref_jvo_45 = Ref_jvo_45 / Ref_jvo_45

    R_WJC2_GK_pwba_45 =  read_JWVO_data(ithnq=45, theory='WJC2', fofa='GKex05', model='PWBA',pmavg=pmiss_avg_45)[1] /  f_red_pwiaXsec_CD_45(pmiss_avg_45)
    R_WJC2_GK_dwba_45 =  read_JWVO_data(ithnq=45, theory='WJC2', fofa='GKex05', model='DWBA',pmavg=pmiss_avg_45)[1] /  f_red_pwiaXsec_CD_45(pmiss_avg_45)
    R_WJC2_AMT_pwba_45 =  read_JWVO_data(ithnq=45, theory='WJC2', fofa='AMT', model='PWBA',pmavg=pmiss_avg_45)[1] /  f_red_pwiaXsec_CD_45(pmiss_avg_45)
    R_WJC2_AMT_dwba_45 =  read_JWVO_data(ithnq=45, theory='WJC2', fofa='AMT', model='DWBA',pmavg=pmiss_avg_45)[1] / f_red_pwiaXsec_CD_45(pmiss_avg_45)

    R_AV18_GK_pwba_45 =  read_JWVO_data(ithnq=45, theory='AV18', fofa='GKex05', model='PWBA',pmavg=pmiss_avg_45)[1] /  f_red_pwiaXsec_CD_45(pmiss_avg_45)
    R_AV18_GK_dwba_45 =  read_JWVO_data(ithnq=45, theory='AV18', fofa='GKex05', model='DWBA',pmavg=pmiss_avg_45)[1] /  f_red_pwiaXsec_CD_45(pmiss_avg_45)
    R_AV18_AMT_pwba_45 =  read_JWVO_data(ithnq=45, theory='AV18', fofa='AMT', model='PWBA',pmavg=pmiss_avg_45)[1] /  f_red_pwiaXsec_CD_45(pmiss_avg_45)
    R_AV18_AMT_dwba_45 =  read_JWVO_data(ithnq=45, theory='AV18', fofa='AMT', model='DWBA',pmavg=pmiss_avg_45)[1] /  f_red_pwiaXsec_CD_45(pmiss_avg_45)

    R_CD_GK_pwba_45 =  read_JWVO_data(ithnq=45, theory='CD', fofa='GKex05', model='PWBA',pmavg=pmiss_avg_45)[1] / f_red_pwiaXsec_CD_45(pmiss_avg_45)
    R_CD_GK_dwba_45 =  read_JWVO_data(ithnq=45, theory='CD', fofa='GKex05', model='DWBA',pmavg=pmiss_avg_45)[1] / f_red_pwiaXsec_CD_45(pmiss_avg_45)
    R_CD_AMT_pwba_45 =  read_JWVO_data(ithnq=45, theory='CD', fofa='AMT', model='PWBA',pmavg=pmiss_avg_45)[1] /  f_red_pwiaXsec_CD_45(pmiss_avg_45)
    R_CD_AMT_dwba_45 =  read_JWVO_data(ithnq=45, theory='CD', fofa='AMT', model='DWBA',pmavg=pmiss_avg_45)[1] / f_red_pwiaXsec_CD_45(pmiss_avg_45)
    
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

    #RATIO using the JWV Orden calculations
    #pmavg_75 = read_JWVO_data(ithnq=75, theory='CD', fofa='GKex05', model='PWBA')[0]
    #Ref_jvo_75 = read_JWVO_data(ithnq=75, theory='CD', fofa='GKex05', model='PWBA')[1]  #reference

    #R_Ref_jvo_75 = Ref_jvo_75 / Ref_jvo_75

    R_WJC2_GK_pwba_75 =  read_JWVO_data(ithnq=75, theory='WJC2', fofa='GKex05', model='PWBA',pmavg=pmiss_avg_75)[1] /  f_red_pwiaXsec_CD_75(pmiss_avg_75)
    R_WJC2_GK_dwba_75 =  read_JWVO_data(ithnq=75, theory='WJC2', fofa='GKex05', model='DWBA',pmavg=pmiss_avg_75)[1] /  f_red_pwiaXsec_CD_75(pmiss_avg_75)
    R_WJC2_AMT_pwba_75 =  read_JWVO_data(ithnq=75, theory='WJC2', fofa='AMT', model='PWBA',pmavg=pmiss_avg_75)[1] / f_red_pwiaXsec_CD_75(pmiss_avg_75)
    R_WJC2_AMT_dwba_75 =  read_JWVO_data(ithnq=75, theory='WJC2', fofa='AMT', model='DWBA',pmavg=pmiss_avg_75)[1] /  f_red_pwiaXsec_CD_75(pmiss_avg_75)

    R_AV18_GK_pwba_75 =  read_JWVO_data(ithnq=75, theory='AV18', fofa='GKex05', model='PWBA',pmavg=pmiss_avg_75)[1] /  f_red_pwiaXsec_CD_75(pmiss_avg_75)
    R_AV18_GK_dwba_75 =  read_JWVO_data(ithnq=75, theory='AV18', fofa='GKex05', model='DWBA',pmavg=pmiss_avg_75)[1] /  f_red_pwiaXsec_CD_75(pmiss_avg_75)
    R_AV18_AMT_pwba_75 =  read_JWVO_data(ithnq=75, theory='AV18', fofa='AMT', model='PWBA',pmavg=pmiss_avg_75)[1] /  f_red_pwiaXsec_CD_75(pmiss_avg_75)
    R_AV18_AMT_dwba_75 =  read_JWVO_data(ithnq=75, theory='AV18', fofa='AMT', model='DWBA',pmavg=pmiss_avg_75)[1] /  f_red_pwiaXsec_CD_75(pmiss_avg_75)

    R_CD_GK_pwba_75 =  read_JWVO_data(ithnq=75, theory='CD', fofa='GKex05', model='PWBA',pmavg=pmiss_avg_75)[1] /  f_red_pwiaXsec_CD_75(pmiss_avg_75)
    R_CD_GK_dwba_75 =  read_JWVO_data(ithnq=75, theory='CD', fofa='GKex05', model='DWBA',pmavg=pmiss_avg_75)[1] /  f_red_pwiaXsec_CD_75(pmiss_avg_75)
    R_CD_AMT_pwba_75 =  read_JWVO_data(ithnq=75, theory='CD', fofa='AMT', model='PWBA',pmavg=pmiss_avg_75)[1] /  f_red_pwiaXsec_CD_75(pmiss_avg_75)
    R_CD_AMT_dwba_75 =  read_JWVO_data(ithnq=75, theory='CD', fofa='AMT', model='DWBA',pmavg=pmiss_avg_75)[1] /  f_red_pwiaXsec_CD_75(pmiss_avg_75)
    
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
    #B.plot_exp(pm_avg[thnq==35], R_ref35,  marker='None', linestyle='--', color='#ff00ff', label='CD-Bonn PWIA',zorder=2)
    B.pl.axhline(y=1.0, xmin = 0.0, xmax=1.2, color='#ff00ff', linestyle='--', label='CD-Bonn PWIA',zorder=2)
    
    B.plot_exp(pm_ha35, R_HAdata35, R_HAdata_err35, marker='s', color='#ff0000', markerfacecolor='white',  label='Hall A Data', capsize=0, zorder=3)
    B.plot_exp(pm_avg[thnq==35], R_data35, R_data_err35, marker='o', color='k', label='Hall C Data', capsize=0, zorder=4)
    B.plot_exp(pm_avg[thnq==35], R_CDBonn_fsi35, marker='None', linestyle='-', color='#ff00ff', label='CD-Bonn FSI', zorder=1)
    B.plot_exp(pm_avg[thnq==35], R_Paris_pwia35, marker='None', linestyle='--', color='#0000ff', label='Paris PWIA', zorder=1)
    B.plot_exp(pm_avg[thnq==35], R_Paris_fsi35, marker='None', linestyle='-', color='#0000ff', label='Paris FSI', zorder=1)
    
    B.plot_exp(pm_avg[thnq==35], R_AV18_pwia35, marker='None', linestyle='--', color='#009000', label='AV18 PWIA', zorder=1)
    B.plot_exp(pm_avg[thnq==35], R_AV18_fsi35, marker='None', linestyle='-', color='#009000', label='AV18 FSI', zorder=1)

    #Plot J.W.V.Orden Calculations RATIO
    #B.plot_exp(pm_avg[thnq==35], R_Ref_jvo_35, marker='None', linestyle='--', color='purple', label='CD-Bonn (GKex05) PWBA', zorder=2)

    B.plot_exp(pm_avg[thnq==35], R_WJC2_GK_pwba_35, marker='None', linestyle='--', color='orange', label='WJC2 (GKex05) PWBA', zorder=1)
    B.plot_exp(pm_avg[thnq==35], R_WJC2_GK_dwba_35, marker='None', linestyle='-', color='orange', label='WJC2 (GKex05) DWBA', zorder=1)
    B.plot_exp(pm_avg[thnq==35], R_WJC2_AMT_pwba_35, marker='None', linestyle='--', color='darkgoldenrod', label='WJC2 (AMT) PWBA', zorder=1)
    B.plot_exp(pm_avg[thnq==35], R_WJC2_AMT_dwba_35, marker='None', linestyle='-', color='darkgoldenrod', label='WJC2 (AMT) DWBA', zorder=1)

    B.plot_exp(pm_avg[thnq==35], R_AV18_GK_pwba_35, marker='None', linestyle='--', color='cyan', label='AV18 (GKex05) PWBA', zorder=1)
    B.plot_exp(pm_avg[thnq==35], R_AV18_GK_dwba_35, marker='None', linestyle='-', color='cyan', label='AV18 (GKex05) DWBA', zorder=1)
    B.plot_exp(pm_avg[thnq==35], R_AV18_AMT_pwba_35, marker='None', linestyle='--', color='darkcyan', label='AV18 (AMT) PWBA', zorder=1)
    B.plot_exp(pm_avg[thnq==35], R_AV18_AMT_dwba_35, marker='None', linestyle='-', color='darkcyan', label='AV18 (AMT) DWBA', zorder=1)

    #B.plot_exp(pm_avg[thnq==35], R_Ref_jvo_35, marker='None', linestyle='--', color='purple', label='CD (GKex05) PWBA', zorder=2)   #Reference
    B.plot_exp(pm_avg[thnq==35], R_CD_GK_pwba_35, marker='None', linestyle='-', color='purple', label='CD-Bonn (GKex05) PWBA', zorder=1)
    B.plot_exp(pm_avg[thnq==35], R_CD_GK_dwba_35, marker='None', linestyle='-', color='purple', label='CD-Bonn (GKex05) DWBA', zorder=1)
    B.plot_exp(pm_avg[thnq==35], R_CD_AMT_pwba_35, marker='None', linestyle='--', color='darkviolet', label='CD-Bonn (AMT) PWBA', zorder=1)
    B.plot_exp(pm_avg[thnq==35], R_CD_AMT_dwba_35, marker='None', linestyle='-', color='darkviolet', label='CD-Bonn (AMT) DWBA', zorder=1)

    
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

    
    '''
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
    '''
    
    #THETA_NQ = 45 DEG
    ax1 = B.pl.subplot(gs[1])

    B.pl.text(0.1, 16, r'(b)', fontsize=19)  #original  txt
    #B.pl.text(0.1, 1.8, r'(b)', fontsize=19)   #new txt
    #B.pl.text(1.1, 16, r'(b)', fontsize=19)  #with inset txt label

    #Remove un-necessary X tick labels from subplot
    B.pl.setp(ax1.get_xticklabels(), visible=False)
    
    #Plot the Data (and all models) to CD-Bonn PWIA model
    #B.plot_exp(pm_avg[thnq==45], R_ref45,  marker='None', linestyle='--', color='#ff00ff', label='CD-Bonn PWIA', zorder=2)
    B.pl.axhline(y=1.0, xmin=0.0, xmax=1.2, color='#ff00ff', linestyle='--', label='CD-Bonn PWIA',zorder=2)

    B.plot_exp(pm_ha45, R_HAdata45, R_HAdata_err45, marker='s', color='#ff0000', markerfacecolor='white',  label='Hall A Data', capsize=0, zorder=3)
    B.plot_exp(pm_avg[thnq==45], R_data45, R_data_err45, marker='o', color='k', label='Hall C Data)', capsize=0, zorder=4)

    B.plot_exp(pm_avg[thnq==45], R_CDBonn_fsi45, marker='None', linestyle='-', color='#ff00ff', label='CD-Bonn FSI', zorder=1)

    B.plot_exp(pm_avg[thnq==45], R_Paris_pwia45, marker='None', linestyle='--', color='#0000ff', label='Paris PWIA', zorder=1)
    B.plot_exp(pm_avg[thnq==45], R_Paris_fsi45, marker='None', linestyle='-', color='#0000ff', label='Paris FSI', zorder=1)
    
    B.plot_exp(pm_avg[thnq==45], R_AV18_pwia45, marker='None', linestyle='--', color='#009000', label='AV18 PWIA', zorder=1)
    B.plot_exp(pm_avg[thnq==45], R_AV18_fsi45, marker='None', linestyle='-', color='#009000', label='AV18 FSI', zorder=1)

    #Plot J.W.V.Orden Calculations RATIO
    #B.plot_exp(pm_avg[thnq==45], R_Ref_jvo_45, marker='None', linestyle='--', color='purple', label='CD-Bonn (GKex05) PWBA', zorder=2)

    B.plot_exp(pm_avg[thnq==45], R_WJC2_GK_pwba_45, marker='None', linestyle='--', color='orange', label='WJC2 (GKex05) PWBA', zorder=1)
    B.plot_exp(pm_avg[thnq==45], R_WJC2_GK_dwba_45, marker='None', linestyle='-', color='orange', label='WJC2 (GKex05) DWBA', zorder=1)
    B.plot_exp(pm_avg[thnq==45], R_WJC2_AMT_pwba_45, marker='None', linestyle='--', color='darkgoldenrod', label='WJC2 (AMT) PWBA', zorder=1)
    B.plot_exp(pm_avg[thnq==45], R_WJC2_AMT_dwba_45, marker='None', linestyle='-', color='darkgoldenrod', label='WJC2 (AMT) DWBA', zorder=1)

    B.plot_exp(pm_avg[thnq==45], R_AV18_GK_pwba_45, marker='None', linestyle='--', color='cyan', label='AV18 (GKex05) PWBA', zorder=1)
    B.plot_exp(pm_avg[thnq==45], R_AV18_GK_dwba_45, marker='None', linestyle='-', color='cyan', label='AV18 (GKex05) DWBA', zorder=1)
    B.plot_exp(pm_avg[thnq==45], R_AV18_AMT_pwba_45, marker='None', linestyle='--', color='darkcyan', label='AV18 (AMT) PWBA', zorder=1)
    B.plot_exp(pm_avg[thnq==45], R_AV18_AMT_dwba_45, marker='None', linestyle='-', color='darkcyan', label='AV18 (AMT) DWBA', zorder=1)

    B.plot_exp(pm_avg[thnq==45], R_CD_GK_pwba_45, marker='None', linestyle='-', color='purple', label='CD-Bonn (GKex05) PWBA', zorder=1)
    B.plot_exp(pm_avg[thnq==45], R_CD_GK_dwba_45, marker='None', linestyle='-', color='purple', label='CD-Bonn (GKex05) DWBA', zorder=1)
    B.plot_exp(pm_avg[thnq==45], R_CD_AMT_pwba_45, marker='None', linestyle='--', color='darkviolet', label='CD-Bonn (AMT) PWBA', zorder=1)
    B.plot_exp(pm_avg[thnq==45], R_CD_AMT_dwba_45, marker='None', linestyle='-', color='darkviolet', label='CD-Bonn (AMT) DWBA', zorder=1)

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

    '''
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
    '''

    #THETA_NQ = 75 DEG
    ax2 = B.pl.subplot(gs[2])
    
    B.pl.text(0.1, 4, r'(c)', fontsize=19)   # original 
    #B.pl.text(0.1, 4, r'(c)', fontsize=19)   # new txt
    #B.pl.text(1.1, 15, r'(c)', fontsize=19)  #with inset txt label


    #Plot the Data (and all models) to CD-Bonn PWIA model
    #B.plot_exp(pm_avg[thnq==35], R_ref35,  marker='None', linestyle='--', color='#ff00ff', label='CD-Bonn PWIA')
    B.pl.axhline(y=1.0, xmin=0.0, xmax=0.7, color='#ff00ff', linestyle='--', label='CD-Bonn PWIA',zorder=2)

    B.plot_exp(pm_ha75, R_HAdata75, R_HAdata_err75, marker='s', color='#ff0000', markerfacecolor='white',  label='Hall A Data', capsize=0,zorder=3)
    B.plot_exp(pm_avg[thnq==75], R_data75, R_data_err75, marker='o', color='k', label='Hall C Data', capsize=0,zorder=4)

    B.plot_exp(pm_avg[thnq==75], R_CDBonn_fsi75, marker='None', linestyle='-', color='#ff00ff', label='CD-Bonn FSI',zorder=1)

    B.plot_exp(pm_avg[thnq==75], R_Paris_pwia75, marker='None', linestyle='--', color='#0000ff', label='Paris PWIA',zorder=1)
    B.plot_exp(pm_avg[thnq==75], R_Paris_fsi75, marker='None', linestyle='-', color='#0000ff', label='Paris FSI',zorder=1)
    
    B.plot_exp(pm_avg[thnq==75], R_AV18_pwia75, marker='None', linestyle='--', color='#009000', label='AV18 PWIA',zorder=1)
    B.plot_exp(pm_avg[thnq==75], R_AV18_fsi75, marker='None', linestyle='-', color='#009000', label='AV18 FSI',zorder=1)

    #Plot J.W.V.Orden Calculations RATIO
    #B.plot_exp(pm_avg[thnq==75], R_Ref_jvo_75, marker='None', linestyle='--', color='purple', label='CD-Bonn (GKex05) PWBA', zorder=2)

    B.plot_exp(pm_avg[thnq==75], R_WJC2_GK_pwba_75, marker='None', linestyle='--', color='orange', label='WJC2 (GKex05) PWBA',zorder=1)
    B.plot_exp(pm_avg[thnq==75], R_WJC2_GK_dwba_75, marker='None', linestyle='-', color='orange', label='WJC2 (GKex05) DWBA', zorder=1)
    B.plot_exp(pm_avg[thnq==75], R_WJC2_AMT_pwba_75, marker='None', linestyle='--', color='darkgoldenrod', label='WJC2 (AMT) PWBA', zorder=1)
    B.plot_exp(pm_avg[thnq==75], R_WJC2_AMT_dwba_75, marker='None', linestyle='-', color='darkgoldenrod', label='WJC2 (AMT) DWBA', zorder=1)

    B.plot_exp(pm_avg[thnq==75], R_AV18_GK_pwba_75, marker='None', linestyle='--', color='cyan', label='AV18 (GKex05) PWBA', zorder=1)
    B.plot_exp(pm_avg[thnq==75], R_AV18_GK_dwba_75, marker='None', linestyle='-', color='cyan', label='AV18 (GKex05) DWBA', zorder=1)
    B.plot_exp(pm_avg[thnq==75], R_AV18_AMT_pwba_75, marker='None', linestyle='--', color='darkcyan', label='AV18 (AMT) PWBA', zorder=1)
    B.plot_exp(pm_avg[thnq==75], R_AV18_AMT_dwba_75, marker='None', linestyle='-', color='darkcyan', label='AV18 (AMT) DWBA', zorder=1)

    #B.plot_exp(pm_avg[thnq==75], R_Ref_jvo_75, marker='None', linestyle='--', color='purple', label='CD (GKex05) PWBA', zorder=2)   #Reference
    B.plot_exp(pm_avg[thnq==75], R_CD_GK_pwba_75, marker='None', linestyle='-', color='purple', label='CD (GKex05) PWBA', zorder=1)
    B.plot_exp(pm_avg[thnq==75], R_CD_GK_dwba_75, marker='None', linestyle='-', color='purple', label='CD (GKex05) DWBA', zorder=1)
    B.plot_exp(pm_avg[thnq==75], R_CD_AMT_pwba_75, marker='None', linestyle='--', color='darkviolet', label='CD-Bonn (AMT) PWBA', zorder=1)
    B.plot_exp(pm_avg[thnq==75], R_CD_AMT_dwba_75, marker='None', linestyle='-', color='darkviolet', label='CD-Bonn (AMT) DWBA', zorder=1)
    
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
    B.pl.legend(loc='upper right', fontsize=6, frameon=False)
    

    '''
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

    B.pl.xlabel('')
    B.pl.ylabel('')
    B.pl.title('')
    '''
    
    #B.pl.show()
    B.pl.savefig('./PRL_plot2.pdf')
    
    '''
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
    make_prl_plots()

if __name__=="__main__":
    main()


