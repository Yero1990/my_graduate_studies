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

#Define utility functions for generating PRL plots


#Conversion factor:  1 fm = 1/ (197 MeV),   
#The reduced cross section is in MeV^-3 units
MeV2fm = 197.3**3    #convert MeV^-3 to fm^3

#___________________________________________________
def convert2NaN(arr=np.array([]), value=0):

    #Function to convert a specified value in a array to nan (not a number)

    for i in enumerate(arr):
        if arr[i[0]]==value:
            arr[i[0]] = np.nan
    return arr

#___________________________________________________
def read_halla_data(thnq=0):

    #Function to read data files containing momentum distributions from Hall A

    #Conversion factor: 1 fm = 1 / (0.1973 GeV)  or  1 fm^-1 = 0.1973 GeV
    #fminv2GeV = 1./(1./0.1973)
    ifm2GeV = 0.1973  #inverse fermi to GeV

    thnq_f = '%.1f' % (thnq)

    #Hall A data filename
    fname = './HallA_data/pm_distribution_results_q3_1_%s.data' % (thnq_f) 
    kin = dfile(fname)

    #Pm here is not the central bin value, but an average over that bin
    p_miss = np.array(kin['p_miss_av']) * ifm2GeV    
    red_dataXsec = np.array(kin['rho'])             #units fm^3
    red_dataXsec_err = np.array(kin['delta_rho'])

    return p_miss, red_dataXsec, red_dataXsec_err


#_____________________________________________________
def read_hallc_data(thnq=0, write_output=0):

    #Function to read E12-10-003 experimental red Xsec data

    fname = "../deep_data_files/Q2_4to5GeV/redXsec_combined.txt"
    #fname = "./deep_data_files/redXsec_combined.txt"

    kin = dfile(fname)
    thnq_data = np.array(kin['xb'])
    pm_bin = np.array( (kin['yb'])[thnq_data==thnq] )  #[GeV]
    pm_avg = np.array( (kin['pm_avg'])[thnq_data==thnq] ) #[GeV]

    #Get Combined Final Red Xsec for all kinematics (and convert to fm^3)
    red_dataXsec_avg     = np.array( (kin['red_dataXsec_avg'])[thnq_data==thnq] ) * MeV2fm
    red_dataXsec_avg_err = np.array( (kin['red_dataXsec_avg_err'])[thnq_data==thnq] ) * MeV2fm            #absolute statistical
    red_dataXsec_avg_syst_err = np.array( (kin['red_dataXsec_avg_syst_err'])[thnq_data==thnq] ) * MeV2fm  #absolute systematic
    red_dataXsec_avg_tot_err = np.array( (kin['red_dataXsec_avg_tot_err'])[thnq_data==thnq] ) * MeV2fm    #absolutee total redXsec error (stats + syst)

    #Require better than 50% statistics
    red_dataXsec_avg_masked = np.ma.array(red_dataXsec_avg, mask=(red_dataXsec_avg_err>0.5*red_dataXsec_avg))
    red_dataXsec_avg_masked = np.ma.filled(red_dataXsec_avg_masked.astype(float), np.nan)

    
    #Get total fractional relative errors per (Pm, thnq) bin for plotting (probably also write to table)
    kin_syst_tot = np.array( (kin['kin_syst_tot'])[thnq_data==thnq] )    #kinematic systematics
    norm_syst_tot = np.array( (kin['norm_syst_tot'])[thnq_data==thnq] )  #normalization systematics
    tot_syst_err = np.array( (kin['tot_syst_err'])[thnq_data==thnq] )    #total systematics
    tot_stats_err = np.array( (kin['tot_stats_err'])[thnq_data==thnq] )  #total statistical
    tot_err = np.array( (kin['tot_err'])[thnq_data==thnq] )              #overall error

    #if(write_output==1):
    #    return 

    return pm_bin, pm_avg, red_dataXsec_avg_masked, red_dataXsec_avg_tot_err   #Units: GeV/c, Gev/c, fm^3, fm^3
#___________________________________________________________________
def read_theoretical_models(theory="", model="", thnq=0):

    #This code read the averaged theoretical red. Xsec and returns arrays in pm and interpolated reduced Xsec
    #The interpolation values can be obtained as follows: pm_avg, f_theory = read_theoretical_models(theory="", model="", thnq=0):
    #Then, the interpolated values are:  f_theory( pm ), where pm is any missing momentum value evaluated at the function, f_theory
    #theory: V18, CD-Bonn    model: PWIA, FSI

    thnq_f = "%.2f" %(thnq)
    fname = './theoretical_models/updated_%s_models/theoryXsec_%s%s_thnq%s_combined.data' % (model, theory, model, thnq_f) #this file only reads CD-Bonn/AV18. JML will bw read separately.
    try:
        kin = dfile(fname)
    except:
        print(fname,' \n does not exist. Bypassing error')
    #print(pm_bin)
    if(model=="PWIA" and theory=="V18"):
        red_pwiaXsec_V18 = np.array(kin['red_pwiaXsec_theory']) * MeV2fm
        pm_avg = np.array(kin['pm_avg'])
        #interpolate
        f_red_pwiaXsec_V18 = interp1d(pm_avg, red_pwiaXsec_V18,fill_value='extrapolate', kind='cubic')    #AV18 (M. Sargsian calculation)
        return pm_avg, f_red_pwiaXsec_V18
    
    if(model=="PWIA" and theory=="CD-Bonn"):                                     
        red_pwiaXsec_CD_Bonn = np.array(kin['red_pwiaXsec_theory']) * MeV2fm 
        pm_avg = np.array(kin['pm_avg'])
        #interpolate
        f_red_pwiaXsec_CD = interp1d(pm_avg, red_pwiaXsec_CD_Bonn,fill_value='extrapolate', kind='cubic') 
        return pm_avg, f_red_pwiaXsec_CD

    if(model=="FSI" and theory=="V18"):                                                 
        red_fsiXsec_V18 = np.array(kin['red_fsiXsec_theory']) * MeV2fm
        pm_avg = np.array(kin['pm_avg'])
        #interpolate
        f_red_fsiXsec_V18 = interp1d(pm_avg, red_fsiXsec_V18,fill_value='extrapolate', kind='cubic')
        return pm_avg, f_red_fsiXsec_V18

    if(model=="FSI" and theory=="CD-Bonn"):                                                                    
        red_fsiXsec_CD_Bonn = np.array(kin['red_fsiXsec_theory']) * MeV2fm
        pm_avg = np.array(kin['pm_avg'])
        #interpolate
        print('CD-Bonn FSI:', pm_avg,  red_fsiXsec_CD_Bonn )
        f_red_fsiXsec_CD = interp1d(pm_avg, red_fsiXsec_CD_Bonn,fill_value='extrapolate', kind='cubic')
        return pm_avg, f_red_fsiXsec_CD

    #---Reading J.M. Laget theory from the data file---
    if(theory=="JML"):
        fname = "../deep_data_files/Q2_4to5GeV/redXsec_combined.txt"
        #fname = "./deep_data_files/redXsec_combined.txt"

        kin = dfile(fname)
        thnq_data = np.array(kin['xb'])
        pm_bin = np.array( (kin['yb'])[thnq_data==thnq] )
        pm_avg = np.array( (kin['pm_avg'])[thnq_data==thnq] )
        
        if(model=="PWIA"):
            red_pwiaXsec_JML = np.array( (kin['red_pwiaXsec_avg'])[thnq_data==thnq] ) * MeV2fm
            if(thnq==75):
                #interpolate
                f_red_pwiaXsec_JML = interp1d(pm_avg, red_pwiaXsec_JML, fill_value='extrapolate', kind='linear')
                return pm_avg, f_red_pwiaXsec_JML
            else:
                f_red_pwiaXsec_JML = interp1d(pm_avg, red_pwiaXsec_JML, fill_value='extrapolate', kind='cubic')
                return pm_avg, f_red_pwiaXsec_JML
        
        if(model=="FSI"):
            red_fsiXsec_JML = np.array( (kin['red_fsiXsec_avg'])[thnq_data==thnq] ) * MeV2fm
            if(thnq==75):
                #interpolate
                f_red_fsiXsec_JML = interp1d(pm_avg, red_fsiXsec_JML, fill_value='extrapolate', kind='linear')
                return pm_avg, f_red_fsiXsec_JML
            else:
                f_red_fsiXsec_JML = interp1d(pm_avg, red_fsiXsec_JML, fill_value='extrapolate', kind='cubic')
                return pm_avg, f_red_fsiXsec_JML
            
#______________________________________________________________________________
def read_JWVO_theory(ithnq=0, theory='', fofa='', model=''):
    
    #Read theoretical cross sections from W.V.Orden
    #The cross sections used the WJC2, AV18 and CD-Bonn models

    #fname = '../deep_data_files/JWVOrden_calculations/JWV_Orden_redXsec_%i_deg_AVERAGE.txt' % (ithnq)
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
    
    #Perform Interpolation on theoretical calculations
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
        return[pm_avg, f_WJC2_GKex05_PWBA_red]
    elif(theory=='WJC2' and fofa=='GKex05' and model=='DWBA'):
        return[pm_avg, f_WJC2_GKex05_DWBA_red]
    elif(theory=='WJC2' and fofa=='AMT' and model=='PWBA'):
        return[pm_avg, f_WJC2_AMT_PWBA_red]
    elif(theory=='WJC2' and fofa=='AMT' and model=='DWBA'):
        return[pm_avg, f_WJC2_AMT_DWBA_red]

    elif(theory=='AV18' and fofa=='GKex05' and model=='PWBA'):
        return[pm_avg, f_AV18_GKex05_PWBA_red]
    elif(theory=='AV18' and fofa=='GKex05' and model=='DWBA'):
        return[pm_avg, f_AV18_GKex05_DWBA_red]
    elif(theory=='AV18' and fofa=='AMT' and model=='PWBA'):
        return[pm_avg, f_AV18_AMT_PWBA_red]
    elif(theory=='AV18' and fofa=='AMT' and model=='DWBA'):
        return[pm_avg, f_AV18_AMT_DWBA_red]

    elif(theory=='CD' and fofa=='GKex05' and model=='PWBA'):
        return[pm_avg, f_CD_GKex05_PWBA_red]
    elif(theory=='CD' and fofa=='GKex05' and model=='DWBA'):
        return[pm_avg, f_CD_GKex05_DWBA_red]
    elif(theory=='CD' and fofa=='AMT' and model=='PWBA'):
        return[pm_avg, f_CD_AMT_PWBA_red]
    elif(theory=='CD' and fofa=='AMT' and model=='DWBA'):
        return[pm_avg, f_CD_AMT_DWBA_red] 
