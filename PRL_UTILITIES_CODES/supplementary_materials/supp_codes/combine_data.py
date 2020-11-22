import LT.box as B
from LT.datafile import dfile
import numpy as np
import numpy.ma as ma
import sys           
import os                                                      
from sys import argv  

#This Code:
#This code assumes the systematics errors have been calculated on the Xsec and red Xsec, and reads them 
#to combine the data. When combining data, ONLY used statistical weight.  To combine the systematic errors
#from different data sets, take the normal average of the systematic errors for now (Mark Jones)
#When plotting combined reduced Xsec, plot statistical and total error bars (two error bars per data point)

# **NOTE** ONLY Reduced Cross Sections (Momentum Distributions) can be combined. Cross Sections CANNOT be combined, 
#          as they depend on spectrometer the kinematic setting, whereas in momentum distribution, that depenency has
#          been factored out by dividing by K*sig_cc1(GEp, GMp)
# ** combines (takes average) different data sets of a given kinematic setting and writes to file
# ** combines different kinematic settings

# some constants                                                                                                                                                     
dtr = np.pi/180.                                                                                                                                                                               
#GeV                                                                                                                                                                             
MP = 938.272 / 1000.                                                                                                                                                                       
MN = 939.566 / 1000.                                                                                                                                                 
MD = 1875.61 / 1000.                                                                                                                                                                           
me = 0.51099 / 1000. 

def convert2NaN(arr=np.array([]), value=0):
    #method to convert a specified value in a array to nan (not a number)

    for i in enumerate(arr):
        if arr[i[0]]==value:
            arr[i[0]] = np.nan
    return arr

MeV2fm = 197.3**3    #convert MeV^-3 to fm^3

#===========User Input (Dir. Name to store output)==================
#ipython combine_data.py Em_final40MeV
sys_ext = sys.argv[1]   
#===================================================================

#check if directory exists, else creates it.
if not os.path.exists(sys_ext):
    os.makedirs(sys_ext)

output_file='./%s/redXsec_combined.txt'%(sys_ext)

fout = open(output_file, 'w') 

#------------------------------------------------------------                   
# header information for the output file    
header = """        
#This file contains combined Reduced Cross Sections from all data sets of the same kinematic bin. pm80 + pm580 (sets 1 and 2) +  pm750 (sets 1, 2 and 3)                                                                                                                                                                                           
#\\ xb = th_nq                                                                                                    
#\\ yb = pm                                                                        
#relative errors df/f=dsig/sig: kin_syst_tot (kinematic systematics),  norm_syst_tot (constant, norm. syst.)
# current header line: all averaged kinematics are either in GeV or deg : Ei_avg (beam energy), kf_avg(final e- momentum),  the_avg(final e- angle), pf_avg(final proton momentum),
# nu_avg(energy transfer, Ei-kf), Q2_avg(4-momentum transfer, GeV2), q_avg(|qlab|, magnitude of 3-momentum transfer), cthpq_cm_avg(c.o.mass angle between proton and qlab), cphi_avg(out-of-plane angle between proton and qlab)
#! i_b[i,0]/ i_x[i,1]/ i_y[i,2]/ xb[f,3]/ yb[f,4]/  pm_avg[f,5]/  kin_syst_tot[f,6]/  norm_syst_tot[f,7]/  tot_syst_err[f,8]/  tot_stats_err[f,9]/  tot_err[f,10]/   red_dataXsec_avg[f,11]/   red_dataXsec_avg_err[f,12]/  red_dataXsec_avg_syst_err[f,13]/   red_dataXsec_avg_tot_err[f,14]/  red_pwiaXsec_avg[f,15]/ red_fsiXsec_avg[f,16]/   
"""
fout.write(header)

def get_sys_fname(pm_set, data_set):
    if(pm_set==80):
        fname = './average_kinematics/%s/pm%i_fsi_norad_avgkin_systematics.txt' %(sys_ext, pm_set)
    else:
        fname = './average_kinematics/%s/pm%i_fsi_norad_avgkin_set%i_systematics.txt' %(sys_ext, pm_set, data_set)

    return fname

def get_fname(pm_set, data_set):
    if(pm_set==80):
        fname = './bin_centering_corrections/%s/pm%i_laget_bc_corr.txt' %(sys_ext, pm_set)
    else:
        fname = './bin_centering_corrections/%s/pm%i_laget_bc_corr_set%i.txt' %(sys_ext, pm_set, data_set)

    return fname

def get_pm_avg(pm_set, data_set):
    #Code to get the average momentum from the avg. kin. file
    if(pm_set==80):
        fname = './average_kinematics/%s/pm%i_fsi_norad_avgkin.txt' %(sys_ext, pm_set)
    else:
        fname = './average_kinematics/%s/pm%i_fsi_norad_avgkin_set%i.txt' %(sys_ext, pm_set, data_set)

    kin = dfile(fname)
    pm = kin['pm']
    pm_bin = kin['yb']
    lowEdge = pm_bin - 20.  #lower edged of Pm bin
    upEdge = pm_bin + 20.   #upper edge of Pm bin
    return pm_bin, pm

def get_avg_kin(header_name, pm_set, data_set):                                                                                                                                                  
    #Code to get any header and average kinematics dara array                                                                                              
    if(pm_set==80):
        fname = './average_kinematics/%s/pm%i_fsi_norad_avgkin.txt' %(sys_ext, pm_set) 
    else:
        fname = './average_kinematics/%s/pm%i_fsi_norad_avgkin_set%i.txt' %(sys_ext, pm_set, data_set)
    kin = dfile(fname)
    data = kin[header_name]                                                                                                                                                          
    return data 

def get_sig_red(header_name, pm_set, data_set):
    #Code to get any header and data array from bin-center corrected data files

    #Get file containing red. Xsec and its stats. error
    f = dfile(get_fname(pm_set, data_set))
    
    #Get data array with desired header name
    data = f[header_name]
    return data

def get_kin_syst(header_name, pm_set, data_set):
    #Code to get any header and data array from the kin. systematic data files

    #Get file containing kin. systematics
    f = dfile(get_sys_fname(pm_set, data_set))
    
    #Get data array with desired header name (kin. syst are in %)
    data = f[header_name] / 100.   #coonvert to fractional relative error
    return data

def get_norm_syst(header_name):
    #Code to get any header and data array from the normalization. systematic data file

    #Get file containing kin. systematics
    #f = dfile('../average_kinematics/Em_final40MeV/normalization_systematics.txt')
    f = dfile('./average_kinematics/%s/normalization_systematics_summary.txt' % (sys_ext))

    #Get data array with desired header name (norm. syst are in fractional rel. error)
    data = f[header_name]   
    return data

#Target wall and spectrometer aceptance contributions
#relative error on cross sections will be added quadratically to normalization (October, 19, 2020)
tgt_wall_err = 0.029   #target wall contributes at most 2.9 % to the yield (integrated over all thnq,)
spec_acc_err = 0.014   #spectrometer acceptance contributes 1.4 % to the total cross section (study by M.K.Jones)


#Get bin information from each file
f1  = dfile(get_fname(80, 1))
i_b = f1['i_b']
i_x = f1['i_x']
i_y = f1['i_y']
xb  = f1['xb']   #th_nq bin
yb  = f1['yb']   #pm_bin

#Get Average Pmiss
pm1_b, pm80 = get_pm_avg(80, 1)
pm2_b, pm580_set1 = get_pm_avg(580, 1)
pm3_b, pm580_set2 = get_pm_avg(580, 2)
pm4_b, pm750_set1 = get_pm_avg(750, 1)
pm5_b, pm750_set2 = get_pm_avg(750, 2)
pm6_b, pm750_set3 = get_pm_avg(750, 3)

#Get Average initial beam energy [MeV]     Get avg final e- energy [MeV]              Get final e- angle [deg]                
Ei_80 = get_avg_kin('Ei', 80, 1)      ;    kf_80 = get_avg_kin('kf', 80, 1)     ;     the_80 = get_avg_kin('th_e', 80, 1)
Ei_580s1 = get_avg_kin('Ei', 580, 1) ;    kf_580s1 = get_avg_kin('kf', 580, 1) ;     the_580s1 = get_avg_kin('th_e', 580, 1)
Ei_580s2 = get_avg_kin('Ei', 580, 2) ;    kf_580s2 = get_avg_kin('kf', 580, 2) ;     the_580s2 = get_avg_kin('th_e', 580, 2)
Ei_750s1 = get_avg_kin('Ei', 750, 1) ;    kf_750s1 = get_avg_kin('kf', 750, 1) ;     the_750s1 = get_avg_kin('th_e', 750, 1)
Ei_750s2 = get_avg_kin('Ei', 750, 2) ;    kf_750s2 = get_avg_kin('kf', 750, 2) ;     the_750s2 = get_avg_kin('th_e', 750, 2)
Ei_750s3 = get_avg_kin('Ei', 750, 3) ;    kf_750s3 = get_avg_kin('kf', 750, 3) ;     the_750s3 = get_avg_kin('th_e', 750, 3)
#omega, Q2_calc, q_lab, theta_pq_calc, th_pq_cm, cos_phi
#Get final proton momentum [MeV]          #Get omega average [MeV]                   #Get averaged Q2 [MeV2]
pf_80 = get_avg_kin('pf', 80, 1)     ;    nu_80 = get_avg_kin('omega', 80, 1)     ;  Q2_80 = get_avg_kin('Q2_calc', 80, 1)         ;
pf_580s1 = get_avg_kin('pf', 580, 1) ;    nu_580s1 = get_avg_kin('omega', 580, 1) ;  Q2_580s1 = get_avg_kin('Q2_calc', 580, 1)     ;
pf_580s2 = get_avg_kin('pf', 580, 2) ;    nu_580s2 = get_avg_kin('omega', 580, 2) ;  Q2_580s2 = get_avg_kin('Q2_calc', 580, 2)     ;
pf_750s1 = get_avg_kin('pf', 750, 1) ;    nu_750s1 = get_avg_kin('omega', 750, 1) ;  Q2_750s1 = get_avg_kin('Q2_calc', 750, 1)     ;
pf_750s2 = get_avg_kin('pf', 750, 2) ;    nu_750s2 = get_avg_kin('omega', 750, 2) ;  Q2_750s2 = get_avg_kin('Q2_calc', 750, 2)     ;
pf_750s3 = get_avg_kin('pf', 750, 3) ;    nu_750s3 = get_avg_kin('omega', 750, 3) ;  Q2_750s3 = get_avg_kin('Q2_calc', 750, 3)     ;

#Get final |q_lab| MeV                    #Get th_pq_cm                                         #Get cos(phi_pq)
q_80 = get_avg_kin('q_lab', 80, 1)     ;  thpq_cm_80 = get_avg_kin('th_pq_cm', 80, 1)     ;     cphi_80 = get_avg_kin('cos_phi', 80, 1)     ;
q_580s1 = get_avg_kin('q_lab', 580, 1) ;  thpq_cm_580s1 = get_avg_kin('th_pq_cm', 580, 1) ;     cphi_580s1 = get_avg_kin('cos_phi', 580, 1)     ;
q_580s2 = get_avg_kin('q_lab', 580, 2) ;  thpq_cm_580s2 = get_avg_kin('th_pq_cm', 580, 2) ;     cphi_580s2 = get_avg_kin('cos_phi', 580, 2)     ;
q_750s1 = get_avg_kin('q_lab', 750, 1) ;  thpq_cm_750s1 = get_avg_kin('th_pq_cm', 750, 1) ;     cphi_750s1 = get_avg_kin('cos_phi', 750, 1)     ;
q_750s2 = get_avg_kin('q_lab', 750, 2) ;  thpq_cm_750s2 = get_avg_kin('th_pq_cm', 750, 2) ;     cphi_750s2 = get_avg_kin('cos_phi', 750, 2)     ;
q_750s3 = get_avg_kin('q_lab', 750, 3) ;  thpq_cm_750s3 = get_avg_kin('th_pq_cm', 750, 3) ;     cphi_750s3 = get_avg_kin('cos_phi', 750, 3)     ;


 
#Get reduced Xsec
sig_red_80 = get_sig_red('red_dataXsec', 80, 1)
sig_red_580_set1 = get_sig_red('red_dataXsec', 580, 1)
sig_red_580_set2 = get_sig_red('red_dataXsec', 580, 2)
sig_red_750_set1 = get_sig_red('red_dataXsec', 750, 1)
sig_red_750_set2 = get_sig_red('red_dataXsec', 750, 2)
sig_red_750_set3 = get_sig_red('red_dataXsec', 750, 3)

#Get reduced Xsec stats. err
sig_red_80_err = get_sig_red('red_dataXsec_err', 80, 1)
sig_red_580_set1_err = get_sig_red('red_dataXsec_err', 580, 1)
sig_red_580_set2_err = get_sig_red('red_dataXsec_err', 580, 2)
sig_red_750_set1_err = get_sig_red('red_dataXsec_err', 750, 1)
sig_red_750_set2_err = get_sig_red('red_dataXsec_err', 750, 2)
sig_red_750_set3_err = get_sig_red('red_dataXsec_err', 750, 3)

#Get theoretical red. Xsec
red_pwiaXsec_80 = get_sig_red('red_pwiaXsec', 80, 1)
red_pwiaXsec_580_set1 = get_sig_red('red_pwiaXsec', 580, 1)
red_pwiaXsec_580_set2 = get_sig_red('red_pwiaXsec', 580, 2)
red_pwiaXsec_750_set1 = get_sig_red('red_pwiaXsec', 750, 1)
red_pwiaXsec_750_set2 = get_sig_red('red_pwiaXsec', 750, 2)
red_pwiaXsec_750_set3 = get_sig_red('red_pwiaXsec', 750, 3)

red_fsiXsec_80 = get_sig_red('red_fsiXsec', 80, 1)
red_fsiXsec_580_set1 = get_sig_red('red_fsiXsec', 580, 1)
red_fsiXsec_580_set2 = get_sig_red('red_fsiXsec', 580, 2)
red_fsiXsec_750_set1 = get_sig_red('red_fsiXsec', 750, 1)
red_fsiXsec_750_set2 = get_sig_red('red_fsiXsec', 750, 2)
red_fsiXsec_750_set3 = get_sig_red('red_fsiXsec', 750, 3)

#Get Kinematic Systematics Relative error
kin_syst_80 = get_kin_syst('sig_kin_tot', 80, 1)
kin_syst_580_set1 = get_kin_syst('sig_kin_tot', 580, 1)
kin_syst_580_set2 = get_kin_syst('sig_kin_tot', 580, 2)
kin_syst_750_set1 = get_kin_syst('sig_kin_tot', 750, 1)
kin_syst_750_set2 = get_kin_syst('sig_kin_tot', 750, 2)
kin_syst_750_set3 = get_kin_syst('sig_kin_tot', 750, 3)


pm_set = get_norm_syst('pm')
data_set = get_norm_syst('set')

#Get norm. syst. errors that are the same for all data sets (these should be added as an overall norm. factor, and not in quadrature, as they affect all data sets the same)
pAbs_syst = get_norm_syst('pAbs_syst')[0]   #proton absorption systematic relative error (dsig / sig), (0.4951 %)
tLT_syst = get_norm_syst('tLT_syst')[0]   #total live time systematics (3 %)
Qtot_syst = get_norm_syst('Qtot_syst')[0]   #total charge systematics (2 %)
const_norm_syst = np.sqrt(pAbs_syst**2 + tLT_syst**2 + Qtot_syst**2 + tgt_wall_err**2 + spec_acc_err**2)   #constant norm. systematics to be added as a single value later

print('pAbs_syst = ',pAbs_syst)
print('tLT_syst = ', tLT_syst)
print('Qtot_syst = ', Qtot_syst)
print('const_norm_syst = ', const_norm_syst)

#Get norm. syst. errors that change per data set 
htrk_eff_syst = get_norm_syst('htrk_eff_syst')
etrk_eff_syst = get_norm_syst('etrk_eff_syst')
tgtBoil_syst = get_norm_syst('tgtBoil_syst')
print('htrk_eff_syst = ',htrk_eff_syst)
print('etrk_eff_syst = ',etrk_eff_syst)
print('tgtBoil_syst = ',tgtBoil_syst)

#norm_syst_tot = get_norm_syst('total_norm_syst')  #new syst. array (systematcis for each separate set, added in quadrature),  norm_syst_tot[i] -> i=0: 80 MeV setting, i=1: 580 (set1) setting, . . . 

norm_syst_80 = np.sqrt(htrk_eff_syst[0]**2 + etrk_eff_syst[0]**2 + tgtBoil_syst[0]**2)
norm_syst_580_set1 = np.sqrt(htrk_eff_syst[1]**2 + etrk_eff_syst[1]**2 + tgtBoil_syst[1]**2)
norm_syst_580_set2 = np.sqrt(htrk_eff_syst[2]**2 + etrk_eff_syst[2]**2 + tgtBoil_syst[2]**2)
norm_syst_750_set1 = np.sqrt(htrk_eff_syst[3]**2 + etrk_eff_syst[3]**2 + tgtBoil_syst[3]**2)
norm_syst_750_set2 = np.sqrt(htrk_eff_syst[4]**2 + etrk_eff_syst[4]**2 + tgtBoil_syst[4]**2)
norm_syst_750_set3 = np.sqrt(htrk_eff_syst[5]**2 + etrk_eff_syst[5]**2 + tgtBoil_syst[5]**2)
print('norm_syst_80 = ', norm_syst_80 )
print('norm_syst_580_set1 = ', norm_syst_580_set1 ) 
print('norm_syst_580_set2 = ', norm_syst_580_set2 )       
print('norm_syst_750_set1 = ', norm_syst_750_set1 )                                                                                                                                       
print('norm_syst_750_set2 = ', norm_syst_750_set2 ) 
print('norm_syst_750_set3 = ', norm_syst_750_set3 )                                                                                                                                                     

#Loop over all 2d bins
for ib in range(len(i_b)):

    
    #----------Calculate the Data Reduced Xsec Weighted Average for a given (Pm, th_nq) bin of a given kin. set------------
    red_dataXsec_arr = np.array([sig_red_80[ib], sig_red_580_set1[ib], sig_red_580_set2[ib], sig_red_750_set1[ib], sig_red_750_set2[ib], sig_red_750_set3[ib]])
    red_dataXsec_m = ma.masked_values(red_dataXsec_arr, -1.)  #maske invalid values (-1)
    
    #Define the weights (Use ONLY statistical)
    red_dataXsec_weights = np.array([1./sig_red_80_err[ib]**2, 1./sig_red_580_set1_err[ib]**2, 1./sig_red_580_set2_err[ib]**2, 1./sig_red_750_set1_err[ib]**2, 1./sig_red_750_set2_err[ib]**2, 1./sig_red_750_set3_err[ib]**2])
    red_dataXsec_weights_m = ma.masked_values(red_dataXsec_weights, 1.)


    #=======DATA REDUCED CROSS SECTION WEIGHTED AVERAGE==========
    red_dataXsec_avg = ma.average(red_dataXsec_m, weights=red_dataXsec_weights_m)
    red_dataXsec_avg_err = 1. / np.sqrt(np.sum(red_dataXsec_weights_m))    #combined data sets statistical uncertainty per (Pm, thnq) bin

    #Construct array of kinematic uncertainty for each bin
    kin_syst_arr = np.array([kin_syst_80[ib], kin_syst_580_set1[ib], kin_syst_580_set2[ib], kin_syst_750_set1[ib], kin_syst_750_set2[ib], kin_syst_750_set3[ib]])
    
    #mask elements corresponding to masked dataXsec (We do not want to add in quadrature elements for which there is no Xsec)
    kin_syst_arr_m = ma.masked_array(kin_syst_arr, ma.getmask(red_dataXsec_m))

    #Add masked kinematic systematics array in quadrature
    kin_syst_tot2 = np.sum(kin_syst_arr_m**2)
    kin_syst_tot = np.sqrt( kin_syst_tot2 )


    #March 6, 2020
    #--Construct array of normalization uncertainties--
    norm_syst_arr = np.array([norm_syst_80, norm_syst_580_set1, norm_syst_580_set2, norm_syst_750_set1, norm_syst_750_set2, norm_syst_750_set3])
    #print('norm_syst_arr = ', norm_syst_arr)
    #print('const_norm_syst = ', const_norm_syst)
    #mask elements corresponding to masked dataXsec (We do not want to add in quadrature elements for which there is no Xsec)
    norm_syst_arr_m = ma.masked_array(norm_syst_arr, ma.getmask(red_dataXsec_m))
    #Add masked norm systematics array in quadrature
    norm_syst_final2 = np.sum(norm_syst_arr_m**2 + const_norm_syst**2)
    norm_syst_final = np.sqrt( norm_syst_final2 )

    norm_syst_final2_corr = np.sum(norm_syst_arr_m**2) + const_norm_syst**2
    norm_syst_final_corr = np.sqrt( norm_syst_final2_corr )        

    if norm_syst_final > 0.08:
        print('------WRONG WAY-------------')
        print('norm_syst_arr_m = ',norm_syst_arr_m)
        print('const_norm_syst = ',const_norm_syst)
        print('norm_syst_arr_m**2 + const_norm_syst**2 = ', norm_syst_arr_m**2 + const_norm_syst**2)
        print('norm_syst_final2 = ', norm_syst_final2 )
        print('norm_syst_final = ', norm_syst_final )     

        print('------RIGHT WAY-------------')                                                                                                                                                       
        print('norm_syst_arr_m = ',norm_syst_arr_m)                                                                                                                                                                      
        print('const_norm_syst = ',const_norm_syst)                                                                                                                                                                      
        print('norm_syst_arr_m**2 = ', norm_syst_arr_m**2)                                                                                                                    
        print('np.sum(norm_syst_arr_m**2) + const_norm_syst**2 = ', norm_syst_final2_corr )                                                                                                                          
        print('norm_syst_final_correct = ', norm_syst_final_corr )  

    #print('norm_syst_final=',norm_syst_final)
    #Add total systematic error in quadrature
    tot_syst_err = np.sqrt(kin_syst_tot**2 + norm_syst_final_corr**2)

    #Define relative statistical relative error
    tot_stats_err = red_dataXsec_avg_err/red_dataXsec_avg 
    
    #Add statistical and systematics relative error in quadrature
    tot_err = np.sqrt(tot_stats_err**2 + tot_syst_err**2)
    
    #===Calculate Absolute Error on the Reduced Cross Section===
    red_dataXsec_avg_syst_err = red_dataXsec_avg*tot_syst_err
    red_dataXsec_avg_tot_err = red_dataXsec_avg*tot_err


    #================
    # THEORY AVERAGE
    #================
    
    red_pwiaXsec_arr = np.array([red_pwiaXsec_80[ib], red_pwiaXsec_580_set1[ib], red_pwiaXsec_580_set2[ib], red_pwiaXsec_750_set1[ib], red_pwiaXsec_750_set2[ib], red_pwiaXsec_750_set3[ib]])
    red_pwiaXsec_arr_m = np.ma.masked_values(red_pwiaXsec_arr, -1.)
    red_pwiaXsec_avg = np.ma.average(red_pwiaXsec_arr_m)
   
    red_fsiXsec_arr = np.array([red_fsiXsec_80[ib], red_fsiXsec_580_set1[ib], red_fsiXsec_580_set2[ib], red_fsiXsec_750_set1[ib], red_fsiXsec_750_set2[ib], red_fsiXsec_750_set3[ib]])
    red_fsiXsec_arr_m = np.ma.masked_values(red_fsiXsec_arr, -1.)
    red_fsiXsec_avg = np.ma.average(red_fsiXsec_arr_m)

    #==================================
    # Get the Average Missing Momentum
    #==================================

    pm_arr = np.array([pm80[ib], pm580_set1[ib], pm580_set2[ib], pm750_set1[ib], pm750_set2[ib], pm750_set3[ib]])     
    pm_arr_m = ma.masked_array(pm_arr, ma.getmask(red_dataXsec_m))     #mask elements for which there is no cross section
    pm_avg = np.ma.average(pm_arr_m[pm_arr_m!=0.]) / 1000.  # average :: Convert to GeV  this is the true averaged recoil momentum
                                                                                                                                                         
    
    l="%i  %i  %i  %f   %f   %f   %f  %f   %f   %f   %f   %.12e  %.12e  %.12e  %.12e  %.12e  %.12e\n" % (i_b[ib], i_x[ib], i_y[ib], xb[ib], yb[ib], pm_avg, kin_syst_tot, norm_syst_final_corr, tot_syst_err, tot_stats_err, tot_err, red_dataXsec_avg, red_dataXsec_avg_err, red_dataXsec_avg_syst_err, red_dataXsec_avg_tot_err, red_pwiaXsec_avg, red_fsiXsec_avg)
    fout.write(l)
    
    

fout.close()


