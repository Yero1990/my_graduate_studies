import numpy as np
from LT.datafile import dfile
import LT.box as B
import shutil as sh
import os

#This code has utilities functions that calculates the systematics 
#due to the correction factors in the data yield for each data set

#If, FullWeight = f_tgt * f_pAbs * f_tLT  * f_eTrkEff * f_pTrkEff * f_Qtot
#Experimental Xsec: sig_exp = Y_corr / P.S. = (Y_uncorr / P.S.) *  f_tgt * f_pAbs * f_tLT  * f_eTrkEff * f_pTrkEff * f_Qtot

#Each of the corrections has an associated uncertainty, df/f, therefore, the
#systematic error on the cross section (sig_exp ) is:  dsig_exp_syst2 = [sig_exp**2 * (df/f)**2 ]_tgtBoil + . . . +  [sig_exp**2 * (df/f)**2 ]_pAbs + ...
#A derivative has to be takem to determine the uncertainty in f.

#***IMPORTANT***: Each of the utilities functions saves a copy of the average kinematics file with the dsig_i (contributions to the cross section systematics error)
#                 The errors are written to file as relative error, df/f. To get the absolute error, one would multiply by the cross section to get: df = (df/f)*f
#                 A special utility function calculates the total dsig2 from the sum in quadrature of the norm systematics: dsig2_norm_tot = dsig2_hEtrk + dsig2_Qtot + dsig2_tLT + . . .
#                 The utilities functions also return their respective the dsig_2 contributions as well as the relative error df/f, which gives an idea of how large is the error relative to the Xsec
                  #(The model='' input parameter refers to the model used for the radiative corrections. We use Laget FSI, as it closely represents data in all kinematical regions)

#Tracking efficiencies
def get_trkEff_syst(pm_set=0, model='', data_set=0):
    
    #This code reads the tracking efficiency and its uncertainty
    #for all runs of a given data set, and calculates the weighted 
    #average. The systematic uncertanity on the data Xsec is determined from the relative
    #uncertainty on the tracking efficiencies.
    
    report_file = get_filenames(pm_set, model, data_set)

    
    #open report 
    r = dfile(report_file)                 #open report file
    
    #Get Tracking Eff. and its error
    htrk_eff = np.array(r['hTrkEff'])
    htrk_eff_err = np.array(r['hTrkEff_err'])
    
    etrk_eff = np.array(r['eTrkEff'])
    etrk_eff_err = np.array(r['eTrkEff_err'])

    #Get Average Correction Factors and Average Error
    htrk_eff_avg = np.average(htrk_eff[htrk_eff_err!=0])
    htrk_eff_avg_err = np.average(htrk_eff_err[htrk_eff_err!=0])
    
    etrk_eff_avg = np.average(etrk_eff[etrk_eff_err!=0])
    etrk_eff_avg_err = np.average(etrk_eff_err[etrk_eff_err!=0])
    
    #average of relative systematic error (per data set)
    df_f_htrk_avg = htrk_eff_avg_err / htrk_eff_avg
    df_f_etrk_avg = etrk_eff_avg_err / etrk_eff_avg
  


    return[htrk_eff_avg, htrk_eff_avg_err, etrk_eff_avg, etrk_eff_avg_err, df_f_htrk_avg, df_f_etrk_avg]

def get_tgtboil_syst(pm_set=0, model='', data_set=0):

    #This code calculated the systematic uncertainty on the data Xsec from the uncertainty in
    #the target boiling factor.

    
    report_file = get_filenames(pm_set, model, data_set)

    #open report 
    r = dfile(report_file)                 #open report file

    I = np.array(r['avg_current'])     #in uA  (is an array if there are multiple runs per data set)
    tgt_boil = np.array(r['tgtBoil_factor'])

    #Code to determine the target boiling factor systematic effect on the cross section.
    m =  0.00080029          #LD2 slope from fit [fract. yield loss / uA]
    sig_m = 0.00007037       #uncertainty for LD2 target boiling fit slope
    dI_I = 2.0/100.          #relative uncertainty in average current ( 2 % to be conservative. Suggested by D. Mack)

    sig_I = dI_I * I             

    #fcorr = (1. - m*I)       #equivalent to tgt_boil 
    tgt_boil_err = np.sqrt((I*sig_m)**2 + (m*sig_I)**2)

    #Get Average tgt_boil factor and error
    tgt_boil_avg = np.average(tgt_boil)
    tgt_boil_avg_err = np.average(tgt_boil_err)
    
    
    #average of relative systematic error (per data set)
    df_f_avg = tgt_boil_avg_err / tgt_boil_avg
    
    return [tgt_boil_avg, tgt_boil_avg_err, df_f_avg]

def get_pT_syst(pm_set=0, model='', data_set=0):
    
    #code to get systematics on proton absorption
    #The proton absorption for E12-10-003 was determined to 
    #be: 4.66 +/- 0.472 %,  the transmission is: (100-4.66 or 95.34) +/- 0.472%
    #It is assumed that this factor is relatively constant over a range of angles/momenta

    report_file = get_filenames(pm_set, model, data_set)

    #open report 
    r = dfile(report_file)                 #open report file
    
    pT = 0.9534           #proton transmission factor (what fraction of coin. protons made it to form the trigger) --See docDB proton absorption Under Yero.
    d_pT = 0.472/100.     

    
    #average of relative systematic error (per data set)
    df_f_avg = d_pT / pT

    return [pT, d_pT, df_f_avg]
    
def get_tLT_syst(pm_set=0, model='', data_set=0):

    #code to get the total live time (from EDTM) systematics
    
    report_file = get_filenames(pm_set, model, data_set)

    #open report 
    r = dfile(report_file)                 #open report file

    tLT = np.array(r['tLT'])     #total live time [fraction or %/100]
    shms_3of4 = np.array(r['ptrig1_rate'])    #shms 3/4 trigger tae in kHz.  Used to crudely estimate relative error on Xsec
    
    df_f_approx = shms_3of4 * 1000. * (100 * 1e-9)   #this is prob. of how many 3/4 triggers would be blocked given 100 ns window.
    df_f_approx_avg = np.average(df_f_approx)

    print('total_live_time_rel_err=',df_f_approx_avg)

    d_tLT = 0.03 * tLT          #assume uncertainty in total live time contributes to 3 % systematics (conservative estimate)

    tLT_avg = np.average(tLT)
    tLT_avg_err = np.average(d_tLT)
    
    
    #average of relative systematic error (per data set)
    df_f_avg = tLT_avg_err / tLT_avg

    return [tLT_avg, tLT_avg_err, round(df_f_avg,4)]
    
def get_Qtot_syst(pm_set=0, model='', data_set=0):
    
    #code to get the systematics of the total accumulated charge (related to BCMs uncertainty)
      
    report_file = get_filenames(pm_set, model, data_set)

    #open report 
    r = dfile(report_file)                 #open report file

    Q = np.array(r['charge'])          #total accumulated charge (mC)
    dQ_Q = 0.02           #uncertainty in total charge (take 2% for now, conservative)

    dQ = dQ_Q * Q

    Q_tot = Q.sum()
    dQ_tot = np.sum(dQ)
    
    
    #average of relative systematic error (per data set)
    df_f_avg = dQ_tot/Q_tot    #ASSUME  2% for now, conservative)

    return [Q_tot, dQ_tot, round(df_f_avg, 4)]


def convert2NaN(arr=np.array([]), value=0):
    #method to convert a specified value in a array to nan (not a number)
    
    for i in enumerate(arr):
        if arr[i[0]]==value:
            arr[i[0]] = np.nan
    return arr


def get_filenames(pm_set, model, data_set):


    report_file='report_deep_pm%i_set%i.dat'%(pm_set, data_set)


    return report_file


def get_average_systematics():
    
    #This code calls functions to get systematics (that have been already averaged per run)
    #and averages them over different data sets to get an overall final systeamtic value

    #Pm = 80 (set1)
    htrk_eff_avg80, htrk_eff_avg_err80, etrk_eff_avg80, etrk_eff_avg_err80, df_f_htrk_avg80, df_f_etrk_avg80 = get_trkEff_syst(80, 'fsi', 1)
    tgt_boil_avg80, tgt_boil_avg_err80, df_f_avg_boil80 = get_tgtboil_syst(80, 'fsi', 1)
    pT_80, d_pT_80, df_f_avg_pT80 = get_pT_syst(80, 'fsi', 1)
    tLT_80, d_tLT_80, df_f_avg_tLT80 = get_tLT_syst(80, 'fsi', 1)
    Q_tot_80, dQ_tot_80, df_f_avg_Q80 = get_Qtot_syst(80, 'fsi', 1)

    #Pm = 580 (set1)
    htrk_eff_avg580s1, htrk_eff_avg_err580s1, etrk_eff_avg580s1, etrk_eff_avg_err580s1, df_f_htrk_avg580s1, df_f_etrk_avg580s1 = get_trkEff_syst(580, 'fsi', 1)
    tgt_boil_avg580s1, tgt_boil_avg_err580s1, df_f_avg_boil580s1 = get_tgtboil_syst(580, 'fsi', 1)
    pT_580s1, d_pT_580s1, df_f_avg_pT580s1 = get_pT_syst(580, 'fsi', 1)
    tLT_580s1, d_tLT_580s1, df_f_avg_tLT580s1 = get_tLT_syst(580, 'fsi', 1)
    Q_tot_580s1, dQ_tot_580s1, df_f_avg_Q580s1 = get_Qtot_syst(580, 'fsi', 1)

    #Pm = 580 (set2)
    htrk_eff_avg580s2, htrk_eff_avg_err580s2, etrk_eff_avg580s2, etrk_eff_avg_err580s2, df_f_htrk_avg580s2, df_f_etrk_avg580s2 = get_trkEff_syst(580, 'fsi', 2)
    tgt_boil_avg580s2, tgt_boil_avg_err580s2, df_f_avg_boil580s2 = get_tgtboil_syst(580, 'fsi', 2)
    pT_580s2, d_pT_580s2, df_f_avg_pT580s2 = get_pT_syst(580, 'fsi', 2)
    tLT_580s2, d_tLT_580s2, df_f_avg_tLT580s2 = get_tLT_syst(580, 'fsi', 2)
    Q_tot_580s2, dQ_tot_580s2, df_f_avg_Q580s2 = get_Qtot_syst(580, 'fsi', 2)

    #Pm = 750 (set1)
    htrk_eff_avg750s1, htrk_eff_avg_err750s1, etrk_eff_avg750s1, etrk_eff_avg_err750s1, df_f_htrk_avg750s1, df_f_etrk_avg750s1 = get_trkEff_syst(750, 'fsi', 1)
    tgt_boil_avg750s1, tgt_boil_avg_err750s1, df_f_avg_boil750s1 = get_tgtboil_syst(750, 'fsi', 1)
    pT_750s1, d_pT_750s1, df_f_avg_pT750s1 = get_pT_syst(750, 'fsi', 1)
    tLT_750s1, d_tLT_750s1, df_f_avg_tLT750s1 = get_tLT_syst(750, 'fsi', 1)
    Q_tot_750s1, dQ_tot_750s1, df_f_avg_Q750s1 = get_Qtot_syst(750, 'fsi', 1)

    #Pm = 750 (set2)
    htrk_eff_avg750s2, htrk_eff_avg_err750s2, etrk_eff_avg750s2, etrk_eff_avg_err750s2, df_f_htrk_avg750s2, df_f_etrk_avg750s2 = get_trkEff_syst(750, 'fsi', 2)
    tgt_boil_avg750s2, tgt_boil_avg_err750s2, df_f_avg_boil750s2 = get_tgtboil_syst(750, 'fsi', 2)
    pT_750s2, d_pT_750s2, df_f_avg_pT750s2 = get_pT_syst(750, 'fsi', 2)
    tLT_750s2, d_tLT_750s2, df_f_avg_tLT750s2 = get_tLT_syst(750, 'fsi', 2)
    Q_tot_750s2, dQ_tot_750s2, df_f_avg_Q750s2 = get_Qtot_syst(750, 'fsi', 2)

    #Pm = 750 (set3)
    htrk_eff_avg750s3, htrk_eff_avg_err750s3, etrk_eff_avg750s3, etrk_eff_avg_err750s3, df_f_htrk_avg750s3, df_f_etrk_avg750s3 = get_trkEff_syst(750, 'fsi', 3)
    tgt_boil_avg750s3, tgt_boil_avg_err750s3, df_f_avg_boil750s3 = get_tgtboil_syst(750, 'fsi', 3)
    pT_750s3, d_pT_750s3, df_f_avg_pT750s3 = get_pT_syst(750, 'fsi', 3)
    tLT_750s3, d_tLT_750s3, df_f_avg_tLT750s3 = get_tLT_syst(750, 'fsi', 3)
    Q_tot_750s3, dQ_tot_750s3, df_f_avg_Q750s3 = get_Qtot_syst(750, 'fsi', 3)
    

    #Create array of relative error, to be ablw to write to file, as well as take an overall average  
    df_f_htrkEff_arr = np.array([df_f_htrk_avg80, df_f_htrk_avg580s1, df_f_htrk_avg580s2, df_f_htrk_avg750s1, df_f_htrk_avg750s2, df_f_htrk_avg750s3])
    df_f_etrkEff_arr = np.array([df_f_etrk_avg80, df_f_etrk_avg580s1, df_f_etrk_avg580s2, df_f_etrk_avg750s1, df_f_etrk_avg750s2, df_f_etrk_avg750s3])
    df_f_tgtBoil_arr = np.array([df_f_avg_boil80, df_f_avg_boil580s1, df_f_avg_boil580s2, df_f_avg_boil750s1, df_f_avg_boil750s2, df_f_avg_boil750s3])
    df_f_pT_arr = np.array([df_f_avg_pT80, df_f_avg_pT580s1, df_f_avg_pT580s2, df_f_avg_pT750s1, df_f_avg_pT750s2, df_f_avg_pT750s3])
    df_f_tLT_arr = np.array([df_f_avg_tLT80, df_f_avg_tLT580s1, df_f_avg_tLT580s2, df_f_avg_tLT750s1, df_f_avg_tLT750s2, df_f_avg_tLT750s3])
    df_f_Qtot_arr = np.array([df_f_avg_Q80, df_f_avg_Q580s1, df_f_avg_Q580s2, df_f_avg_Q750s1, df_f_avg_Q750s2, df_f_avg_Q750s3])

    #Average of relative systematics (df/f) on correction factors over all momentum settings (to quote on PRL).
    df_f_htrkEff_avg = np.average(df_f_htrkEff_arr)
    df_f_etrkEff_avg = np.average(df_f_etrkEff_arr)
    df_f_tgtBoil_avg = np.average(df_f_tgtBoil_arr)
    df_f_pT_avg = np.average(df_f_pT_arr)
    df_f_tLT_avg = np.average(df_f_tLT_arr)
    df_f_Qtot_avg = np.average(df_f_Qtot_arr)

    #Add (averaged syst. err over all settings) quadratically (to quote on PRL)
    df_f_total_avg = np.sqrt( df_f_htrkEff_avg**2 +  df_f_etrkEff_avg**2 + df_f_tgtBoil_avg**2 + df_f_pT_avg**2 + df_f_tLT_avg**2 + df_f_Qtot_avg**2 )

    
    #Find averages of correction factors over all data sets (to quote on PRL)
    htrk_eff_avg = (htrk_eff_avg80 + htrk_eff_avg580s1 + htrk_eff_avg580s2 + htrk_eff_avg750s1 + htrk_eff_avg750s2 + htrk_eff_avg750s3) / 6.
    etrk_eff_avg = (etrk_eff_avg80 + etrk_eff_avg580s1 + etrk_eff_avg580s2 + etrk_eff_avg750s1 + etrk_eff_avg750s2 + etrk_eff_avg750s3) / 6.
    tgt_boil_avg = (tgt_boil_avg80 + tgt_boil_avg580s1 + tgt_boil_avg580s2 + tgt_boil_avg750s1 + tgt_boil_avg750s2 + tgt_boil_avg750s3) / 6.
    pT_avg = (pT_80 + pT_580s1 + pT_580s2 + pT_750s1 + pT_750s2 +pT_750s3) / 6.
    tLT_avg = (tLT_80 + tLT_580s1 + tLT_580s2 + tLT_750s1 + tLT_750s2 + tLT_750s3) / 6.
    Q_tot = Q_tot_80 + Q_tot_580s1 + Q_tot_580s2 + Q_tot_750s1 + Q_tot_750s2 + Q_tot_750s3

    htrk_eff_avg_err = (htrk_eff_avg_err80 + htrk_eff_avg_err580s1 + htrk_eff_avg_err580s2 + htrk_eff_avg_err750s1 + htrk_eff_avg_err750s2 + htrk_eff_avg_err750s3) / 6.
    etrk_eff_avg_err = (etrk_eff_avg_err80 + etrk_eff_avg_err580s1 + etrk_eff_avg_err580s2 + etrk_eff_avg_err750s1 + etrk_eff_avg_err750s2 + etrk_eff_avg_err750s3) / 6.
    tgt_boil_avg_err = (tgt_boil_avg_err80 + tgt_boil_avg_err580s1 + tgt_boil_avg_err580s2 + tgt_boil_avg_err750s1 + tgt_boil_avg_err750s2 + tgt_boil_avg_err750s3) / 6.
    pT_avg_err = (d_pT_80 + d_pT_580s1 + d_pT_580s2 + d_pT_750s1 + d_pT_750s2 + d_pT_750s3) / 6.
    tLT_avg_err = (d_tLT_80 + d_tLT_580s1 + d_tLT_580s2 + d_tLT_750s1 + d_tLT_750s2 + d_tLT_750s3) / 6.
    Q_tot_err = (dQ_tot_80 + dQ_tot_580s1 + dQ_tot_580s2 + dQ_tot_750s1 + dQ_tot_750s2 + dQ_tot_750s3) / 6.


    #---Add normalization syst. err quadratically per data set (to be added to statistical error later on)
    df_f_80 = np.sqrt( df_f_htrkEff_arr[0]**2 + df_f_etrkEff_arr[0]**2 + df_f_tgtBoil_arr[0]**2 + df_f_pT_arr[0]**2 + df_f_tLT_arr[0]**2 + df_f_Qtot_arr[0]**2)
    df_f_580s1 = np.sqrt( df_f_htrkEff_arr[1]**2 + df_f_etrkEff_arr[1]**2 + df_f_tgtBoil_arr[1]**2 + df_f_pT_arr[1]**2 + df_f_tLT_arr[1]**2 + df_f_Qtot_arr[1]**2)
    df_f_580s2 = np.sqrt( df_f_htrkEff_arr[2]**2 + df_f_etrkEff_arr[2]**2 + df_f_tgtBoil_arr[2]**2 + df_f_pT_arr[2]**2 + df_f_tLT_arr[2]**2 + df_f_Qtot_arr[2]**2)
    df_f_750s1 = np.sqrt( df_f_htrkEff_arr[3]**2 + df_f_etrkEff_arr[3]**2 + df_f_tgtBoil_arr[3]**2 + df_f_pT_arr[3]**2 + df_f_tLT_arr[3]**2 + df_f_Qtot_arr[3]**2)
    df_f_750s2 = np.sqrt( df_f_htrkEff_arr[4]**2 + df_f_etrkEff_arr[4]**2 + df_f_tgtBoil_arr[4]**2 + df_f_pT_arr[4]**2 + df_f_tLT_arr[4]**2 + df_f_Qtot_arr[4]**2)
    df_f_750s3 = np.sqrt( df_f_htrkEff_arr[5]**2 + df_f_etrkEff_arr[5]**2 + df_f_tgtBoil_arr[5]**2 + df_f_pT_arr[5]**2 + df_f_tLT_arr[5]**2 + df_f_Qtot_arr[5]**2)


    
    fout_path='./normalization_systematics_summary.txt'
    fout = open(fout_path, 'w')
    fout.write('#This file contains the run-by-run averaged (for all data sets) correction factor and systematics on the cross section\n')
    fout.write('#These are relative errors, dsig / dig. Units: Multiply by 100 to convert to %%\n')
    fout.write('The last header is the quadrature sum of the headers, 0-6.  \n')
    fout.write('#! pm[i,0]/  set[i,1]/   htrk_eff[f,2]/  htrk_eff_err[f,3]/ etrk_eff[f,4]/ etrk_eff_err[f,5]/  tgtBoil[f,6]/ tgtBoil_err[f,7]/  pAbs[f,8]/   pAbs_err[f,9]/  tLT[f,10]/   tLT_err[f,11]/  Qtot[f,12]/   Qtot_err[f,13]/  htrk_eff_syst[f,14]/  etrk_eff_syst[f,15]/  tgtBoil_syst[f,16]/  pAbs_syst[f,17]/  tLT_syst[f,18]/   Qtot_syst[f,19]/  total_norm_syst[f,20]/ \n')
    fout.write('%i  %i  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f %f  %f  %f  %f  %f  %f  %f \n'%(80, 1, htrk_eff_avg80, htrk_eff_avg_err80, etrk_eff_avg80, etrk_eff_avg_err80, tgt_boil_avg80, tgt_boil_avg_err80, pT_80, d_pT_80, tLT_80, d_tLT_80, Q_tot_80, dQ_tot_80, df_f_htrk_avg80, df_f_etrk_avg80, df_f_avg_boil80, df_f_avg_pT80, df_f_avg_tLT80, df_f_avg_Q80, df_f_80))
    fout.write('%i  %i  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f %f  %f  %f  %f  %f  %f  %f \n'%(580, 1, htrk_eff_avg580s1, htrk_eff_avg_err580s1, etrk_eff_avg580s1, etrk_eff_avg_err580s1, tgt_boil_avg580s1, tgt_boil_avg_err580s1, pT_580s1, d_pT_580s1, tLT_580s1, d_tLT_580s1, Q_tot_580s1, dQ_tot_580s1, df_f_htrk_avg580s1, df_f_etrk_avg580s1, df_f_avg_boil580s1, df_f_avg_pT580s1, df_f_avg_tLT580s1, df_f_avg_Q580s1, df_f_580s1))
    fout.write('%i  %i  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f %f  %f  %f  %f  %f  %f  %f \n'%(580, 2, htrk_eff_avg580s2, htrk_eff_avg_err580s2, etrk_eff_avg580s2, etrk_eff_avg_err580s2, tgt_boil_avg580s2, tgt_boil_avg_err580s2, pT_580s2, d_pT_580s2, tLT_580s2, d_tLT_580s2, Q_tot_580s2, dQ_tot_580s2, df_f_htrk_avg580s2, df_f_etrk_avg580s2, df_f_avg_boil580s2, df_f_avg_pT580s2, df_f_avg_tLT580s2, df_f_avg_Q580s2, df_f_580s2))
    fout.write('%i  %i  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f %f  %f  %f  %f  %f  %f  %f \n'%(750, 1, htrk_eff_avg750s1, htrk_eff_avg_err750s1, etrk_eff_avg750s1, etrk_eff_avg_err750s1, tgt_boil_avg750s1, tgt_boil_avg_err750s1, pT_750s1, d_pT_750s1, tLT_750s1, d_tLT_750s1, Q_tot_750s1, dQ_tot_750s1, df_f_htrk_avg750s1, df_f_etrk_avg750s1, df_f_avg_boil750s1, df_f_avg_pT750s1, df_f_avg_tLT750s1, df_f_avg_Q750s1, df_f_750s1))
    fout.write('%i  %i  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f %f  %f  %f  %f  %f  %f  %f \n'%(750, 2, htrk_eff_avg750s2, htrk_eff_avg_err750s2, etrk_eff_avg750s2, etrk_eff_avg_err750s2, tgt_boil_avg750s2, tgt_boil_avg_err750s2, pT_750s2, d_pT_750s2, tLT_750s2, d_tLT_750s2, Q_tot_750s2, dQ_tot_750s2, df_f_htrk_avg750s2, df_f_etrk_avg750s2, df_f_avg_boil750s2, df_f_avg_pT750s2, df_f_avg_tLT750s2, df_f_avg_Q750s2, df_f_750s2))
    fout.write('%i  %i  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f %f  %f  %f  %f  %f  %f  %f \n'%(750, 3, htrk_eff_avg750s3, htrk_eff_avg_err750s3, etrk_eff_avg750s3, etrk_eff_avg_err750s3, tgt_boil_avg750s3, tgt_boil_avg_err750s3, pT_750s3, d_pT_750s3, tLT_750s3, d_tLT_750s3, Q_tot_750s3, dQ_tot_750s3, df_f_htrk_avg750s3, df_f_etrk_avg750s3, df_f_avg_boil750s3, df_f_avg_pT750s3, df_f_avg_tLT750s3, df_f_avg_Q750s3, df_f_750s3))
    fout.write('%i  %i  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f %f  %f  %f  %f  %f  %f  %f \n'%(-1, -1, htrk_eff_avg,      htrk_eff_avg_err,      etrk_eff_avg,      etrk_eff_avg_err,      tgt_boil_avg,      tgt_boil_avg_err,      pT_avg,   pT_avg_err, tLT_avg,   tLT_avg_err, Q_tot,       Q_tot_err,    df_f_htrkEff_avg,   df_f_etrkEff_avg,   df_f_tgtBoil_avg,   df_f_pT_avg,      df_f_tLT_avg,      df_f_Qtot_avg,   df_f_total_avg))
    fout.close()

    


def main():

    print('Main')
    #Call method to calculate and combine normalization systematic errors

    #combine_norm_systematics(80, 'fsi', 1)
    #get_radbc_syst(fname, 80, 1)
    #plot_relative_error(fname, 80, 1)
    #plot_Xsec_Ratio(fname, 80, 1)
    get_average_systematics()

if __name__=="__main__":
    main()

