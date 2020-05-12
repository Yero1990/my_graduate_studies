#This code reads the theory data files (grouped in th_nq bins) and divided by
#the corresponding Ksig_cc1 to get the reduced Xsec

import LT.box as B
from LT.datafile import dfile
import numpy as np
import sys           
import os                                                      
from sys import argv  


#Conversion Factors
dtr = np.pi / 180.
nb2ub = 1./1000.    # 1ub = 1000 nb
GeV2MeV = 1000.  # 1 GeV = 1000 MeV

#-----------------USER INPUT-------------------------------------
model_dir = "./FSI_models/"    #PWIA_models or FSI_models
#-----------------------------------------------------------------

if(model_dir=="./PWIA_models/"):
    model_name = "2_1_1_0_1"  #2_1_1_0_1 = V18 Potential, 3_1_1_0_1 = CD-Bonn     #!!!!! USER SET BY HAND
    
    #Create Base Name & Directory to Store Update Files
    if(model_name=="2_1_1_0_1"):
        updated_model_name = "V18"
    if(model_name=="3_1_1_0_1"):
        updated_model_name = "CD-Bonn"
    updated_base="theoryXsec_" + updated_model_name + "PWIA"

    output_dir_name = "./updated_PWIA_models/"
    if not os.path.exists(output_dir_name):
        os.mkdir(output_dir_name)

elif(model_dir == "./FSI_models/"):
    model_name = "2_1_1_0_12"  #2_1_1_0_12 = V18 Potential, 3_1_1_0_12 = CD-Bonn,  12 -> PWIA+ FSI    #USER SET BY HAND

    #Create Base Name & Directory to Store Update Files
    if(model_name=="2_1_1_0_12"):
        updated_model_name = "V18"
    if(model_name=="3_1_1_0_12"):
        updated_model_name = "CD-Bonn"
    updated_base="theoryXsec_" + updated_model_name + "FSI"
    output_dir_name = "./updated_FSI_models/"
    if not os.path.exists(output_dir_name):
        os.mkdir(output_dir_name)

base = "csec_calc_theta_pm_"
base2="_fsi_norad_avgkin_"


thnq_arr = np.array([5, 15, 25, 35, 45, 55, 65, 75, 85, 95, 105])

#Loop over theta_nq values (aka distinct theory files)
for i, ithnq in enumerate(thnq_arr):
 
    thnq_f = "%.2f"%(ithnq)  #convert to float

    #OPEND UPDATED THEORY FILES
    updated_theory_fname =  output_dir_name + updated_base + '_thnq%s.data' % (thnq_f)

    #Set header for updated Theory File Names
    fout = open(updated_theory_fname, "w")
    fout.write('#This file contains updated Xsec and Reduced Xsec\n')
    fout.write('#Xsec Units: ub * MeV^-1 *sr^-2\n')
    fout.write('#Ksigcc1 Units: (ub * MeV^2 * sr^-2\n')
    fout.write('#Red Xsec Units: MeV^-3\n')
    fout.write("#!xb[i,0]/  yb[i,1]/   thnq_bin[f,2]/   pm_bin[f,3]/  pm_avg[f,4]/   Ksigcc1[f,5]/   pwiaXsec_theory[f,6]/    fsiXsec_theory[f,7]/   red_pwiaXsec_theory[f,8]/    red_fsiXsec_theory[f,9]/   setting[s,10]/\n")
    fout.close()

    print('Analyzing ',updated_theory_fname)
    
    #Open ORIGINAL THEORY FILES
    theory_fname = model_dir + base +  thnq_f + base2 + model_name + ".data"
    print(theory_fname)

    
    #Open ORIGINAL file and read data into arrays
    kin = dfile(theory_fname)
    pm_setting = kin['kin']
    pm = kin['pm_b']         #central Pmiss value
    pm_avg = kin['pr']         #averaged Pmiss value
    ix = kin['ix']           #x-bin number
    iy = kin['iy']           #y-bin number
    pwia_theoryXsec = kin['crs0']   #crs0: PWIA
    fsi_theoryXsec = kin['crs12']   #crs12: PWIA+FSI
    
    #Loop over all bins of theory datafile
    dlen = len(kin)
    for j in range(dlen):

        if (pm_setting[j]=="80"):
            pm_set = 80
            data_set = 1
        elif (pm_setting[j]=="580_set1"):
            pm_set = 580
            data_set = 1
        elif (pm_setting[j]=="580_set2"):
            pm_set = 580
            data_set = 2
        elif (pm_setting[j]=="750_set1"):
            pm_set = 750
            data_set = 1
        elif (pm_setting[j]=="750_set2"):
            pm_set = 750
            data_set = 2
        elif (pm_setting[j]=="750_set3"):
            pm_set = 750
            data_set = 3

        #Open the corresponding file with the K*sigcc1, and match (ix,iy) bins before dividing by K*sigcc1
        avg_kin_dir = "../average_kinematics/Em_final40MeV/"
        if(model_dir=="./PWIA_models/"):
            if(pm_set==80):
                avg_kin_fname = avg_kin_dir + 'pm80_pwia_norad_avgkin.txt'
            else:
                avg_kin_fname = avg_kin_dir + 'pm%i_pwia_norad_avgkin_set%i.txt' % (pm_set, data_set)
        elif(model_dir=="./FSI_models/"):
            if(pm_set==80):
                avg_kin_fname = avg_kin_dir + 'pm80_fsi_norad_avgkin.txt'
            else:
                avg_kin_fname = avg_kin_dir + 'pm%i_fsi_norad_avgkin_set%i.txt' % (pm_set, data_set)

        avg_kin = dfile(avg_kin_fname)
        dklen = len(avg_kin['i_b'])  #get length of array of avgkin file
        pm_k = avg_kin['yb']
        pm_k_avg = avg_kin['pm'] / 1000.   #convert to GeV/c
        thnq_k = avg_kin['xb']
        ix_k = avg_kin['i_x']
        iy_k = avg_kin['i_y']
        Ksig_cc1 = avg_kin['Ksig_cc1']

        #Loop over all bins in avg. kin file to find matching bins between avg kin file and theory
        for k in range(dklen):
            if(ix_k[k]==ix[j] and iy_k[k]==iy[j]):
                #print('Found Matching Bins')
                
                #Debug/TEST
                #print('datafile: pm=',pm[j],' thnq=',ithnq)
                #print('avgkin_datafile: pm=',pm_k[k],' thnq=',thnq_k[k])
                #print('KSig_cc1 = ', Ksig_cc1[k])
                #Do the Calculations Here ! ! !
                if(Ksig_cc1[k]==0.0):
                    continue
                #Convert Xsec units from nb * GeV ^-1 * sr^-2 to ub * MeV^-1 *sr^-2
                pwia_theoryXsec[j] = pwia_theoryXsec[j] * nb2ub * (1./GeV2MeV)
                fsi_theoryXsec[j]  = fsi_theoryXsec[j] * nb2ub * (1./GeV2MeV)
                
                #Calculate Reduced Xsec
                red_pwiaXsec = pwia_theoryXsec[j] / Ksig_cc1[k]
                red_fsiXsec = fsi_theoryXsec[j] / Ksig_cc1[k]
                
                fout = open(updated_theory_fname, "a")   #open in append mode
                
                #Update Theory Files
                fout.write("%i   %i   %f   %f   %f   %f  %.12e   %.12e    %.12e   %.12e   %s\n"%(ix_k[k], iy_k[k], thnq_k[k], pm_k[k], pm_k_avg[k],  Ksig_cc1[k], pwia_theoryXsec[j], fsi_theoryXsec[j], red_pwiaXsec, red_fsiXsec, pm_setting[j]))


