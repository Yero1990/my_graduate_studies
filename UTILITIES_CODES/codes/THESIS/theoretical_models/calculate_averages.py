import numpy as np
from LT.datafile import dfile
import averages as av

#This code uses the averages utilities code to calculate averages based on grouped values in a given array input
#It is assumed that the reduced Xsec have been calculated, so their average can be taken.

#----------USER INPUT-------------
theory = "V18"   #V18 or CD-Bonn
model = "FSI"   #FSI or PWIA

thnq_arr = np.array([5, 15, 25, 35, 45, 55, 65, 75, 85, 95, 105])

#Loop over theta_nq values (aka distinct theory files)
for i, ithnq in enumerate(thnq_arr):
 
    thnq_f = "%.2f"%(ithnq)  #convert to float

    print('thnq=',thnq_f)

    #Set Theory Filename to be read
    fname = './updated_%s_models/theoryXsec_%s%s_thnq%s.data' %(model, theory, model, thnq_f )
    kin = dfile(fname)

    #Get arrays to take average
    pm_bin = np.array(kin['pm_bin'])
    pm_avg = np.array(kin['pm_avg'])
    thnq_bin = np.array(kin['thnq_bin'])
    xb = np.array(kin['xb'])
    yb = np.array(kin['yb'])
    red_pwiaXsec = np.array(kin['red_pwiaXsec_theory'])
    red_fsiXsec = np.array(kin['red_fsiXsec_theory'])
    Xsec_err = np.ones(len(kin['xb']))   #Set error to 1.0 as this is used for the weight 1/sig**2, which should be 1 (or no preference)


    #pm_avg1, red_pwiaXsec_avg, dsig_avg1, dsig1_avg = av.get_matched_array_average(pm_bin, red_pwiaXsec, Xsec_err, same = True, match = None)
    #pm_avg2, red_fsiXsec_avg, dsig_avg2, dsig2_avg = av.get_matched_array_average(pm_bin, red_fsiXsec, Xsec_err, same = True, match = None)
    
    pm_avg1, red_pwiaXsec_avg, dsig_avg1, dsig1_avg = av.get_matched_array_average(pm_avg, red_pwiaXsec, Xsec_err, same = False, match = pm_bin)
    pm_avg2, red_fsiXsec_avg, dsig_avg2, dsig2_avg = av.get_matched_array_average(pm_avg, red_fsiXsec, Xsec_err, same = False, match = pm_bin)


    #Create Output file name
    fout_name = './updated_%s_models/theoryXsec_%s%s_thnq%s_combined.data' %(model, theory, model, thnq_f )
    fout = open(fout_name, "w")
    fout.write('#This file contains combined (averaged) Reduced Xsec\n')
    fout.write('#Red Xsec Units: MeV^-3\n')
    if model=="PWIA":
        fout.write("#!thnq_bin[f,0]/  pm_avg[f,1]/   red_pwiaXsec_theory[f,2]/\n")
        for j in range(len(pm_avg1)):
            fout.write('%f   %f   %.12e\n' % (thnq_bin[j], pm_avg1[j], red_pwiaXsec_avg[j]))
    else:
        fout.write("#!thnq_bin[f,0]/  pm_avg[f,1]/   red_fsiXsec_theory[f,2]/\n")
        for j in range(len(pm_avg1)):
            fout.write('%f   %f   %.12e\n' % (thnq_bin[j], pm_avg1[j], red_fsiXsec_avg[j]))

    
    
    fout.close()
