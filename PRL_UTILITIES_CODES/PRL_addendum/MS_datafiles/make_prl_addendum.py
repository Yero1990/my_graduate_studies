#This script prepares the PRL addendum data files in the correct formats / units


import LT.box as B
from LT.datafile import dfile
import numpy as np
import matplotlib.pyplot as plt
import numpy.ma as ma
import sys                                     
import os                                                                                                       
from sys import argv  
import matplotlib
from matplotlib import rc
from matplotlib import *

def make_prl_addendum(model='', thnq=0):

    #input file to read
    fname = 'theoryXsec_%sFSI_thnq%.2f.data' % (model, thnq)

    #output file to write
    fname_out = 'theory_MS%s_thnq%d_test.txt' % (model, thnq)
    
    #open file to write results:
    ofile = open(fname_out, 'w')
    header = """ 
    # Misak Sargsian (MS %s) theoretical cross sections and reduced cross sections as a function of missing momentum
    # 
    # theta_nq = %d (deg)
    #
    # header definitions
    #
    # pm_bin            : missig momentum bin center (GeV/c) (bin width from center is +/- 0.02 GeV)
    # pm_avg            : average missing momentum over pm_bin (GeV/c)
    # theory_pwiaXsec      : theoretical cross section using the MS CD-Bonn PWIA model ( nb / (MeV Sr^2) ) 
    # theory_fsiXsec       : theoretical cross section using the MS CD-Bonn FSI model ( nb / (MeV Sr^2) )
    # theory_red_pwiaXsec  : theoretical reduced cross section using the MS CD-Bonn PWIA model (fm^3) 
    # theory_red_fsiXsec   : theoretical reduced cross section using the MS CD-Bonn FSI model (fm^3)
    
    #! pm_bin[f,0]/ pm_avg[f,1]/ theory_pwiaXsec[f,2]/ theory_fsiXsec[f,3]/  theory_red_pwiaXsec[f,4]/ theory_red_fsiXsec[f,5]/  setting[s,6]/ 
    """

    ofile.write(header % (model, thnq))

    # Conversion factor:  1 fm = 1/ (197 MeV),   
    # The reduced cross section is in MeV^-3 units
    MeV2fm = 197.3**3    #convert MeV^-3 to fm^3
    ub2nb = 1000.        # 1 ub == 1000 nb

    #open file and read
    kin = dfile(fname)
    
    #read arrays and find matching thnq bins ONLY
    pm_bin = np.array(kin['pm_bin'])       # GeV/c
    pm_avg = np.array(kin['pm_avg'])       # GeV/c
    theory_pwiaXsec = np.array(kin['pwiaXsec_theory']) * ub2nb              # cross section: nb / (MeV * Sr^2)
    theory_fsiXsec = np.array(kin['fsiXsec_theory']) * ub2nb                # cross section: nb / (MeV * Sr^2)
    theory_red_pwiaXsec = np.array(kin['red_pwiaXsec_theory']) * MeV2fm     # red. Xsec (fm^3)
    theory_red_fsiXsec = np.array(kin['red_fsiXsec_theory']) * MeV2fm       # red. Xsec (fm^3)
    pm_setting = np.array(kin['setting'])
        
    for i, ival in enumerate(pm_bin):
        if(pm_setting[i]=='80'):
            pm_setting[i] = '80_set1'
        ofile.write("%.5f  %.5f  %.5E  %.5E  %.5E  %.5E %s\n" % (pm_bin[i], pm_avg[i], theory_pwiaXsec[i], theory_fsiXsec[i], theory_red_pwiaXsec[i], theory_red_fsiXsec[i], pm_setting[i] ))
        
def main():
    print('Entering Main')

    #make_prl_addendum('CD-Bonn', thnq=35)  #V18 or CD-Bonn
    #make_prl_addendum('CD-Bonn', thnq=45)  #V18 or CD-Bonn
    #make_prl_addendum('CD-Bonn', thnq=75)  #V18 or CD-Bonn

    make_prl_addendum('V18', thnq=35)  #V18 or CD-Bonn
    make_prl_addendum('V18', thnq=45)  #V18 or CD-Bonn
    make_prl_addendum('V18', thnq=75)  #V18 or CD-Bonn

if __name__ == "__main__":
    main()
