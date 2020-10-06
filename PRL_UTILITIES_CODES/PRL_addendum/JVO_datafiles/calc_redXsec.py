import LT.box as B
from LT.datafile import dfile
import numpy as np
from deForest_Xsec import *   #import deForest sig_cc1, mott, GMp, GEp functs.
from itertools import groupby

# some constants
dtr = np.pi/180.
hbarc = 197.327053       #MeV * fm

#particle masses (MeV)
MP = 938.272 
MN = 939.566 
MD = 1875.61 
me = 0.51099 


thnq = 75
#read data from Wally V. Orden calculations
f = dfile('./%i_deg.data' % (thnq) )

#read avg. kinematics to calculate the recoil factors to determine the reduced cross sections
ib = np.array(f['Nr'])
pm_avg = np.array(f['pm_avg']) * 1000.  #recoil neutron momentum [MeV/c]
Eb = np.array(f['Eb']) #beam energy, MeV
Q2 = np.array(f['Q2']) * 1e6   #Q2, MeV^2
x = np.array(f['x'])   #X-Bjorken
nu = Q2 / (2*MP*x)    #nu = Eb - kf, energy transfer MeV
kf = Eb - nu          #final electron energy (or momentum), MeV
th_e = np.array(f['th_e'])  #e- scattering angle (deg)
q = np.sqrt(Q2 + nu*nu)  #magnitude of 3-momentum transfer, MeV
pf = np.array(f['pf'])    #final proton momentum , MeV
ph_pq = np.array(f['phi'])   #deg

cphi_pq = np.cos(ph_pq * dtr)
theta_q = np.arccos((Eb**2 + q**2 - kf**2) / (2.*Eb*q))   #radians
theta_pq = np.arccos((pf**2 + q**2 - pm_avg**2)/(2.*pf*q))   #radians
theta_nq = np.arccos((q - pf*np.cos(theta_pq))/pm_avg)       #radians
theta_p = theta_q + theta_pq                             #radians

th_p = theta_p / dtr                                    # proton angle (deg)
th_pq = theta_pq / dtr                                  # angle between proton and q-vec (deg)

#Read the theoretical cross sections (nb * MeV^-1 * sr^-2) 
#( nb * MeV^-1 * sr^-2 ) * 0.001 ub/nb * (1. fm^2 / 1.e4 ub) * 1 fm / (1/hbarc MeV^-1) --> fm^3 * sr^-2
sig_exp = np.array(f['sig_exp']) * 0.001 * (1. / 1.e4) * (1.  / (1/hbarc) )  #units of 5-fold Xsec: fm^3 sr^-2 
dsig_exp = np.array(f['dsig_exp']) * 0.001 * (1. / 1.e4) * (1.  / (1/hbarc) ) 
WJC2_GKex05_PWBA = np.array(f['WJC2_GKex05_PWBA']) * 0.001 * (1. / 1.e4) * (1.  / (1/hbarc) ) 
WJC2_GKex05_DWBA = np.array(f['WJC2_GKex05_DWBA']) * 0.001 * (1. / 1.e4) * (1.  / (1/hbarc) ) 
WJC2_AMT_PWBA = np.array(f['WJC2_AMT_PWBA']) * 0.001 * (1. / 1.e4) * (1.  / (1/hbarc) ) 
WJC2_AMT_DWBA = np.array(f['WJC2_AMT_DWBA']) * 0.001 * (1. / 1.e4) * (1.  / (1/hbarc) ) 
AV18_GKex05_PWBA = np.array(f['AV18_GKex05_PWBA']) * 0.001 * (1. / 1.e4) * (1.  / (1/hbarc) ) 
AV18_GKex05_DWBA = np.array(f['AV18_GKex05_DWBA']) * 0.001 * (1. / 1.e4) * (1.  / (1/hbarc) ) 
AV18_AMT_PWBA = np.array(f['AV18_AMT_PWBA']) * 0.001 * (1. / 1.e4) * (1.  / (1/hbarc) ) 
AV18_AMT_DWBA = np.array(f['AV18_AMT_DWBA']) * 0.001 * (1. / 1.e4) * (1.  / (1/hbarc) ) 
CD_GKex05_PWBA = np.array(f['CD_GKex05_PWBA']) * 0.001 * (1. / 1.e4) * (1.  / (1/hbarc) ) 
CD_GKex05_DWBA = np.array(f['CD_GKex05_DWBA']) * 0.001 * (1. / 1.e4) * (1.  / (1/hbarc) ) 
CD_AMT_PWBA = np.array(f['CD_AMT_PWBA']) * 0.001 * (1. / 1.e4) * (1.  / (1/hbarc) ) 
CD_AMT_DWBA = np.array(f['CD_AMT_DWBA']) * 0.001 * (1. / 1.e4) * (1.  / (1/hbarc) ) 

#Write to file
#fout_name = 'JWV_Orden_redXsec_%i_deg.txt' % (thnq)
#fout = open(fout_name, 'w')
#fout.write('#This file contains the reduced cross section from J.W.V.Orden theoretical calculations using some of the models described in 2014 Deuteron Momentum Dist. Article\n')
#fout.write('#The units are: pm_avg(GeV/c), red_Xsec (fm^3)\n')
#fout.write('#! Nr[i,0]/ pm_avg[f,1]/ WJC2_GKex05_PWBA_red[f,2]/ WJC2_GKex05_DWBA_red[f,3]/ WJC2_AMT_PWBA_red[f,4]/ WJC2_AMT_DWBA_red[f,5]/ AV18_GKex05_PWBA_red[f,6]/ AV18_GKex05_DWBA_red[f,7]/  AV18_AMT_PWBA_red[f,8]/ AV18_AMT_DWBA_red[f,9]/  CD_GKex05_PWBA_red[f,10]/ CD_GKex05_DWBA_red[f,11]/ CD_AMT_PWBA_red[f,12]/ CD_AMT_DWBA_red[f,13]/ \n')

#output file to write
fname_out = 'theory_JVO_WJC2_thnq%d.txt' % (thnq)

#open file to write results:
ofile = open(fname_out, 'w')
header = """ 
# Sabine Jeschonnek, J.W.V.Orden, W.P.Ford (JVO WJC2) theoretical cross sections and reduced cross sections as a function of missing momentum
# 
# theta_nq = %d (deg)
#
# header definitions
#
# pm_bin            : missig momentum bin center (GeV/c) (bin width from center is +/- 0.02 GeV)
# pm_avg            : average missing momentum over pm_bin (GeV/c)
# theory_pwbaXsec      : theoretical cross section using the JVO WJC2 PWBA model ( nb / (MeV Sr^2) ) 
# theory_fsiXsec       : theoretical cross section using the JVO WJC2 FSI model ( nb / (MeV Sr^2) )
# theory_red_pwbaXsec  : theoretical reduced cross section using the JVO WJC2 PWBA model (fm^3) 
# theory_red_fsiXsec   : theoretical reduced cross section using the JVO WJC2 FSI model (fm^3)

#! pm_bin[f,0]/ pm_avg[f,1]/ theory_pwbaXsec[f,2]/ theory_fsiXsec[f,3]/  theory_red_pwbaXsec[f,4]/ theory_red_fsiXsec[f,5]/
"""

ofile.write(header)

#Loop over each bin
for i in range(len(ib)):
    
    #Calculate the deForest Cross Section Factor  K * f_rec * sig_cc1,  where K = Pf * Ep , units: ub * MeV^2 / sr^2
    sig_Mott =  sigMott(kf[i], th_e[i], Q2[i])   #ub / sr
    GE_p = GEp(Q2[i], 'JRA')  #assumes Q2 is in MeV2
    GM_p = GMp(Q2[i], 'JRA')  #assumes Q2 is in MeV2

    #All units are in (deg) and (MeV)
    Kfact, f_rec, sig_eN, de_Forest = deForest(kf[i], Q2[i], q[i], pf[i], pm_avg[i], th_e[i], th_p[i], cphi_pq[i], th_pq[i], sig_Mott, GE_p, GM_p)  #deForest -> K * f_rec * sig_eN [ ub * MeV^2 / sr^2]


    #Useful conversion factors
    #1 fm2 = 1e4 ub
    #1 fm = 1/hbarc MeV^-1
    #1 fm3 = (1/hbarc)**3 MeV^-3
    #1 GeV3 = 1e9 MeV3

    #Convert deForest to fermi
    de_Forest_final = de_Forest * (1e-4) * (1./hbarc)**2  #sr^-2

    
    #Calculate Reduced Xsec
    sig_exp_red =           sig_exp[i] / de_Forest_final        
    dsig_exp_red =          dsig_exp[i] / de_Forest_final        
    WJC2_GKex05_PWBA_red =  WJC2_GKex05_PWBA[i] / de_Forest_final 
    WJC2_GKex05_DWBA_red =  WJC2_GKex05_DWBA[i] / de_Forest_final 
    WJC2_AMT_PWBA_red =     WJC2_AMT_PWBA[i] / de_Forest_final   
    WJC2_AMT_DWBA_red =     WJC2_AMT_DWBA[i] / de_Forest_final   
    AV18_GKex05_PWBA_red =  AV18_GKex05_PWBA[i] / de_Forest_final 
    AV18_GKex05_DWBA_red =  AV18_GKex05_DWBA[i] / de_Forest_final 
    AV18_AMT_PWBA_red =     AV18_AMT_PWBA[i] / de_Forest_final   
    AV18_AMT_DWBA_red =     AV18_AMT_DWBA[i] / de_Forest_final   
    CD_GKex05_PWBA_red =    CD_GKex05_PWBA[i] / de_Forest_final  
    CD_GKex05_DWBA_red =    CD_GKex05_DWBA[i] / de_Forest_final  
    CD_AMT_PWBA_red =       CD_AMT_PWBA[i] / de_Forest_final     
    CD_AMT_DWBA_red =       CD_AMT_DWBA[i] / de_Forest_final     
    
    pm = pm_avg[i]/1000.
    
    #fout.write("%i  %f  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E \n" % (ib[i], pm, WJC2_GKex05_PWBA_red, WJC2_GKex05_DWBA_red, WJC2_AMT_PWBA_red, WJC2_AMT_DWBA_red, AV18_GKex05_PWBA_red, AV18_GKex05_DWBA_red, AV18_AMT_PWBA_red, AV18_AMT_DWBA_red, CD_GKex05_PWBA_red, CD_GKex05_DWBA_red, CD_AMT_PWBA_red, CD_AMT_DWBA_red))
    ofile.write("%.5f  %.5f  %.5E  %.5E  %.5E  %.5E \n" % (ib[i], pm,  WJC2_GKex05_PWBA[i], WJC2_GKex05_DWBA[i], WJC2_GKex05_PWBA_red, WJC2_GKex05_DWBA_red))
    
ofile.close()

