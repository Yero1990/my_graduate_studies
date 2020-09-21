import LT.box as B
from LT.datafile import dfile
import numpy as np
from deForest_Xsec import *   #import deForest sig_cc1, mott, GMp, GEp functs.

# some constants
dtr = np.pi/180.
hbarc = 197.327053       #MeV * fm

#particle masses (MeV)
MP = 938.272 
MN = 939.566 
MD = 1875.61 
me = 0.51099 

model = "PWIA"

#Digitized data of PWIA 24 calculations from Wally V. Orden (2014 article)
f = dfile('Digitize_Deepn/%s_digitized.data' % (model) )

Pm = np.array(f['xpw'])  #recoil neutron momentum [GeV/c]
theory_Xsec = np.array(f['ypw'])  #5-fold Xsec [GeV^-3 * sr^-2]

#B.plot_exp(x[:13], y[:13], marker='None', linestyle='--', color='b', logy=True)  #lower PWIA
#B.plot_exp(x[14:], y[14:], marker='None', linestyle='--', color='r', logy=True)  #upper PWIA

#----next steps----
#Use E12-10-003 kinematics as states in the W.V. Orden (2014) to evalueate K*f_rec*sig_cc1
#Given the e- kinematics, the proton kinematics is completely defined by the choice of missing momentum
#Eb = 12 GeV, Q2 = 4.25 GeV2, x = 1.35

#-----FROM  ARTICLE-------
Eb = 12.*1000.     #incident energy in MeV
Q2 = 4.25*1e6    #Q2 = 4*Eb*kf*sin^2(theta_e/2) [MeV2]
x = 1.35     #x = Q2 / 2*MD*nu
#-------------------------

nu = Q2 / (2*MP*x)    #nu = Eb - kf, energy transfer
kf = Eb - nu          #final electron energy (or momentum)
theta_e = 2.*np.arcsin( np.sqrt(Q2 / (4.*Eb*kf)) )
th_e = theta_e/dtr  #convert to degrees
q = np.sqrt(Q2 + nu*nu)  #magnitude of 3-momentum transfer


print('Eb = ',Eb)
print('Q2 = ',Q2)
print('x = ',x)
print('nu = ', nu)
print('kf = ',kf)
print('th_e = ', theta_e/dtr)

#--------------Missing Momentum (determine the proton side)------
# Pm2 = (nu + MD - Ep)**2 - MN**2   square of missing momentum
Pr = Pm*1000  #convert Pm from GeV to MeV
Ep = nu + MD - np.sqrt(Pr**2 + MN**2)
Pf = np.sqrt(Ep*Ep - MP*MP)

theta_q = np.arccos((Eb**2 + q**2 - kf**2) / (2.*Eb*q))
theta_pq = np.arccos((Pf**2 + q**2 - Pr**2)/(2.*Pf*q))
#theta_nq = np.arccos( (q**2 + Pr**2 - Pf**2) / (2.*q*Pr) ) 
theta_nq = np.arccos((q - Pf*np.cos(theta_pq))/Pr)
theta_p = theta_q + theta_pq

th_p = theta_p / dtr
th_pq = theta_pq / dtr

ph_pq = 180. #from 2014 article kinematics ( cos(180)=-1 )
cphi_pq = np.cos(180*dtr)
print(cphi_pq)

#Calculate the deForest Cross Section Factor  K * f_rec * sig_cc1,  where K = Pf * Ep , units: ub * MeV^2 / sr^2
sig_Mott =  sigMott(kf, th_e, Q2)   #ub / sr
GE_p = GEp(Q2, 'JRA')  #assumes Q2 is in MeV2
GM_p = GMp(Q2, 'JRA')  #assumes Q2 is in MeV2
Kfact, f_rec, sig_eN, de_Forest = deForest(kf, Q2, q, Pf, Pr, th_e, th_p, cphi_pq, th_pq, sig_Mott, GE_p, GM_p)  #deForest -> K * f_rec * sig_eN [ ub * MeV^2 / sr^2]


#Useful conversion factors
#1 fm2 = 1e4 ub
#1 fm = 1/hbarc MeV^-1
#1 fm3 = (1/hbarc)**3 MeV^-3
#1 GeV3 = 1e9 MeV3

#Convert deForest to fermi
de_Forest_final = de_Forest * (1e-4) * (1./hbarc)**2  #sr^2

#Convert theory Xsec[GeV^-3 * sr^-2] to fm
theory_Xsec_final = theory_Xsec * 1e-9 * (1. / (1/hbarc)**3 )   #fm^3 * sr^2
red_Xsec = theory_Xsec_final / de_Forest_final


#Plot
B.plot_exp(Pm[:19], red_Xsec[:19], marker='None', linestyle='--', color='b', logy=True)  #lower PWIA
B.plot_exp(Pm[20:], red_Xsec[20:], marker='None', linestyle='--', color='r', logy=True)  #upper PWIA
B.pl.show()

#Write to file
fout_name = '%s_digitized_final.txt' % (model)
fout = open(fout_name, 'w')
fout.write('#This file contains the digitized 5-fold diff. cross section, and corresponsing reduced Xsec for the %s calculations plot from W.V.Orden 2014 Deuteron Momentum Dist. Article\n' % (model))
fout.write('#The electron kinematics used were: Eb = 12 GeV, Q2 = 4.25 GeV2, x_bj = 1.35, phi_pq = 180 deg:: theoryXsec[fm^3 sr^2], redXsec[fm^3], K_frec_sigcc1 [sr^2]\n ')
fout.write('#! Pm[f,0]/ theory_Xsec[f,1]/  red_theoryXsec[f,2]/  K_frec_sigcc1[f,3]/\n')

for i in range(len(Pm)):
    line = "{:<15.5f}{:<15.5E}{:<15.5E}{:<15.5E}\n".format( (float(Pm[i])), theory_Xsec_final[i], red_Xsec[i], de_Forest_final[i]) 

    fout.write(line)
fout.close()

