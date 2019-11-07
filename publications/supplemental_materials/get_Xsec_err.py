import numpy as np
from interpolate_ff import *

#This code defined the H(e,e'p) elastic cross section in terms of
#the electric and magnetic form factors. Then apply error propagation
#in terms of the form factors to get uncertainty in cross section at thw
#gven beam energy, elecrtron angle and momentum.
#Finally, use the error in the Xsec as the error as the simc yield error
#Define Physical Constnats


hbarc = 197.327053 / 1000.       #GeV * fm
alpha = 1./137.0359895    #structure constant
dtr = np.pi / 180.

#(1 GeV)^2 = (1000 MeV)^2 = 1e6 MeV^2

#Masses in GeV
MP = 938.272 / 1000.
MN = 939.566 / 1000.
MD = 1875.61 / 1000.
me = 0.51099 / 1000.

#Elastic Cross Section and its error
def sig(kf, th_e, Q2, GEp_GD, GMp_muGD, dGEp_GD, dGMp_muGD):

    #Units: kf ~ Ef [GeV]   electron momentum  
    #th_e[deg] electron angle  
    #Q2[GeV2]   4-momentum transfer

    mu_N = 3.1524512326 *1e-14 / 1000.   #nuclear magneton [GeV T^-1]
    mu_p = 2.7928473508 * mu_N   #proton magnetic moment [GeV T^-1]

    #Mott cross section
    th_e = th_e * dtr    #convert degree to radians
    sigM = (2.*alpha*hbarc*kf*np.cos(th_e/2.)/Q2 )**2   #GeV^2 * fm^2 *GeV^2/GeV^4--> fm2 / sr
    sigMott = sigM * 1e4    #convert fm^2 -> microbarn

    tau =  Q2 /(4.*MP)**2
    L2 = 0.71   #Lambda parameter [GeV^2]
    GD = 1. / (1 + Q2/L2)**2     #dipole form factor
    muGD = mu_p * GD   #[GeV T^-1]
    
    #Get form factors
    GEp = GEp_GD * GD
    GMp = GMp_muGD * muGD

    #Elastic H(e,e'p) Cross Section
    sig = sigMott * ( (GEp**2 + tau*GMp**2)/(1. + tau) + 2.*tau*GMp**2*np.tan(th_e/2.)**2 )

    #Error Propagation
    #derivatives of cross sectio w.r.to form factors
    dsig_dGEp = sigMott * ( (2.*GEp + tau*GMp**2) / (1. + tau))
    dsig_dGMp = sigMott * ( (GEp**2 + 2.*tau*GMp) / (1. + tau) + 4.*tau*GMp*np.tan(th_e/2.)**2 )
    #uncertainty in form factors
    dGEp = dGEp_GD * GD
    dGMp = dGMp_muGD * muGD

    dsig2 =  (dsig_dGEp* dGEp)**2 + (dsig_dGMp*dGMp)**2
    dsig = np.sqrt(dsig2)
    
    return sig, dsig


#Heep Kinematics
dtr = np.pi / 180.
Eb = 10.6005   #beam
Ef = 8.5342488  #e- momentum
th_e = np.array([12.194, 13.93, 9.928, 8.495])  #e- angle
Q2 = 4.*Eb*Ef*np.sin(th_e*dtr/2.)**2   #central Q2 for given kinematics
run = np.array([3288, 3371, 3374, 3377])

for i in range(len(th_e)):
    
    ith_e = th_e[i]
    kf = Ef
    iQ2 = Q2[i] 

    
    GEp_GD = get_GEp_GD(iQ2)
    dGEp_GD = get_dGEp_GD(iQ2)
    GMp_muGD = get_GMp_muGD(iQ2)
    dGMp_muGD = get_dGMp_muGD(iQ2)
     
    Xsec, Xsec_err = sig(kf, ith_e, iQ2, GEp_GD, GMp_muGD, dGEp_GD, dGMp_muGD)

    print('Run Number: %i'%(run[i]))

    print('Beam Energy: ', Eb,' GeV')
    print('electron momentum: ', Ef, ' GeV')
    print('electron angle: ',ith_e, ' deg')
    print('Q2: ', iQ2, ' GeV2')
    print('Xsec = ',Xsec,' +/- ',Xsec_err, 'microbarn/sr')
    print('Xsec_rel_err=',Xsec_err / Xsec)
