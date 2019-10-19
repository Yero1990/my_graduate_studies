import numpy as np
import LT.box as B
import LT.datafile as dfile
import LT_Fit.parameters as P
import LT_Fit.linear_fit as gen_linfit 
#import LT_fit.gen_fit as G

dtr = np.pi / 180.

Mp = 938.272
kin = B.get_file('hms_kin-2.txt')
run = kin['Run']

Eb = np.array(kin['beam_e']) * 1e3
Ef = np.array(kin['nmr_P']) * 1e3
th_e = np.array(kin['hms_Angle']) * dtr 
dataW = np.array(kin['data_W_mean']) * 1e3
simcW = np.array(kin['simc_W_mean']) * 1e3

Q2 = 4.*Eb*Ef*np.sin(th_e/2.)**2
omega = Eb - Ef
W = np.sqrt(Mp**2 + 2.*Mp*omega - Q2)


dW_obs = dataW - simcW
dW_obs_err = Eb * 8000. *2.5e-5 / Ef



#x = np.arange(16) + 1  #array of 1-16
#derivative of W with respect to beam energy
dW_dEb = Ef / Eb *Eb      
#derivative of W with respect to final e- energy
dW_dEf = -Eb / Ef * Ef    
#derivative of W with respect to e- angle
dW_dthe = -2.*Eb*Ef/Mp * np.sin(th_e/2.) * np.cos(th_e/2.)

#Get Initial Relative uncertainties based on dW observations
#These are used to estimate the boundary, or how large should 
#we set the rel. uncertainty. (pmin, pmax) in the chi2Min code
F0 = dW_dEb        #p1 = dEb / Eb    
F1 = dW_dEf 
F2 = dW_dthe    


#--------------------------------------------------------------------------------------------------------------------
#FIT FUNCTION
# dW = a0*dW_dEb + a1*dE_dEf + a2*dW_dthe
#
# the a0, a1 and a2 are the common offsets that get the best fit for all data

#dW/dEb * Eb
def f0(x):
    return F0[np.floor((x+0.5)).astype(int)]

#dW/dE' * E'
def f1(x):
    return F1[np.floor((x+0.5)).astype(int)]

#dW/dth_e 
def f2(x):
    return F2[np.floor((x+0.5)).astype(int)]


#dEm / E'

#Define Initial Parameters
a0 = P.Parameter(0., 'Eb')
a1 = P.Parameter(0., 'Ef')
a2 = P.Parameter(0., 'theta_e')

x = np.arange(len(run))

# select which runs to fit

#sel1 = run < 7000
sel1 = (run >= 3288) & (run<3377)

#sel2 = run > 7000
sel2 = (run >= 6621) & (run<=6632)

#sel3 = (run < 6633) | (7860< run)
sel3 = (run >= 6869) & (run<=6976)

#sel4 = (run >= 6633) & (7860>= run)
sel4 = (run >= 7848) & (run<=7860)

sel5 = (run >= 7962) & (run<=7970)

sel23 = (run >= 6621) & (run<=6976)

sel45 = (run >= 7848) & (run<=7970)

sel_all = (run > 3377)

# set the final selection
sel = sel_all


# fit selected functions and data
#fit = B.gen_linfit([f0, f1, f2], x[sel],  dW_obs[sel], yerr=dW_obs_err[sel], parameters=[a0, a1, a2])
print(dW_obs)
'''
B.plot_exp(x, dW_obs, dW_obs_err)
B.plot_exp(x[sel], fit.func(x[sel]), color = 'r', markersize = 15)
B.pl.xlabel('run index')
B.pl.ylabel('dW diff exp- SIMC')
B.pl.title('')
B.pl.show()

print a0
print a1
print a2

print 'covariance matrix:'
print fit.cov

C = np.zeros_like(fit.cov)
for i,si in enumerate(fit.sig_par):
    for j,sj in enumerate(fit.sig_par):
        sig_i = np.sqrt(fit.cov[i,i])
        sig_j = np.sqrt(fit.cov[j,j])
        C[i,j] = fit.cov[i,j]/(sig_i*sig_j)
print 'Correlation matrix:'
print C
'''

