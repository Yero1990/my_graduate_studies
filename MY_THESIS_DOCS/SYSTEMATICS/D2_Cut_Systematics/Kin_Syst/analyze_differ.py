import numpy  as np
dtr = np.pi/180.


def get_sig_tot(differ_file='',  ftable_name='', ftable_bin = np.array([]), print_all = False):
    data = open(differ_file).readlines()
    found = 0
    
    # This is where you enter the estimated uncertainties
    # sigma in mr 

    #Values were taken from SHMS Coin. Chi2 Min. MODEL 01, where Eb, Ef, th_e, th_f are determined simultaneously
    #electron arm uncertainties
    sig_the = 0.1659    #electron theta_e (in-plane, horizontal) uncertainty in mr
    sig_phe = 0.    #electron phi_e (out-of-plane, vertical) uncertainty in mr
    
    #proton (or hadron) arm uncertainties
    sig_thp = 0.2369    #proton theta_p (in-plane, horizontal) uncertainty mr 
    sig_php = 0.    #proton phi_p (out-of-plane, vertical) uncertainty mr
    
    #beam direction uncertainties (need to ask Mark Jones)
    sig_thb = 0.     #beam theta_b (out-of-plane, vertical, YES! IT IS OPPOSITE) uncertainty
    sig_phb = 0.     #beam phi_b (in-plane, horizontal, YES! IT IS OPPOSITE) uncertainty
    
    # relative sigma
    #final electron energy relative uncertainty (dEf/Ef)
    sig_ef = 9.132e-4
    #beam energy relative uncertainty (dEb / Eb)
    sig_E = 7.498e-4   
    
    #Covariance Errors (Errors due to the correlation between different kinematic variables)
    #For now, I consider only uncertainties in th_e, th_p, dEf/Ef, dEb/Eb. Angle units are in [rad] or [rad**2]
    sig_the_thp = 8.43255e-9     #dth_e * dth_p [rad**2]
    sig_the_Ef = -1.48849e-7      #dth_e * dEf/Ef [rad]
    sig_the_Eb = -1.21324e-7      #dth_e * dEb/Eb [rad]
    sig_thp_Ef = -7.01424e-8      #dth_p * dEf/Ef [rad]
    sig_thp_Eb = -6.26743e-8      #dth_p * dEb/Eb [rad]
    sig_Ef_Eb = 6.83851e-7       #dEf/Ef * dEb/Eb [no units]

    # This is where you enter the estimated uncertainties
    # sigma in mr 
    
    ds_dthe = 0.
    ds_dphe = 0.
    ds_def = 0.
    ds_dthp = 0.
    ds_dphp = 0.
    ds_dthb = 0.
    ds_dphb = 0.
    ds_dE = 0.
    # these derivatives are in %/deg or %/MeV ->convert to %/rad
    for d in data:
        p = d.split(':')
        if len(p) >= 2:
            # there are values separated by colons
            if p[0].find('incident energy') >= 0:
                E_inc = float(p[1])
                found += 1
            if p[0].find('scattering angle') >= 0:
                ef = float(p[2])
                found += 1
        f = d.split()
        if len(f) == 0 :
            continue
        if f[0].find('dth_e') >= 0:
            ds_dthe = float(f[1])/dtr
            found += 1
        if f[0].find('dph_e') >= 0:
            ds_dphe = float(f[1])/dtr
            found += 1
        if f[0].find('def') >= 0:
            ds_def = float(f[1])
            found += 1
        if f[0].find('dth_p') >= 0:
            ds_dthp = float(f[1])/dtr
            found += 1
        if f[0].find('dph_p') >= 0:
            ds_dphp = float(f[1])/dtr
            found += 1
        if f[0].find('dth_b') >= 0:
            ds_dthb = float(f[1])/dtr
            found += 1
        if f[0].find('dph_b') >= 0:
            ds_dphb = float(f[1])/dtr
            found += 1
        if f[0].find('dE') >= 0:
            ds_dE = float(f[1])
            found += 1
    if found != 10 :
        print 'did not find all data !'
    # calculate total errors
    sigma_the = ds_dthe*sig_the*1.e-3    #the 1e-3 is to convert relative error sig_the from mr to rad
    sigma_phe = ds_dphe*sig_phe*1.e-3
    
    sigma_thp = ds_dthp*sig_thp*1.e-3
    sigma_php = ds_dphp*sig_php*1.e-3
    
    sigma_thb = ds_dthb*sig_thb*1.e-3
    sigma_phb = ds_dphb*sig_phb*1.e-3
    
    sigma_ef = ds_def*sig_ef*ef
    sigma_dE = ds_dE*sig_E*E_inc

    #Covariance elements [angles are already in radians]
    #Since the covariance for E_inc, Ef are in relative errors, dE_inc/E_inc, dEf/Ef
    #multiply by Ef or E_inc get absolute errors
    sigma_the_thp = 2.*(ds_dthe*ds_dthp)*sig_the_thp
    sigma_the_Ef = 2.*(ds_dthe*ds_def)*sig_the_Ef * ef
    sigma_the_Eb = 2.*(ds_dthe*ds_dE)*sig_the_Eb *E_inc
    sigma_thp_Ef = 2.*(ds_dthp*ds_def)*sig_thp_Ef * ef
    sigma_thp_Eb = 2.*(ds_dthp*ds_dE)*sig_thp_Eb * E_inc
    sigma_Ef_Eb = 2.*(ds_def*ds_dE)*sig_Ef_Eb * ef * E_inc

    sigma_tot = np.sqrt(sigma_the**2 + 
                        sigma_phe**2 + 
                        sigma_thp**2 + 
                        sigma_php**2 + 
                        sigma_thb**2 + 
                        sigma_phb**2 + 
                        sigma_ef**2 + 
                        sigma_dE**2 +
                        #Covariance contribution
                        sigma_the_thp +
                        sigma_the_Ef  +
                        sigma_the_Eb +
                        sigma_thp_Ef +
                        sigma_thp_Eb +
                        sigma_Ef_Eb)

    problem = np.isnan(sigma_tot)

    #array to store individual systematics errors (relative error to Xsec in %)
    sig_arr = np.array([sigma_the, sigma_phe, sigma_thp, sigma_php, sigma_thb, sigma_phb, sigma_ef, sigma_dE, sigma_tot]) 
    
    #array to store correlated systematic errors (% change in cross section)
    sig_corr_arr = np.array([sigma_the_thp, sigma_the_Ef, sigma_the_Eb, sigma_thp_Ef, sigma_thp_Eb, sigma_Ef_Eb])

    #array to store derivatives
    derv_arr = np.array([ds_dthe, ds_dphe, ds_dthp, ds_dphp, ds_dthb, ds_dphb, ds_def, ds_dE])
    
    #Write derivatives to table
    write_table(sig_arr, sig_corr_arr, derv_arr, ftable_bin, ftable_name)


    if (print_all or problem):
        print 'sigma_the  = ', sigma_the
        print 'sigma_phe = ', sigma_phe
        
        print 'sigma_thp = ', sigma_thp
        print 'sigma_php = ', sigma_php
        
        print 'sigma_thb = ', sigma_thb
        print 'sigma_phb = ', sigma_phb
        
        print 'sigma_ef = ' , sigma_ef
        print 'sigma_dE = ', sigma_dE
        
        print 'sigma_tot = ', sigma_tot
    return sigma_tot
# end of it all


#Define a function to write derivatives to a .txt file
def write_table(sig_arr = np.array([]), sig_corr_arr = np.array([]), derv_arr = np.array([]), bin_arr=np.array([]), ftable_name=''):
    
    #open file to be writte in append mode
    fout = open(ftable_name, 'a')
    
    ib = bin_arr[0]
    pm = bin_arr[1]
    thnq = bin_arr[2]

    #Derivatives are in %/rad or %/MeV
    ds_dthe = derv_arr[0] * 1.e-3    #convert %/rad to %/mr (1 mr = 1e-3 rad)
    ds_dphe = derv_arr[1] * 1.e-3
    ds_dthp = derv_arr[2] * 1.e-3
    ds_dphp = derv_arr[3] * 1.e-3
    ds_dthb = derv_arr[4] * 1.e-3
    ds_dphb = derv_arr[5] * 1.e-3 
    ds_def = derv_arr[6]
    ds_dE = derv_arr[7]

    #Get individual systematic errors [in %]
    sig_the = sig_arr[0] 
    sig_phe = sig_arr[1] 
    sig_thp = sig_arr[2] 
    sig_php = sig_arr[3] 
    sig_thb = sig_arr[4] 
    sig_phb = sig_arr[5]  
    sig_ef = sig_arr[6]
    sig_E = sig_arr[7]
    
    #Get Correlated Errors
    sig_the_thp =  sig_corr_arr[0]
    sig_the_Ef  =  sig_corr_arr[1]
    sig_the_Eb  =  sig_corr_arr[2]
    sig_thp_Ef  =  sig_corr_arr[3]
    sig_thp_Eb  =  sig_corr_arr[4]
    sig_Ef_Eb   =  sig_corr_arr[5]

    sig_tot = sig_arr[8]   #total systematic error (added in quadrature)

    fout.write('%i   %f   %f   %.3e   %.3e   %.3e   %.3e   %.3e   %.3e   %.3e   %.3e   %.3e   %.3e   %.3e   %.3e   %.3e   %.3e   %.3e   %.3e  %.3e   %.3e   %.3e   %.3e   %.3e   %.3e   %.3e\n'%(ib, pm, thnq, ds_dthe, ds_dphe, ds_def, ds_dthp, ds_dphp, ds_dthb, ds_dphb, ds_dE, sig_the, sig_phe, sig_ef, sig_thp, sig_php, sig_thb, sig_phb, sig_E, sig_the_thp, sig_the_Ef, sig_the_Eb, sig_thp_Ef, sig_thp_Eb, sig_Ef_Eb, sig_tot))
    fout.close()
