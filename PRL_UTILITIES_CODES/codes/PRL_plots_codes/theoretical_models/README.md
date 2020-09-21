*********************
 theoretical_models
*********************

In this directory are the cross sections from other theretical models other 
than Laget theory.

The cross-sections were determined in a similar manner as the Laget theory,
by using code external to SIMC and taking the average kinematics as input
parameters.

The models are:

Argonne V18 Potential: Calculation by Misak Sargsian
CD-Bonn Potential: Calculation by Misak Sargsian
Paris Potential: Calculation by J.M. Laget


The Units of Cross Sections are as follows:

K*sigcc1:  [ub * MeV^-1 * sr^-2]  #I calculate based on avg. kin.

Laget theory Xsec:  
From the Laget gfortran code, the Xsec are returned in [fm^2 * sr^-2 * MeV^-1] 
In dir ../theory_Xsec/calc_theory_Xsec.py, I convert to: ub * MeV^-1 * sr^-2, which is the same as SIMC

For V18 and CD-Bonn:
Werner does the calculation, but Im not sure in which units it is.

-----------------------

STEPS:

Use calc_redXsec.py code
1) Convert Xsec from [nb * GeV^-1 * sr^-2] to [ub * MeV^-1 * sr^-2]
2) Divide the Xsec by Ksig cc1.  Make sure the [xbin,ybin] from Xsec matches Ksig_cc1

Use calculate_averages.py
3) Use average.py (mehtod get_average_array(), to combine the reduced Xsec data)

**NOTE**
The *fofaPBosted* directories contain the reduced cross sections where Ksig_cc1 uses P. Bosted
Form Factor Parametrizations

Otherwise, the (updated_FSI_models, updated_PWIA_models) directories used the JRA parametrization

