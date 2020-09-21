import LT.box as B
from LT.datafile import dfile
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy.ma as ma
import sys                                     
import os                                                                                                       
from sys import argv  
import matplotlib
from matplotlib import rc



#Use latex commands (e.g. \textit ot \textbf)
rc('text', usetex=True)
#Set default font to times new roman
font = {'family' : 'Times New Roman',
        'weight' : 'normal',
        'size'   : 12
}
plt.rc('font', **font)

#Set font
csfont = {'fontname':'Times New Roman'}

fname = 'FOFA_ratio_HallC_thnq35_deg.txt'
fname2 = 'FOFA_ratio_HallC_thnq45_deg.txt'


pm35 = np.array(dfile(fname)['pm_avg'])
pm45 = np.array(dfile(fname2)['pm_avg'])   

R_av18_GK_pwia = np.array(dfile(fname)['R_av18_GK_pwia'])   
R_av18_amt_pwia = np.array(dfile(fname)['R_av18_amt_pwia'])   
R_cd_GK_pwia = np.array(dfile(fname)['R_cd_GK_pwia'])   
R_cd_amt_pwia = np.array(dfile(fname)['R_cd_amt_pwia'])   

R_av18_GK2amt_pwia = np.array(dfile(fname)['R_av18_GK2amt_pwia'])   
R_cd_GK2amt_pwia = np.array(dfile(fname)['R_cd_GK2amt_pwia'])   
R_wjc2_GK2amt_pwia = np.array(dfile(fname)['R_wjc2_GK2amt_pwia'])   

R_av18_GK_pwia2 = np.array(dfile(fname2)['R_av18_GK_pwia'])   
R_av18_amt_pwia2 = np.array(dfile(fname2)['R_av18_amt_pwia'])   
R_cd_GK_pwia2 = np.array(dfile(fname2)['R_cd_GK_pwia'])   
R_cd_amt_pwia2 = np.array(dfile(fname2)['R_cd_amt_pwia'])   

R_av18_GK2amt_pwia2 = np.array(dfile(fname2)['R_av18_GK2amt_pwia'])   
R_cd_GK2amt_pwia2 = np.array(dfile(fname2)['R_cd_GK2amt_pwia'])   
R_wjc2_GK2amt_pwia2 = np.array(dfile(fname2)['R_wjc2_GK2amt_pwia'])   

B.pl.figure(1)
B.plot_exp(pm35, R_av18_GK2amt_pwia, marker='o', color='k', label='R (35 deg) = (AV18;GKex05 - AV18;AMT) / AV18;AMT')
B.plot_exp(pm35, R_cd_GK2amt_pwia, marker='s', color='r', label='R (35 deg) = (CD-Bonn;GKex05 - CD-Bonn;AMT) / CD-Bonn;AMT')
B.plot_exp(pm35, R_wjc2_GK2amt_pwia, marker='<', color='m', label='R (35 deg) = (WJC2;GKex05 - WJC2;AMT) / WJC2;AMT')

B.plot_exp(pm45, R_av18_GK2amt_pwia2, marker='^', color='b', label='R (45 deg) = (AV18;GKex05 - AV18;AMT) / AV18;AMT')
B.plot_exp(pm45, R_cd_GK2amt_pwia2, marker='v', color='g', label='R (45 deg) = (CD-Bonn;GKex05 - CD-Bonn;AMT) / CD-Bonn;AMT')
B.plot_exp(pm45, R_wjc2_GK2amt_pwia2, marker='>', color='c', label='R (45 deg) = (WJC2;GKex05 - WJC2;AMT) / WJC2;AMT')

B.pl.ylim(-8, -5.6)

B.pl.title(r'Percent deviation of $\sigma_{red}$ using different form-factor parametrization ')
B.pl.xlabel(r'Neutron Recoil Momenta')
B.pl.ylabel(r'Percent Deviation')
B.pl.legend()
B.pl.show()



'''
B.pl.figure(1)
B.plot_exp(pm35, R_av18_GK_pwia, marker='o', color='k', label='R = (AV18;GKex05 - AV18;JJK) / AV18;JJK')
B.plot_exp(pm35, R_av18_amt_pwia, marker='s', color='r', label='R = (AV18;AMT - AV18;JJK) / AV18;JJK')

B.plot_exp(pm35, R_cd_GK_pwia, marker='^', color='b', label='R = (CD-Bonn;GKex05 - CD-Bonn;JJK) / CD-Bonn;JJK')
B.plot_exp(pm35, R_cd_amt_pwia, marker='v', color='g', label='R =  (CD-Bonn;AMT - CD-Bonn;JJK) / CD-Bonn;JJK ')

B.pl.title(r'Percent deviation of $\sigma_{red}$ using different form-factor parametrization (35 deg) ')
B.pl.xlabel(r'Neutron Recoil Momenta')
B.pl.ylabel(r'Percent Deviation')
B.pl.legend()
B.pl.show()

B.pl.figure(2)
B.plot_exp(pm45, R_av18_GK_pwia2, marker='o', color='k', label='R = (AV18;GKex05 - AV18;JJK) / AV18;JJK')
B.plot_exp(pm45, R_av18_amt_pwia2, marker='s', color='r', label='R = (AV18;AMT - AV18;JJK) / AV18;JJK')

B.plot_exp(pm45, R_cd_GK_pwia2, marker='^', color='b', label='R = (CD-Bonn;GKex05 - CD-Bonn;JJK) / CD-Bonn;JJK')
B.plot_exp(pm45, R_cd_amt_pwia2, marker='v', color='g', label='R =  (CD-Bonn;AMT - CD-Bonn;JJK) / CD-Bonn;JJK ')

B.pl.title(r'Percent deviation of $\sigma_{red}$ using different form-factor parametrization (45 deg) ')
B.pl.xlabel(r'Neutron Recoil Momenta')
B.pl.ylabel(r'Percent Deviation')
B.pl.legend()
B.pl.show()
'''
