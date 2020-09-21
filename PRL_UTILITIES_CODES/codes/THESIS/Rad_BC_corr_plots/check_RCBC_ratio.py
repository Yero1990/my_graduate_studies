#This code check the ratio of norad / rad for PWIA and FSI

from LT.datafile import dfile
import LT.box as B
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc
import matplotlib.ticker as tic
import sys
import os
import numpy.ma as ma

'''
#Set Global Param
params = {'legend.fontsize': 14,
          'lines.markersize': 9.0,
           'legend.markerscale':1.,
          'figure.figsize': (9, 7),
         'axes.labelsize': 20,
         'axes.labelweight': 'bold',          
         'axes.titlesize':22,
          'axes.titleweight':'bold',
          'axes.linewidth':2,
         'xtick.major.width':2,
          'ytick.major.width':2,          
         'xtick.labelsize':20,
          'ytick.labelsize':20} 
           
plt.rcParams.update(params)
'''
#Use latex commands (e.g. \textit ot \textbf)
rc('text', usetex=True)
#Set default font to times new roman
font = {'family' : 'Times New Roman',
        'weight' : 'bold',
        'size'   : 14
}
plt.rc('font', **font)

def convert2NaN(arr=np.array([]), value=0):
    #method to convert a specified value in a array to nan (not a number)
    
    for i in enumerate(arr):
        if arr[i[0]]==value:
            arr[i[0]] = np.nan
    return arr


#User Input (Usage: python check_RC_ratio.py 80 1 )
pm_set = int(sys.argv[1])
data_set = int(sys.argv[2])

#csfont = {'fontname':'Times New Roman'}

plt.rcParams["font.family"] = "Times New Roman"

#Read data file
if pm_set==80:
    fname = './pm%i_laget_RCratio.txt' % (pm_set)
else:
    fname = './pm%i_laget_RCratio_set%i.txt' %(pm_set, data_set)

f = dfile(fname)

#Get Relevant headers to plot
pm = np.array(f['yb'])    #missing momentum bin
thnq = np.array(f['xb'])  #theta_nq bin



#define function to get RC ratio
def get_RC_ratio(pm_set, data_set, model=''):

    
    #-----Read radiative correction ratio data file----
    if pm_set==80:
        fname = './pm%i_laget_RCratio.txt' % (pm_set)
    else:
        fname = './pm%i_laget_RCratio_set%i.txt' %(pm_set, data_set)

    f = dfile(fname)    

    #Get Relevant headers to plot
    pm = np.array(f['yb'])    #missing momentum bin
    thnq = np.array(f['xb'])  #theta_nq bin 

    if(model=='fsi'):
        
        #Get Ratio and Error
        RC_fact_fsi = np.array(f['fsiRC_ratio'])
        RC_fact_fsi_err = np.array(f['fsiRC_ratio_err'])
        fsiRC_dataXsec = np.array(f['fsiRC_dataXsec'])
        fsiRC_dataXsec_err = np.array(f['fsiRC_dataXsec_err'])
        
        RC_fact_fsi_m = np.ma.array(RC_fact_fsi, mask=((fsiRC_dataXsec==-1.)))
        RC_fact_fsi_m = np.ma.filled(RC_fact_fsi_m.astype(float), np.nan)


        return [RC_fact_fsi_m, RC_fact_fsi_err]
    
    elif(model=='pwia'):
        #Get Ratio and Error
        RC_fact_pwia = np.array(f['pwiaRC_ratio'])
        RC_fact_pwia_err = np.array(f['pwiaRC_ratio_err'])
        pwiaRC_dataXsec = np.array(f['pwiaRC_dataXsec'])
        pwiaRC_dataXsec_err = np.array(f['pwiaRC_dataXsec_err'])

        RC_fact_pwia_m = np.ma.array(RC_fact_pwia, mask=((pwiaRC_dataXsec==-1.)))
        RC_fact_pwia_m = np.ma.filled(RC_fact_pwia_m.astype(float), np.nan)

        return [RC_fact_pwia_m, RC_fact_pwia_err]

    
#define function to get BC ratio
def get_BC_ratio(pm_set, data_set, model=''):

    
    #-----Read radiative correction ratio data file----
    if pm_set==80:
        fname = './pm%i_laget_bc_corr.txt' % (pm_set)
    else:
        fname = './pm%i_laget_bc_corr_set%i.txt' %(pm_set, data_set)

    f = dfile(fname)    

    #Get Relevant headers to plot
    pm = np.array(f['yb'])    #missing momentum bin
    thnq = np.array(f['xb'])  #theta_nq bin 

    if(model=='fsi'):
        
        #Get Ratio and Error
        BC_fact_fsi = np.array(f['bc_fact_fsi'])
        BC_fact_fsi_err = np.array(f['bc_fact_fsi_err'])
        fsiRC_dataXsec = np.array(f['fsiRC_dataXsec'])
        fsiRC_dataXsec_err = np.array(f['fsiRC_dataXsec_err'])
        
        BC_fact_fsi_m = np.ma.array(BC_fact_fsi, mask=((fsiRC_dataXsec==-1.)))
        BC_fact_fsi_m = np.ma.filled(BC_fact_fsi_m.astype(float), np.nan)


        return [BC_fact_fsi_m, BC_fact_fsi_err]
    
    elif(model=='pwia'):
        #Get Ratio and Error
        BC_fact_pwia = np.array(f['bc_fact_pwia'])
        BC_fact_pwia_err = np.array(f['bc_fact_pwia_err'])
        pwiaRC_dataXsec = np.array(f['pwiaRC_dataXsec'])
        pwiaRC_dataXsec_err = np.array(f['pwiaRC_dataXsec_err'])

        BC_fact_pwia_m = np.ma.array(BC_fact_pwia, mask=((pwiaRC_dataXsec==-1.)))
        BC_fact_pwia_m = np.ma.filled(BC_fact_pwia_m.astype(float), np.nan)

        return [BC_fact_pwia_m, BC_fact_pwia_err]
    
    
#Get Corr. Factor for each Pm set

#RC (80 MeV)
RC_fact_fsi_80 = get_RC_ratio(pm_set=80, data_set=1, model='fsi')[0]
RC_fact_fsi_80_err = get_RC_ratio(pm_set=80, data_set=1, model='fsi')[1]

RC_fact_pwia_80 = get_RC_ratio(pm_set=80, data_set=1, model='pwia')[0]
RC_fact_pwia_80_err = get_RC_ratio(pm_set=80, data_set=1, model='pwia')[1]

#BC (80 MeV)
BC_fact_fsi_80 = get_BC_ratio(pm_set=80, data_set=1, model='fsi')[0]
BC_fact_fsi_80_err = get_BC_ratio(pm_set=80, data_set=1, model='fsi')[1]

BC_fact_pwia_80 = get_BC_ratio(pm_set=80, data_set=1, model='pwia')[0]
BC_fact_pwia_80_err = get_BC_ratio(pm_set=80, data_set=1, model='pwia')[1]
#----------------

#RC (580 set1)
RC_fact_fsi_580s1 = get_RC_ratio(pm_set=580, data_set=1, model='fsi')[0]
RC_fact_fsi_580s1_err = get_RC_ratio(pm_set=580, data_set=1, model='fsi')[1]

RC_fact_pwia_580s1 = get_RC_ratio(pm_set=580, data_set=1, model='pwia')[0]
RC_fact_pwia_580s1_err = get_RC_ratio(pm_set=580, data_set=1, model='pwia')[1]

#BC (580 set1)
BC_fact_fsi_580s1 = get_BC_ratio(pm_set=580, data_set=1, model='fsi')[0]
BC_fact_fsi_580s1_err = get_BC_ratio(pm_set=580, data_set=1, model='fsi')[1]

BC_fact_pwia_580s1 = get_BC_ratio(pm_set=580, data_set=1, model='pwia')[0]
BC_fact_pwia_580s1_err = get_BC_ratio(pm_set=580, data_set=1, model='pwia')[1]
#----------------

#RC (580 set2)
RC_fact_fsi_580s2 = get_RC_ratio(pm_set=580, data_set=2, model='fsi')[0]
RC_fact_fsi_580s2_err = get_RC_ratio(pm_set=580, data_set=2, model='fsi')[1]

RC_fact_pwia_580s2 = get_RC_ratio(pm_set=580, data_set=2, model='pwia')[0]
RC_fact_pwia_580s2_err = get_RC_ratio(pm_set=580, data_set=2, model='pwia')[1]

#BC (580 set2)
BC_fact_fsi_580s2 = get_BC_ratio(pm_set=580, data_set=2, model='fsi')[0]
BC_fact_fsi_580s2_err = get_BC_ratio(pm_set=580, data_set=2, model='fsi')[1]

BC_fact_pwia_580s2 = get_BC_ratio(pm_set=580, data_set=2, model='pwia')[0]
BC_fact_pwia_580s2_err = get_BC_ratio(pm_set=580, data_set=2, model='pwia')[1]
#----------------

#RC (750 set1)
RC_fact_fsi_750s1 = get_RC_ratio(pm_set=750, data_set=1, model='fsi')[0]
RC_fact_fsi_750s1_err = get_RC_ratio(pm_set=750, data_set=1, model='fsi')[1]

RC_fact_pwia_750s1 = get_RC_ratio(pm_set=750, data_set=1, model='pwia')[0]
RC_fact_pwia_750s1_err = get_RC_ratio(pm_set=750, data_set=1, model='pwia')[1]

#BC (750 set1)
BC_fact_fsi_750s1 = get_BC_ratio(pm_set=750, data_set=1, model='fsi')[0]
BC_fact_fsi_750s1_err = get_BC_ratio(pm_set=750, data_set=1, model='fsi')[1]

BC_fact_pwia_750s1 = get_BC_ratio(pm_set=750, data_set=1, model='pwia')[0]
BC_fact_pwia_750s1_err = get_BC_ratio(pm_set=750, data_set=1, model='pwia')[1]
#----------------

#RC (750 set2)
RC_fact_fsi_750s2 = get_RC_ratio(pm_set=750, data_set=2, model='fsi')[0]
RC_fact_fsi_750s2_err = get_RC_ratio(pm_set=750, data_set=2, model='fsi')[1]

RC_fact_pwia_750s2 = get_RC_ratio(pm_set=750, data_set=2, model='pwia')[0]
RC_fact_pwia_750s2_err = get_RC_ratio(pm_set=750, data_set=2, model='pwia')[1]

#BC (750 set2)
BC_fact_fsi_750s2 = get_BC_ratio(pm_set=750, data_set=2, model='fsi')[0]
BC_fact_fsi_750s2_err = get_BC_ratio(pm_set=750, data_set=2, model='fsi')[1]

BC_fact_pwia_750s2 = get_BC_ratio(pm_set=750, data_set=2, model='pwia')[0]
BC_fact_pwia_750s2_err = get_BC_ratio(pm_set=750, data_set=2, model='pwia')[1]
#----------------

#RC (750 set3)
RC_fact_fsi_750s3 = get_RC_ratio(pm_set=750, data_set=3, model='fsi')[0]
RC_fact_fsi_750s3_err = get_RC_ratio(pm_set=750, data_set=3, model='fsi')[1]

RC_fact_pwia_750s3 = get_RC_ratio(pm_set=750, data_set=3, model='pwia')[0]
RC_fact_pwia_750s3_err = get_RC_ratio(pm_set=750, data_set=3, model='pwia')[1]

#BC (750 set3)
BC_fact_fsi_750s3 = get_BC_ratio(pm_set=750, data_set=3, model='fsi')[0]
BC_fact_fsi_750s3_err = get_BC_ratio(pm_set=750, data_set=3, model='fsi')[1]

BC_fact_pwia_750s3 = get_BC_ratio(pm_set=750, data_set=3, model='pwia')[0]
BC_fact_pwia_750s3_err = get_BC_ratio(pm_set=750, data_set=3, model='pwia')[1]
#----------------

#label_size = 15
#plt.rcParams['xtick.labelsize'] = label_size
#plt.rcParams['ytick.labelsize'] = label_size



thnq_arr = [5, 15, 25, 35, 45, 55, 65, 75, 85, 95, 105]
#thnq_arr = [55]

for i, ithnq in enumerate(thnq_arr):
    
    print('thnq=',ithnq)
         
    B.pl.clf()
    B.pl.figure(i)

    #Plot RC ratio PWIA and FSI
    B.plot_exp(pm[thnq==ithnq], RC_fact_pwia_80[thnq==ithnq], RC_fact_pwia_80_err[thnq==ithnq], marker='o', color='k', label='PWIA')
    B.plot_exp(pm[thnq==ithnq], RC_fact_pwia_580s1[thnq==ithnq], RC_fact_pwia_580s1_err[thnq==ithnq], marker='s', color='b', logy=False, label='PWIA')
    B.plot_exp(pm[thnq==ithnq], RC_fact_pwia_580s2[thnq==ithnq], RC_fact_pwia_580s2_err[thnq==ithnq], marker='^', color='r', logy=False, label='PWIA')
    B.plot_exp(pm[thnq==ithnq], RC_fact_pwia_750s1[thnq==ithnq], RC_fact_pwia_750s1_err[thnq==ithnq], marker='v', color='g', logy=False, label='PWIA')
    B.plot_exp(pm[thnq==ithnq], RC_fact_pwia_750s2[thnq==ithnq], RC_fact_pwia_750s2_err[thnq==ithnq], marker='p', color='m', logy=False, label='PWIA')
    B.plot_exp(pm[thnq==ithnq], RC_fact_pwia_750s3[thnq==ithnq], RC_fact_pwia_750s3_err[thnq==ithnq], marker='D', color='c', logy=False, label='PWIA')

     
    B.plot_exp(pm[thnq==ithnq], RC_fact_fsi_80[thnq==ithnq],  RC_fact_fsi_80_err[thnq==ithnq], marker='o', markerfacecolor='w', color='k', label='FSI 80 (set1)')
    B.plot_exp(pm[thnq==ithnq], RC_fact_fsi_580s1[thnq==ithnq],  RC_fact_fsi_580s1_err[thnq==ithnq], marker='s',  markerfacecolor='w',color='b', label='FSI 580 (set1)')
    B.plot_exp(pm[thnq==ithnq], RC_fact_fsi_580s2[thnq==ithnq],  RC_fact_fsi_580s2_err[thnq==ithnq], marker='^',  markerfacecolor='w',color='r', label='FSI 580 (set2)')
    B.plot_exp(pm[thnq==ithnq], RC_fact_fsi_750s1[thnq==ithnq],  RC_fact_fsi_750s1_err[thnq==ithnq], marker='v',  markerfacecolor='w',color='g', label='FSI 750 (set1)')
    B.plot_exp(pm[thnq==ithnq], RC_fact_fsi_750s2[thnq==ithnq],  RC_fact_fsi_750s2_err[thnq==ithnq], marker='p',  markerfacecolor='w',color='m', label='FSI 750 (set2)')
    B.plot_exp(pm[thnq==ithnq], RC_fact_fsi_750s3[thnq==ithnq],  RC_fact_fsi_750s3_err[thnq==ithnq], marker='D',  markerfacecolor='w',color='c', label='FSI 750 (set3)')

    B.pl.legend(ncol=2, loc='upper right', prop={'size': 10}) 

    B.pl.ylim(0.0, 5.)
    B.pl.title(r'Radiative Correction Factor, $\theta_{\mathrm{nq}}=%i \pm 5^{\circ}$'%(ithnq), fontsize=15)                          #Set common title
    B.pl.xlabel(r'$p_{\mathrm{r}}$ [GeV/c]', fontsize=14)                           #Set common x-label
    B.pl.ylabel(r'f$_{\mathrm{rad}}$=$\sigma_{\mathrm{norad}}/\sigma_{\mathrm{rad}}$', fontsize=14)  #Set common y-label
    B.pl.xticks(fontsize=14)
    B.pl.yticks(fontsize=14)
    B.pl.tight_layout()
    #B.pl.show()

    B.pl.savefig('./RCratio_thnq%i.pdf'%(ithnq))

    
    B.pl.clf()
    B.pl.figure(i)

    #Plot BC ratio PWIA and FSI
    B.pl.hlines(1.1, 0., 1.2, linestyles='dashed', label=r'$\pm 10\%$' )
    B.pl.hlines(0.9,0., 1.2, linestyles='dashed')
    B.plot_exp(pm[thnq==ithnq], BC_fact_pwia_80[thnq==ithnq], BC_fact_pwia_80_err[thnq==ithnq], marker='o', color='k', label='PWIA')
    B.plot_exp(pm[thnq==ithnq], BC_fact_pwia_580s1[thnq==ithnq], BC_fact_pwia_580s1_err[thnq==ithnq], marker='s', color='b', logy=False, label='PWIA')
    B.plot_exp(pm[thnq==ithnq], BC_fact_pwia_580s2[thnq==ithnq], BC_fact_pwia_580s2_err[thnq==ithnq], marker='^', color='r', logy=False, label='PWIA')
    B.plot_exp(pm[thnq==ithnq], BC_fact_pwia_750s1[thnq==ithnq], BC_fact_pwia_750s1_err[thnq==ithnq], marker='v', color='g', logy=False, label='PWIA')
    B.plot_exp(pm[thnq==ithnq], BC_fact_pwia_750s2[thnq==ithnq], BC_fact_pwia_750s2_err[thnq==ithnq], marker='p', color='m', logy=False, label='PWIA')
    B.plot_exp(pm[thnq==ithnq], BC_fact_pwia_750s3[thnq==ithnq], BC_fact_pwia_750s3_err[thnq==ithnq], marker='D', color='c', logy=False, label='PWIA')

    B.pl.hlines(1.2, 0., 1.2, linestyles='dashdot', color='r', label=r'$\pm 20\%$')
    B.pl.hlines(0.8,0., 1.2, linestyles='dashdot', color='r')
    B.plot_exp(pm[thnq==ithnq], BC_fact_fsi_80[thnq==ithnq],  BC_fact_fsi_80_err[thnq==ithnq], marker='o', markerfacecolor='w', color='k', label='FSI, 80 (set1)')
    B.plot_exp(pm[thnq==ithnq], BC_fact_fsi_580s1[thnq==ithnq],  BC_fact_fsi_580s1_err[thnq==ithnq], marker='s',  markerfacecolor='w',color='b', label='FSI, 580 (set1)')
    B.plot_exp(pm[thnq==ithnq], BC_fact_fsi_580s2[thnq==ithnq],  BC_fact_fsi_580s2_err[thnq==ithnq], marker='^',  markerfacecolor='w',color='r', label='FSI, 580 (set2)')
    B.plot_exp(pm[thnq==ithnq], BC_fact_fsi_750s1[thnq==ithnq],  BC_fact_fsi_750s1_err[thnq==ithnq], marker='v',  markerfacecolor='w',color='g', label='FSI, 750 (set1)')
    B.plot_exp(pm[thnq==ithnq], BC_fact_fsi_750s2[thnq==ithnq],  BC_fact_fsi_750s2_err[thnq==ithnq], marker='p',  markerfacecolor='w',color='m', label='FSI, 750 (set2)')
    B.plot_exp(pm[thnq==ithnq], BC_fact_fsi_750s3[thnq==ithnq],  BC_fact_fsi_750s3_err[thnq==ithnq], marker='D',  markerfacecolor='w',color='c', label='FSI, 750 (set3)')


    #Plot horizontal lines to denote +/- 10 and 20 %



    
    B.pl.legend(ncol=2, loc='upper right', prop={'size': 10}) 

    B.pl.ylim(0.5, 1.5)
    B.pl.title(r'Bin Centering Correction Factor, $\theta_{\mathrm{nq}}=%i \pm 5^{\circ}$'%(ithnq), fontsize=15)                          #Set common title
    B.pl.xlabel(r'$p_{\mathrm{r}}$ [GeV/c]', fontsize=14)                           #Set common x-label
    B.pl.ylabel(r'f$_{\mathrm{bc}}=\sigma^{\mathrm{model}}/\bar{\sigma}^{\mathrm{model}}$', fontsize=14)  #Set common y-label
    B.pl.xticks(fontsize=14)
    B.pl.yticks(fontsize=14)
    B.pl.tight_layout()
    B.pl.show()

    B.pl.savefig('./BCratio_thnq%i.pdf'%(ithnq))


    
