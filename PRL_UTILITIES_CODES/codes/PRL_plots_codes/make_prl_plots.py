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
from matplotlib import *
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, InsetPosition, mark_inset)
from prl_utilities import *

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

def make_prl_plots():
  
    #Read Hall C (E12-10-003) experimental data (requires <=50% statistical uncertainty, units: fm^3, GeV/c)
    pm_avg35, red_dataXsec_avg_masked35, red_dataXsec_avg_tot_err35 = read_hallc_data(35)
    pm_avg45, red_dataXsec_avg_masked45, red_dataXsec_avg_tot_err45 = read_hallc_data(45)
    pm_avg75, red_dataXsec_avg_masked75, red_dataXsec_avg_tot_err75 = read_hallc_data(75)

    #Read Hall A experimental data
    pm_ha35, red_dataXsec_ha35, red_dataXsec_err_ha35 = read_halla_data(35)
    pm_ha45, red_dataXsec_ha45, red_dataXsec_err_ha45 = read_halla_data(45)  
    pm_ha75, red_dataXsec_ha75, red_dataXsec_err_ha75 = read_halla_data(75) 

    #Read Theoretical 
def main():
    print('Entering Main . . .')

    make_prl_plots()

if __name__=="__main__":
    main()


