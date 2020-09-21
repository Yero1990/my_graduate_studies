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


#1 fm^-1 = 0.1973 GeV

GeV2fm = 0.197**3    #convert GeV^-3 to fm^3

kin = dfile('new_halla/mom_dist_cdbonn_mevc.data')

pm = kin['pm']   #pmiss is in MeV/c
rho = kin['rho'] * GeV2fm  #momentum distributions in 1 / GeV^3

