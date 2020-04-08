import numpy as np
from LT.datafile import dfile
import LT.box as B

kin = dfile('heep_yield.dat')


run = kin['Run']
shms_rate = kin['shms_3of4_rate']   #kHz
R = kin['R']
R_err = kin['R_err']

B.plot_exp(shms_rate, R, R_err, marker='o', color='red', label='Yield Ratio')
B.pl.xlabel('SHMS 3/4 Rate [kHz]')
B.pl.ylabel(r'Y$_{Data}$/Y$_{SIMC}$')
B.pl.title('H(e,e\'p) Elastic Yield Rtaio')
B.pl.show()
