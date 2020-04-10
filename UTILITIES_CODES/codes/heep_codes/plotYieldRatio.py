import numpy as np
from LT.datafile import dfile
import LT.box as B

kin = dfile('elasticsYield.txt')

shms_rate = kin['shms_rate'] #in kHz
dataW = np.array(kin['dataW'])
simcW = np.array(kin['simcW'])
dataW_err = np.sqrt(kin['dataW_err'])
simcW_err = np.sqrt(kin['simcW_err'])
R = dataW/simcW
R_err = R * np.sqrt((dataW_err/dataW)**2 + (simcW_err/simcW)**2)

print('R = ', R, ' +/- ', R_err)
B.plot_exp(shms_rate, R, R_err, marker='o', color='r', label='Full Corrections', markersize=4.5)
B.pl.legend()
B.pl.grid(True)
B.pl.xticks(shms_rate, ('44.64', '77.71', '271.98', '590.26'))
B.pl.ylim(0.84,1.01)
B.pl.axhline(1., color='k', linestyle='--')
B.pl.axvline(120., color='b', linestyle='--')
B.pl.axvline(170., color='b', linestyle='--')
B.pl.xlabel(r'SHMS Rate [kHz]')
B.pl.ylabel(r'Y$_{\mathrm{data}}$/Y$_{\mathrm{SIMC}}$')
B.pl.title(r'H$(e,e\'p)$ Elastics Yield Ratio')
B.pl.show()
