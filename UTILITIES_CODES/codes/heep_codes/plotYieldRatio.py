import numpy as np
from LT.datafile import dfile
import LT.box as B
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib import *
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, InsetPosition, mark_inset)


rc('text', usetex=True)
#Set default font to times new roman
font = {'family' : 'Times New Roman',
        'weight' : 'normal',
        'size'   : 12
}
B.pl.rc('font', **font)


kin = dfile('elasticsYield.txt')
report = dfile('../../root_files/pre_HEEP_ELASTICS/report_heep.dat')
etrkEff = report['eTrkEff']
etrkEff_err = report['eTrkEff_err']
shms_rate = report['ptrig1_rate'] #in kHz
dataW = np.array(kin['dataW'])
simcW = np.array(kin['simcW'])
dataW_err = np.sqrt(kin['dataW_err'])
simcW_err = np.sqrt(kin['simcW_err'])
simcW_fofa_err = np.array([0.11178, 0.05785, 0.02722, 0.0029311]) #relative Xsec error dSig/Sig = dY / Y relative yield error
simcW_tot_err = np.sqrt((simcW_err/simcW)**2 + simcW_fofa_err**2)
R = dataW/simcW
R_err = R * np.sqrt((dataW_err/dataW)**2 + (simcW_err/simcW)**2)  #SIMC Monte Carlo Error
R_err_syst = R * np.sqrt((dataW_err/dataW)**2 + simcW_fofa_err**2)  #FOFA errors
R_err_tot = R * np.sqrt((dataW_err/dataW)**2 + simcW_tot_err**2)  #total error

fig = B.pl.figure(figsize=(7,4))
ax1 = B.pl.subplot(1,1,1)

print('>>>>> R = ', R, ' +/- ', R_err, ' +/- ', R_err_syst )
#B.plot_exp(shms_rate, etrkEff, etrkEff_err, marker='o', color='r', label='SHMS tracking efficinecy', markersize=4.5)

B.plot_exp(shms_rate, R, R_err, marker='o', color='r', label='Statistical Error', markersize=4.5)
B.plot_exp(shms_rate, R, R_err_syst, marker='', ecolor='lightgray', label='Systematics Error (Form Factors)', elinewidth=3)


B.pl.legend()
B.pl.grid(True)
B.pl.xticks(shms_rate, ('44.64', '77.71', '271.98', '590.26'))
B.pl.ylim(0.8,1.1)
B.pl.xlim(-10.0,600.)
B.pl.axvline(0., color='k', linestyle='-')

B.pl.axhline(1., color='k', linestyle='--')
B.pl.axvline(120., color='b', linestyle='--')
B.pl.axvline(170., color='b', linestyle='--')
B.pl.xlabel(r'SHMS Rate [kHz]')
B.pl.ylabel(r'Y$_{\mathrm{data}}$/Y$_{\mathrm{SIMC}}$')
B.pl.title(r'H$(e,e^{\prime}p)$ Elastics Yield Ratio')

#Set second axis
ax2 = ax1.twiny()

newlabel = [4.95, 4.03, 2.87, 2.19]  #Q2 values
newpos = np.array([44.64, 77.71, 271.98, 590.26])  #corresponding to old positions

ax2.set_xticks(newpos)
ax2.set_xticklabels(newlabel)
ax2.xaxis.set_ticks_position('bottom') # set the position of the second x-axis to bottom
ax2.xaxis.set_label_position('bottom') # set the position of the second x-axis to bottom
ax2.spines['bottom'].set_position(('outward', 36))
ax2.set_xlabel('$Q^{2}$ [GeV$^{2}$]')
ax2.set_xlim(ax1.get_xlim())

fig.tight_layout()
'''
# Move twinned axis ticks and label from top to bottom
ax2.xaxis.set_ticks_position("bottom")
ax2.xaxis.set_label_position("bottom")

# Offset the twin axis below the host
ax2.spines["bottom"].set_position(("axes", -0.15))

# Turn on the frame for the twin axis, but then hide all 
# but the bottom spine
ax2.set_frame_on(True)
ax2.patch.set_visible(False)
for sp in ax2.spines.itervalues():
    sp.set_visible(False)
ax2.spines["bottom"].set_visible(True)

ax2.set_xticks(new_tick_locations)
ax2.set_xticklabels(['3.2', '4.5', '3.1', '1.1'])
ax2.set_xlabel(r"Modified x-axis: $1/(1+X)$")
#plt.show()
'''
B.pl.show()

'''
Run Number: 3288
('Beam Energy: ', 10.6005, ' GeV')
('electron momentum: ', 8.4407, ' GeV')
('electron angle: ', 12.194, ' deg')
('Q2: ', 4.037495429489856, ' GeV2')
('Xsec = ', 4.357186983589756e-05, ' +/- ', 2.5207087137165505e-06, 'microbarn/sr')
('Xsec_rel_err=', 0.05785174524779779)
Run Number: 3371
('Beam Energy: ', 10.6005, ' GeV')
('electron momentum: ', 7.951, ' GeV')
('electron angle: ', 13.93, ' deg')
('Q2: ', 4.957523420004123, ' GeV2')
('Xsec = ', 8.014773860717714e-06, ' +/- ', 8.958921366947425e-07, 'microbarn/sr')
('Xsec_rel_err=', 0.1117800891533222)
Run Number: 3374
('Beam Energy: ', 10.6005, ' GeV')
('electron momentum: ', 9.0603, ' GeV')
('electron angle: ', 9.928, ' deg')
('Q2: ', 2.876472983977319, ' GeV2')
('Xsec = ', 0.0005118915789305864, ' +/- ', 1.3935587849740395e-05, 'microbarn/sr')
('Xsec_rel_err=', 0.027223709909144824)
Run Number: 3377
('Beam Energy: ', 10.6005, ' GeV')
('electron momentum: ', 9.425, ' GeV')
('electron angle: ', 8.495, ' deg')
('Q2: ', 2.1922671879570754, ' GeV2')
('Xsec = ', 0.002931134650051497, ' +/- ', 6.711949111983509e-05, 'microbarn/sr')
('Xsec_rel_err=', 0.02289880852749493)
'''
