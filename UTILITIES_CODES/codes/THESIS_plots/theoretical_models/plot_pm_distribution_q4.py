# compare the experimental red. cross sections as 
# a function of pm with the calculated ones 
#
# 
# plot for q1
#
# The script assumes that all calculations have been
# based on the kinematics in the exp. file
#


from LT.datafile import dfile
import numpy as np
import plotting2 as pl

import match_arrays as MA
import match_string_arrays as MSA

import pdb
from itertools import cycle, islice

hbarc = 197.33
Mp = 0.93827

# momentum distributions
"""
d_paris = dfile('../exp_results/nk_paris.data')
d_cdbonn = dfile('../exp_results/nk_bonn.data')
d_v18 = dfile('../exp_results/v18_mom_dist_short.data')
"""

fm_to_mev = 197.33

# control dashes
dashes0 = (10,3)  # 10 points on, 1 point off, 2 point on, 1 point off etc.
dashes1 = (5,3)  # 5 points on, 1 point off
dashes2 = (5,1)  # 5 points on, 1 point off
dashdot1 = (5,1,1,1)  # 5 points on, 2 point off, 2 point on
dashdot2 = (20,5,5,5)  # 
dots1 = (1,1)  #
dots2 = (2,1)  #
dots3 = (5,2)  #
# solid = (None,None)
solid = []

exp_color = 'k'
CD_color = (0.0,0.,1.)
V18_color= (1., 0., 1.)
paris_color = (0.,.8, 0.)


line_w = 1.75


#---------------------------------------------------------------------- 
# colors
my_colors = ['b','g','r','c','m','y','k']

my_markers_dict = {'circle':'o',\
              'tri_down':'v',\
              'tri_left':'<', \
              'tri_right':'>', \
              'tri_up':'^', \
              'square':'s', \
              'star':'*', \
              'pentagon':'p', \
              'hex1':'h', \
              'hex2':'H', \
              'diamond':'D', \
              'diamond1':'d'}

my_markers =  ['o',\
              'v',\
              '<', \
              '>', \
              '^', \
              's', \
              '*', \
              'p', \
              'h', \
              'H', \
              'D', \
              'd']
# markers for the various theta settings
theta_markers = {'35':'o', \
               '45':'^', \
               '75':'s', \
               '85':'*'}
"""
my_color = ['b','g','r','y','m','c','k','b','g','r','y']
color_base = ['r', 'g', 'b', 'm', 'c', 'k']
"""
my_color = ['b','g','r','y','m','c','k','b','g','r','y']

color_base = ['b', 'g']
marker_base = ['o','v']

def take(iterable, start=0, stop=1):
    "Return first n items of the iterable as a list"
    return list(islice(iterable, start, stop))

color_cycle = cycle(color_base)
marker_cycle = cycle(marker_base)

#---------------------------------------------------------------------- 

#---------------------------------------------------------------------- 
def get_bin_id(d):
    kins = d['kin']
    ix = np.array(d['ix']).astype(int)
    iy = np.array(d['iy']).astype(int)
    id_s = []
    for i, k in enumerate(kins): 
        id_s.append( '{:d}_{:d}_{:s}'.format( ix[i], iy[i], k.strip()))
    return np.array(id_s)
#---------------------------------------------------------------------- 
# match arrays and return numpy arrays
def match(a1,a2):
    m1,m2 = MSA.match_arrays(a1,a2)
    return np.array(m1), np.array(m2) 
    
# necessary tools to calculate averages and remove duplicates
#----------------------------------------------------------------------
def remove_duplicates(array):
    na = []
    for a in array:
        try:
            na.index(a)
        except:
            na.append(a)
    return na
        
#----------------------------------------------------------------------
def find_duplicates(np_array):
    # np_array is a numpy array
    # return a list with np_array values that contain no duplicates and a
    # corresponding list of indices pointing to identical values in np_array
    array_no_duplicates = remove_duplicates( list(np_array) )
    # sort the array 
    array_no_duplicates.sort()
    index_list = []
    for a in array_no_duplicates:
        in_dup = np.where(np_array == a)[0]
        index_list.append( list(in_dup) )
    return array_no_duplicates, index_list
#----------------------------------------------------------------------
#----------------------------------------------------------------------
def calc_p_beam_kin(f):
    # calculate protn beam which correponds to final state
    d = dfile(f)
    x = np.array(d['x'])
    pr = np.array(d['pr'])
    pm_b = np.array(d['pm_b'])
    q2 = np.array(d['q2']) 
    Ep = q2/(2.*Mp*x) + Mp
    p = np.sqrt(Ep**2 - Mp**2)
    s = (2.*Mp)**2 + q2*(2./x -1.)
    sqs = np.sqrt(s)
    return x, pr, pm_b, q2, p, sqs
#----------------------------------------------------------------------
def calc_inv_mass(f):
    # calculate inveriant mass off the proton for given kinematics
    d = dfile(f)
    x = np.array(d['x'])
    pr = np.array(d['pr'])
    pm_b = np.array(d['pm_b'])
    q2 = np.array(d['q2'])
    W2 = Mp**2 + q2*(1.-x)/x
    W = np.sqrt(W2)
    return x, pr, pm_b, q2, W

#----------------------------------------------------------------------
def get_averages(pm, sig, dsig):
    # calculate the weighted mean of sig and pm where the weight
    # is given by 1./dsig**2
    weight = 1./dsig**2
    sig_av = np.sum(sig*weight)/np.sum(weight)
    pm_av = np.sum(pm*weight)/np.sum(weight)
    dsig_av_sq = 1./np.sum(weight)
    dsig_av = np.sqrt( dsig_av_sq)
    # include the chi2 of the weighted mean in the error
    if (len(pm) > 1):
        chi2 = (1./(len(pm) - 1.)) * np.sum( (sig-sig_av)**2*weight )
        dsig1_av = np.sqrt( dsig_av_sq * chi2 )
    else:
        dsig1_av = dsig_av
    return pm_av, sig_av, dsig_av, dsig1_av
#----------------------------------------------------------------------
def get_array_average(pm, sig, dsig):
    # calculate the averages for the given arrays, where
    pm_av, pm_ind = find_duplicates(pm)
    pm_av_c = np.array([])
    sig_av = np.array([])
    dsig_av = np.array([])
    dsig1_av = np.array([])
    for p in pm_ind:
        pma, siga, dsiga,  dsiga1 = get_averages (pm[p], sig[p], dsig[p] )
        pm_av_c = np.append(pm_av_c, pma)
        sig_av = np.append(sig_av, siga)
        dsig_av = np.append(dsig_av, dsiga)
        # error including chisq
        dsig1_av = np.append(dsig1_av, dsiga1)
    return pm_av_c, sig_av, dsig_av, dsig1_av
#----------------------------------------------------------------------#----------------------------------------------------------------------
def plot_theory(pm_exp, rho, **kwargs):
    # calculate simple average
    # now calculate the averages
    pm, pm_ind = find_duplicates(pm_exp)
    rho_av = np.array([])
    for p in pm_ind:
        # weight average with cross section
        # pma, rhoav, dr,  dr1 = get_averages (pm_exp[p], rho[p], np.sqrt(rho[p]) )
        # no weight for the average
        pma, rhoav, dr,  dr1 = get_averages (pm_exp[p], rho[p], np.ones_like(rho[p]) )
        rho_av = np.append(rho_av, rhoav)
    # pdb.set_trace()
    return pl.plot_line(pm, rho_av, **kwargs)
#----------------------------------------------------------------------

def plot_all(theta, scale = 1., sigexp_scale = 1.e7, **kwargs):
    # unit conversions for rho
    rho_scale = 1.e-3
    #deutpwia calc
    try:
        dp = dfile(pwia_dir + pwia_file(theta))
    except:
        print "file not found : ", pwia_dir + pwia_file(theta)
        return(None, 0, 0)
    #
    sigpwia = np.array(dp.get_data('sigpwia'))
    sigep = np.array(dp.get_data('sig_cc1'))
    rho = np.array(dp.get_data('rho'))
    k_fact = np.array(dp.get_data('k_fact'))
    f_rec = np.array(dp.get_data('f_rec'))
    ppi = np.array(dp.get_data('p_i'))
    R = 1./(sigep*k_fact*f_rec)
    rho_p = sigpwia*R*hbarc**3
    #
    # bin_corrections     
    dbc = dfile(bc_dir + bc_file(theta))
    # calculated bin corrections    
    sig_laget_FSI = np.array(dbc.get_data('sig2'))
    bc = np.array(dbc.get_data('bc'))
    # 
    if not use_bin_corr:
        bc = np.ones_like(sig_laget_FSI)
    # sig_red_laget_FSI_avkin = sig_red_laget_FSI/bc
    # experimental values
    dexp = dfile(exp_dir + exp_file(theta))

    # find the good events: no left and good cross section & apply bin corrections
    sig_exp_av = np.array( dexp.get_data('sig_exp') )
    # is_ok = np.where(sig_exp_av > 0.)[0]
    is_ok = (sig_exp_av > 0.)
    # exp. bin matching
    ixy_exp = get_bin_id(dexp)[is_ok]    
    # create the experimental arrays and apply bin corrections
    pm_exp = np.array( dexp.get_data('pm_av') )[is_ok]/1000.
    sig_red_exp = np.array(dexp.get_data('sig_exp') )[is_ok]*bc*sigexp_scale
    dsig_red_exp = np.array( dexp.get_data('dsig_exp') )[is_ok]*bc*sigexp_scale
    exp_rel_err = dsig_red_exp/sig_red_exp
    # select good statistics (smaller that 50% rel error)
    stat_ok = (dsig_red_exp/sig_red_exp < 0.4)
    R_ok = R[stat_ok]
    pm_exp_ok = pm_exp[stat_ok]
    # get exp. cross sections
    sig_exp = (np.array(dexp.get_data('csec_exp') )[is_ok]*bc/rho_scale)[stat_ok]
    sig_pwia = (np.array(dexp.get_data('csec_pwia') )[is_ok]*bc/rho_scale)[stat_ok]
    dsig_exp = sig_exp * exp_rel_err[stat_ok]
    rho_exp = sig_exp * R_ok*hbarc**3 
    drho_exp = dsig_exp * R_ok*hbarc**3
    
    
    # calculate reduced cross section like theory
    
    # PARIS calculation --------------------------------------------------
    # calculate averages
    pm_exp_av, pm_ind = find_duplicates(pm_exp)
    # pm_av contains the averaged pm values
    # pm_ind is an array of indices for the identical pm values to be averaged
    # now calculate the averages
    # pm_av_r, sig_red_av, dsig_red_av, dsig1_red_av = get_array_average(pm_exp, sig_red_exp, dsig_red_exp)
    pm_av_r, sig_red_av, dsig_red_av, dsig1_red_av = get_array_average(pm_exp_ok, rho_exp, drho_exp)
    #
    # for legend handling
    legend_txt = []
    graphs = []
    #  PWIA Laget paris.
    if (l_dir_pp != ""):
        try:
            dl = dfile(l_dir_pp + l_file_pp(theta))
            # pm calculated by the code
            # pm bin numbers for matching
            ixy_m = get_bin_id(dl)
            # find matching pm_b for ratio
            m1,m2 = match(ixy_m, ixy_exp)
            # averaged pm for the bins
            pm_m = np.array( dl.get_data('pm_b_av'))
            sig_mp = np.array( dl.get_data('sig2'))[m1]*1.e7 # 
            # rho in fm**3
            rho_m = sig_mp*R[m2]*hbarc**3
            ppl0, = plot_theory(pm_exp[m2], rho_m * scale, color = paris_color , lw = line_w)
            ppl0.set_dashes( dots3 )
            graphs.append(ppl0)
            legend_txt.append('Paris PWIA')
        except:
            print 'Problem plotting ', l_file_pp(theta), ' !'
    
    #  FSI Laget paris.
    if (l_dir_pf != ""):
        try:
            dl = dfile(l_dir_pf + l_file_pf(theta))
            # pm calculated by the code
            # pm bin numbers for matching
            ixy_m = get_bin_id(dl)
            # find matching pm_b for ratio
            m1,m2 = match(ixy_m, ixy_exp)
            # averaged pm for the bins
            pm_m = np.array( dl.get_data('pm_b_av'))
            sig_m = np.array( dl.get_data('sig2'))[m1]*1.e7 # 
            # rho in fm**3
            rho_m = sig_m*R[m2]*hbarc**3
            pfl0, = plot_theory(pm_exp[m2], rho_m * scale, color = paris_color, lw = line_w)
            pfl0.set_dashes( solid )
            graphs.append(pfl0)
            legend_txt.append('Paris FSI')
        except:
            print 'Problem plotting ', l_file_pf(theta), ' !'
    
    #  PWIA misak paris.
    if (m_dir_pp != ""):
        try:
            dm = dfile(m_dir_pp + m_file_pp(theta))
            # pm calculated by the code
            # pm bin numbers for matching
            ixy_m = get_bin_id(dm)
            # find matching pm_b for ratio
            m1,m2 = match(ixy_m, ixy_exp)
            # averaged pm for the bins
            pm_m = np.array( dm.get_data('pm_b_av'))
            sig_m = np.array( dm.get_data('crs12'))[m1]*1.e-3 # 
            # rho in fm**3
            rho_m = sig_m*R[m2]*hbarc**3
            pp0, = plot_theory(pm_exp[m2], rho_m * scale, color = paris_color, lw = line_w)
            pp0.set_dashes( dots3 )
            graphs.append(pp0)
            legend_txt.append('Paris PWIA')
        except:
            print 'Problem plotting ', m_file_pp(theta), ' !'
    if (m_dir_pf != ""):
    #   
    # PARIS FSI calculation
        try:    
            dm = dfile(m_dir_pf + m_file_pf(theta))
            # pm calculated by the code
            # pm bin numbers for matching
            ixy_m = get_bin_id(dm)
            # find matching pm_b for ratio
            m1,m2 = match(ixy_m, ixy_exp)
            # averaged pm for the bins
            pm_m = np.array( dm.get_data('pm_b_av'))
            # apply bin correction
            sig_m = np.array( dm.get_data('crs12'))[m1]*1.e-3 #  * bc
            # rho in fm**3
            rho_m = sig_m*R[m2]*hbarc**3
            pf0, = plot_theory(pm_exp[m2], rho_m * scale, color = paris_color, ls = '-', lw = line_w)
            pf0.set_dashes( solid ) 
            graphs.append(pf0)  
            legend_txt.append('Paris FSI')
        except:
            print 'Problem plotting ', m_file_pf(theta), ' !'
    if (m_dir_cp != ""):
    # CD Bonn calculation --------------------------------------------------
    # PWIA misak CD Bonn
        try:    
            dm = dfile(m_dir_cp + m_file_cp(theta))
            # pm calculated by the code
            # pm bin numbers for matching
            ixy_m = get_bin_id(dm)
            # find matching pm_b for ratio
            m1,m2 = match(ixy_m, ixy_exp)
            # averaged pm for the bins
            pm_m = np.array( dm.get_data('pm_b_av'))
            sig_m = np.array( dm.get_data('crs12'))[m1]*1.e-3 # * bc
            # rho in fm**3
            rho_m = sig_m*R[m2]*hbarc**3
            pp1, = plot_theory(pm_exp[m2], rho_m * scale, color = CD_color, lw = line_w)
            pp1.set_dashes( dots3)
            graphs.append(pp1)
            legend_txt.append('CD-Bonn PWIA')            
        except:
            print 'Problem plotting ', m_file_cp(theta), ' !'        
    #
    # FSI misak CD Bonn
    if (m_dir_cf != ""):
        try:
            dm = dfile(m_dir_cf + m_file_cf(theta))
            # pm calculated by the code
            # pm_m = np.array( dm.get_data('pr'))*1000.
            # pm bin numbers for matching
            ixy_m = get_bin_id(dm)
            # find matching pm_b for ratio
            m1,m2 = match(ixy_m, ixy_exp)
            # averaged pm for the bins
            pm_m = np.array( dm.get_data('pm_b_av'))
            sig_m = np.array( dm.get_data('crs12'))[m1]*1.e-3 #  * bc
            # rho in fm**3
            rho_m = sig_m*R[m2]*hbarc**3
            pf1, = plot_theory(pm_exp[m2], rho_m * scale, color = CD_color, lw = line_w)
            pf1.set_dashes( solid )
            graphs.append(pf1)
            legend_txt.append('CD-Bonn FSI')            
        except:
            print 'Problem plotting ', m_file_cf(theta), ' !'
    # V18 calculation --------------------------------------------------
    # PWIA misak V18
    if (m_dir_vp != ""):
        #try:    
            dm = dfile(m_dir_vp + m_file_vp(theta))
            # pm calculated by the code
            # pm bin numbers for matching
            ixy_m = get_bin_id(dm)
            # find matching pm_b for ratio
            m1,m2 = match(ixy_m, ixy_exp)
            # averaged pm for the bins
            pm_m = np.array( dm.get_data('pm_b_av'))
            sig_m = np.array( dm.get_data('crs12'))[m1]*1.e-3 # * bc
            # rho in fm**3
            rho_m = sig_m*R[m2]*hbarc**3
            pp2, = plot_theory(pm_exp[m2], rho_m * scale, color = V18_color, lw = line_w)
            pp2.set_dashes( dots3)
            graphs.append(pp2)
            legend_txt.append('V18 PWIA')            
        #except:
        #    print 'Problem plotting ', m_file_vp(theta), ' !'        
    #
    # FSI misak V18
    if (m_dir_vf != ""):
        try:
            dm = dfile(m_dir_vf + m_file_vf(theta))
            # pm calculated by the code
            # pm_m = np.array( dm.get_data('pr'))*1000.
            # pm bin numbers for matching
            ixy_m = get_bin_id(dm)
            # find matching pm_b for ratio
            m1,m2 = match(ixy_m, ixy_exp)
            # averaged pm for the bins
            pm_m = np.array( dm.get_data('pm_b_av'))
            sig_m = np.array( dm.get_data('crs12'))[m1]*1.e-3 #  * bc
            # rho in fm**3
            rho_m = sig_m*R[m2]*hbarc**3
            pf2, = plot_theory(pm_exp[m2], rho_m * scale, color = V18_color, lw = line_w)
            pf2.set_dashes( solid )
            graphs.append(pf2)
            legend_txt.append('V18 FSI')            
        except:
            print 'Problem plotting ', m_file_vf(theta), ' !'
    # EXPERIMENT --------------------------------------------------
    # plot experimental_data
    exp_marker = 'o'
    if (pm_av_r.shape[0] == 0):
        return (None, 0, 0)    # skip last exp. point (with huge errorbars)
    if skip_last:
        pm_av_pl = pm_av_r[:-1]
        sig_red_pl = sig_red_av[:-1]
        dsig_red_pl = dsig_red_av[:-1]
    else:
        pm_av_pl = pm_av_r
        sig_red_pl = sig_red_av
        dsig_red_pl = dsig_red_av
    #  
    # first data point
    x0 = pm_av_pl[-1]
    y0 = (sig_red_pl * scale)[-1]
    #      
    pexp = pl.plot_exp(pm_av_pl, sig_red_pl * scale, dy = dsig_red_pl * scale, \
                logy = True, \
                # color = 'k', \
                # marker = exp_marker, \
                capsize = 0., \
                markersize = 4., \
                # below is for empty circles
                # mfc = 'w', \
                # mew = 1.2, \
                **kwargs)
    if skip_legend:
        return (pexp, x0, y0) 
    # make legends
    pl.plt.legend( graphs, legend_txt, loc = 'upper right', frameon = False)
    #pl.plt.legend()
    return (pexp, x0, y0)
#----------------------------------------------------------------------

# plot a panel with 3 angles 35, 45, and 75 degrees

# file names and directories
# ps cut off value to be used
# ps_cut_off = 0.200
# ps_cut_off = 0.100
ps_cut_off = 0.05

pss = '{:.3f}'.format(ps_cut_off)
pss = ''

# deutpwia calc
pwia_dir = "./deutpwia/q4_sig_av_TP/"
# pwia_file = lambda theta: "pwia_calc_theta_%0.2f_"%(theta) + pss + "_oo.data"
pwia_file = lambda theta: "pwia_calc_theta_%0.2f"%(theta) + ".data"

# bin_corrections    
bc_dir = "./laget/q4_sig_avkin_fsi_thnq_pm/"
bc_file = lambda theta: "sig_theta_%0.2f"%(theta) + ".data"  

# experimental values
exp_dir = '../thnq_pm/sig_exp_kinav/'
exp_file = lambda theta: "theta_%0.2f"%(theta) + ".data"

# Laget calculations
#  PWIA  paris.
# set anything non-blank and Laget ratios will be drawn
l_dir_pp = "./laget/q4_sig_avkin_pwia_thnq_pm/"
#l_dir_pp = ""
l_file_pp = lambda theta: "all_sig_theta_%0.2f.data"%(theta)

# FSI calculation
l_dir_pf = "./laget/q4_sig_avkin_fsi_thnq_pm/"
#l_dir_pf = ""
l_file_pf = lambda theta: "all_sig_theta_%0.2f.data"%(theta)


#  PWIA misak paris.
#m_dir_pp = "../q4_TP_calc/q4_sig_avkin_thnq_pm/"
m_dir_pp = ""
m_file_pp = lambda theta: "csec_calc_theta_%0.2f_no_left_1_1_1_1.data"%(theta)
# PARIS FSI calculation
m_dir_pf = ""
m_file_pf = lambda theta: "csec_calc_theta_%0.2f_no_left_1_1_1_123.data"%(theta)

# PWIA misak CD Bonn
m_dir_cp = "../q4_TP_calc/q4_sig_avkin_thnq_pm/"
#m_dir_cp = ""
m_file_cp = lambda theta: "csec_calc_theta_pm_%0.2f_fsi_norad_avgkin_3_1_1_1.data"%(theta)

# FSI misak CD Bonn
m_dir_cf = "../q4_TP_calc/q4_sig_avkin_thnq_pm/"
#m_dir_cf = ""
m_file_cf = lambda theta: "csec_calc_theta_pm_%0.2f_fsi_norad_avgkin_3_1_1_12.data"%(theta)


# PWIA misak V18
m_dir_vp = "../q4_TP_calc/q4_sig_avkin_thnq_pm/"
#m_dir_vp = ""
m_file_vp = lambda theta: "csec_calc_theta_pm_%0.2f_fsi_norad_avgkin_2_1_1_1.data"%(theta)

# FSI misak V18
m_dir_vf = "../q4_TP_calc/q4_sig_avkin_thnq_pm/"
#m_dir_vf = ""
m_file_vf = lambda theta: "csec_calc_theta_pm_%0.2f_fsi_norad_avgkin_2_1_1_12.data"%(theta)

# y limits:
y_limits = (1.0000000000000001e-05, 10.0)
x_limits = (0., 1.1)
y_ticks = np.array([  1.00000000e-05,   1.00000000e-04,   1.00000000e-03,\
                   1.00000000e-02,   1.00000000e-01,   1.00000000e+00,\
                   1.00000000e+01])

x_ticks = np.array([ 0.1,  0.2,  0.3,  0.4,  0.5,  0.6, .7, .8, .9, 1.])

use_bin_corr = True
skip_legend = True
skip_last = True

#-------------------------------------
# set label font size
matplotlib.rc('xtick', labelsize=8)     
matplotlib.rc('ytick', labelsize=8)

matplotlib.rc('axes', linewidth = 1.)
matplotlib.rc('xtick.major', width = 1.)
matplotlib.rc('ytick.major', width = 1.)

# using tex?
matplotlib.rc('text', usetex = True)
#-------------------------------------

# save the figure

# plot all pm distributions scaled on a semilog plot
pl.plt.figure(figsize = (6.,7.5))
f_scale = 10.

theta_v = np.arange(20., 50., 10)
skip_legend = True
legends = []
plots = []
ylim_max = y_limits[1]
ylim_min = y_limits[0]
colors_v = take(color_cycle, stop = len(theta_v))
marker_v = take(marker_cycle, stop = len(theta_v))

text_off_x = 0.02
text_scale_y = 0.3

for i,th in enumerate(theta_v):
    f_scale *= 1.e-2
    legend_txt = '{:.0f}'.format(th) + r'$^\circ$'
    px, x0, y0 = plot_all(th, scale = f_scale, label = legend_txt, color = colors_v[i], marker = marker_v[i])
    if px is None:
        continue
    pl.plt.ylim((f_scale*ylim_min, ylim_max))
    if (px == None):
        f_scale *=10.
    else:
        pl.plt.text(x0 + text_off_x , y0 * text_scale_y, legend_txt, fontsize=8, bbox={'facecolor':'white','edgecolor':'none', 'pad':0})

# pl.plt.title(r'$\rho(p_m)$')
pl.plt.title(r'')
pl.plt.xlabel(r'$p_m$ GeV/c')
pl.plt.ylabel(r'$\sigma_{red}$')
# pl.plt.legend(loc = 'upper left', frameon = False, fontsize = 8. )
# std plot
#pl.plt.xlim((0., 1.2))
#pl.plt.ylim((1e-19, 2.))


# zoomed version
# high pm
#pl.plt.xlim((0.48, 1.2))
#pl.plt.ylim((1e-09, 1e-05))
# low pm
pl.plt.xlim((0., .2))
pl.plt.ylim((1e-07, 1))


pl.plt.subplots_adjust(left = 0.165,\
                right = 0.95,\
                bottom = 0.09,\
                top = 0.95,\
                hspace = 0.0,\
                wspace = 0.0)
pl.plt.show()
