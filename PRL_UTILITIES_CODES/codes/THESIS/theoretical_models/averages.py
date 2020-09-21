#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 10:46:02 2019

Tools for calculating averages

@author: boeglinw
"""

import numpy as np

#---------------------------------------------------------------------- 
def get_bin_id(d):
    kins = d['kin']
    ix = np.array(d['ix']).astype(int)
    iy = np.array(d['iy']).astype(int)
    id_s = []
    for i, k in enumerate(kins): 
        id_s.append( '{:d}_{:d}_{:s}'.format( ix[i], iy[i], k.strip()))
    return np.array(id_s)
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
    # calculate the averages for the given arrays, where there are duplicate pm values
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
 #----------------------------------------------------------------------
 
 #----------------------------------------------------------------------
   
def get_matched_array_average(pm, sig, dsig, same = True, match = None):
    # calculate the averages for the given arrays, where there are duplicate pm values
    # to use another array for matching: same = False and match = 'other array'
    #
    # 
    # example : get_matched_array_average(d35['pr'], d35['crs0'], d35['crs0']*0.05, same = False, match = d35['pm_b'])
    #
    # calculates the average for all common values of pm_b but calculates the average of the corresponding pr values
    #
    if same:
        m_av, m_ind = find_duplicates(pm)
    else:
        m_av, m_ind = find_duplicates(match)
    pm_av_c = np.array([])
    sig_av = np.array([])
    dsig_av = np.array([])
    dsig1_av = np.array([])
    for p in m_ind:
        pma, siga, dsiga,  dsiga1 = get_averages (pm[p], sig[p], dsig[p] )
        pm_av_c = np.append(pm_av_c, pma)
        sig_av = np.append(sig_av, siga)
        dsig_av = np.append(dsig_av, dsiga)
        # error including chisq
        dsig1_av = np.append(dsig1_av, dsiga1)
    return pm_av_c, sig_av, dsig_av, dsig1_av
 #----------------------------------------------------------------------