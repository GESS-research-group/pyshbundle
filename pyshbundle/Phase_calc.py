#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 09:10:27 2022
%UNTITLED3 Summary of this function goes here
% calculates the phase difference between two time series based on the
% hilbert transform method explained by Phillip et al.

% Phillips, T., R. S. Nerem, B. Fox-Kemper, J. S. Famiglietti, and B. Rajagopalan (2012),
% The influence of ENSO on global terrestrial water storage using GRACE, Geophysical
% Research Letters, 39 (16), L16,705, doi:10.1029/2012GL052495.
%--------------------------------------------------------------------------------
% written by Bramha Dutt Vishwakarma, Institute of Geodesy, University of
% Stuttgart.      ----- 12 July 2015 ----
%--------------------------------------------------------------------------------%%
@author: Amin Shakya
"""

def Phase_calc(fts, ffts):
    import scipy as sc
    from scipy import signal as signal
    import numpy as np
    
    c = fts.shape[1]
    
    ps = np.zeros((1, c))
    
    filter_ = ~np.isnan(fts)
    filter__ = ~np.isnan(ffts)
    
    fts_ = fts[filter_] #Extract values and leave Nan
    ffts_ = ffts[filter__] #Extract values and leave Nan
    
    fts = fts_.reshape(int(fts_.shape[0]/c),c)
    ffts = ffts_.reshape(int(ffts_.shape[0]/c),c)
    
    rn = fts.shape[0]
    
    for i in range(c):
        # A = np.concatenate(np.ones((rn,1)), np.real(signal.hilbert(ffts[:, i])), np.imag(signal.hilbert(ffts[:, i]))) #design matrix
        
        A = np.array((np.ones((rn)), np.real(signal.hilbert(ffts[:, i])), np.imag(signal.hilbert(ffts[:, i])))).T
        
        A = A.astype('double')
        B = fts[:,i]
        B = B.astype('double')
        abc = np.linalg.lstsq(A.T @ A, A.T @ B)[0]
        
        ps[0,i] = np.arctan2(abc[3-1],abc[2-1])*(180/np.pi) #check indices and degree/radian
    return ps