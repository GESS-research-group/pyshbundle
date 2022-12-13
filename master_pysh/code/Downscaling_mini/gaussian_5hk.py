#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 13:47:51 2022
% GAUSSIAN delivers the spherical harmonic coefficients of a gaussian
% smoothing filter. The coefficients are calculates according to Wahr et.
% al. (1998) equation (34) and Swenson and Wahr equation (34)
%
% How:      Wl = gaussian(L,cap)
%
% Input:    L    integer    maximum degree
%           cap  integer    half width of Gaussian smoothing function [km]
%
% Output:   Wl   [L+1 x 1]  smoothing coefficients
@author: Amin Shakya
"""

def gaussian(L, cap):
    import numpy as np
    
    #Check input
    if type(L) != int:
        raise Exception('Degree must be integer')
        
    
    if L<2:
        raise Exception('Maximum degree must be higher than 2')
        
    
    if type(cap) != int:
        raise Exception('Cap size must be an integer')
        
    
    #Calculations
    W = np.zeros([L+1, 1])
    b = np.log(2)/(1 - np.cos(cap/6371))
    
    #Recursive calculation of the weighting coefficients
    W[0,0] = 1
    W[1,0] = np.power(10, np.log10( (1+np.exp(-2*b))/(1-np.exp(-2*b)) - (1/b)))
    
    i = 1
    while i<L:        
        j = i + 1
        W[i+1][0] = W[i-1][0] - (2 * (j-1) + 1)/b * W[i][0]
        if W[i+1,0] > W[i] or W[i+1] < 0:
            W[i+1] = 0
        i = i + 1
    
    return W

