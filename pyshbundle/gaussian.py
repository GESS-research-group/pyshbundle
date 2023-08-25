#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Created on Thu Jun 30 13:47:51 2022
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Original codes

# gaussian(L, cap): The program delivers the spherical harmonic coefficients of a gaussian
# smoothing filter. The coefficients are calculates according to Wahr et.
# al. (1998) equation (34) and Swenson and Wahr equation (34)

# How:      Wl = gaussian(L,cap)

# Input:    L    integer    maximum degree
#          cap  integer    half width of Gaussian smoothing function [km]
# Output:   Wl   [L+1 x 1]  smoothing coefficients

#@author: Amin Shakya, Interdisciplinary Center for Water Research (ICWaR), Indian Institute of Science (IISc)

# License:
#    This file is part of PySHbundle.
#    PySHbundle is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# Acknowledgement Statement:
#    Please note that PySHbundle has adapted the following code packages, 
#    both licensed under GNU General Public License
#       1. SHbundle: https://www.gis.uni-stuttgart.de/en/research/downloads/shbundle/

#       2. Downscaling GRACE Total Water Storage Change using Partial Least Squares Regression
#          https://springernature.figshare.com/collections/Downscaling_GRACE_Total_Water_Storage_Change_using_Partial_Least_Squares_Regression/5054564 
    
# Key Papers Referred:
#    1. Vishwakarma, B. D., Horwath, M., Devaraju, B., Groh, A., & Sneeuw, N. (2017). 
#       A data‐driven approach for repairing the hydrological catchment signal damage 
#       due to filtering of GRACE products. Water Resources Research, 
#       53(11), 9824-9844. https://doi.org/10.1002/2017WR021150

#    2. Vishwakarma, B. D., Zhang, J., & Sneeuw, N. (2021). 
#       Downscaling GRACE total water storage change using 
#       partial least squares regression. Scientific data, 8(1), 95.
#       https://doi.org/10.1038/s41597-021-00862-6
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

import numpy as np

def gaussian(L: int, cap: int):
    """The program delivers the spherical harmonic coefficients of a gaussian
    smoothing filter. The coefficients are calculated according to Wahr et. al.(1998)
    equation (34) and Swenson and Wahr equation (34)


    Args:
        L (int): maximum degree
        cap (int): half width of Gaussian smoothing function [km]

    Returns:
        np.ndarray: smoothing coefficients
    
    Raises:
        TypeError: Degree must be integer
        ValueError: Maximum degree must be higher than 2
        TypeError: Cap size must be an integer
    
    References:
        Wahr et.al. (1998) equation (34) and Swenson and Wahr equation (34)
    
    """
    
    #Check input
    if type(L) != int:
        raise TypeError('Degree must be integer')
        
    if L<2:
        raise ValueError('Maximum degree must be higher than 2')
        
    if type(cap) != int:
        raise TypeError('Cap size must be an integer')
        
    #Calculations
    W = np.zeros([L+1, 1])
    b = np.log(2)/(1 - np.cos(cap/6371))
    
    #Recursive calculation of the weighting coefficients
    W[0,0] = 1
    W[1,0] = np.power(10, np.log10( (1 + np.exp(-2*b))/(1-np.exp(-2*b)) - (1/b)))
    
    i = 1
    while i < L:        
        j = i + 1
        W[i+1][0] = W[i-1][0] - (2*(j-1) + 1)/b * W[i][0]
        if W[i+1, 0] > W[i] or W[i+1] < 0:
            W[i+1] = 0
        i = i + 1
    
    return W

