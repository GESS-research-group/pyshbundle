#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Created on Thu Jul 14 09:10:27 2022

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Summary of this function goes here
# calculates the phase difference between two time series based on the
# Hilbert transform method explained by Phillip et al.

# Phillips, T., R. S. Nerem, B. Fox-Kemper, J. S. Famiglietti, and B. Rajagopalan (2012),
# The influence of ENSO on global terrestrial water storage using GRACE, Geophysical
# Research Letters, 39 (16), L16,705, doi:10.1029/2012GL052495.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

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
   
# @author: Amin Shakya, Interdisciplinary Center for Water Research (ICWaR), Indian Institute of Science (IISc)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

import scipy as sc
from scipy import signal as signal
import numpy as np

def Phase_calc(fts, ffts):
    """calculates the phase difference between two time series based on the
    Hilbert transform method explained by Phillip et al.

    Args:
        fts (np.ndarray): _description_
        ffts (np.ndarray): _description_

    Returns:
        _type_: _description_
    
    References:
        1. Phillips, T., R. S. Nerem, B. Fox-Kemper, J. S. Famiglietti, and B. Rajagopalan (2012),
        The influence of ENSO on global terrestrial water storage using GRACE, Geophysical
        Research Letters, 39 (16), L16,705, doi:10.1029/2012GL052495.
    """
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