#!/usr/bin/env python3
#-*- coding: utf-8 -*-

# This function uses cubic interpolation to replace NaNs
# @author: Amin Shakya, Interdisciplinary Center for Water Research (ICWaR), Indian Institute of Science (IISc)

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


from scipy.interpolate import PchipInterpolator
import numpy as np

def naninterp(X):
    """This function uses cubic interpolation to replace NaNs

    Args:
        X (_type_): _description_

    Returns:
        _type_: _description_
    """

    nan = np.nan
    
    ok = ~np.isnan(X)
    xp = ok.ravel().nonzero()[0] #Indices of xs with values
    fp = X[~np.isnan(X)]
    
    x  = np.isnan(X).ravel().nonzero()[0] #Indices of xs without values
    
    pchip = PchipInterpolator(xp,fp) #Initialize scipy PHCIP cubic interpolation
    X[np.isnan(X)] = pchip(x) #Interpolate Nan values in X
    
    return X
    

    