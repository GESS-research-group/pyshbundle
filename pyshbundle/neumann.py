#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Created on Sat Jun 18 14:38:16 2022
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# SHBundle Original Codes

#% NEUMANN returns the weights and nodes for Neumann's numerical integration
#%
#% NB1 In 1st N-method, length(x) should not become too high, 
#% since a linear system of this size is solved. Eg: 500.
#% NB2 No use is made of potential symmetries of nodes.
#%
#% HOW:   w     = neumann(x)   -- 1st Neumann method 
#%        [w,x] = neumann(n)   -- 2nd Neumann method (Gauss quadrature)
#%
#% IN:
#%    x ...... base points (nodes) in the interval [-1;1]
#%    n ...... number of weights and number of base points
#%
#% OUT:
#%    w ...... quadrature weights
#%    x ...... base points (nodes) in the interval [-1;1]
#%
#% USES:
#%    plm, uberall/grule
#%
#% REMARKS:
#%    1st N.-method: see Sneeuw (1994) GJI 118, pp 707-716, eq. 19.5
#%    2nd N.-method: see uberall/GRULE
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

import numpy as np
from pyshbundle import grule
from pyshbundle import plm

def neumann(inn):
    """Returns the weights and nodes for Neumann's numerical integration

    Args:
        inn (int, np.array): base points (nodes) in the interval [-1;1]

    Raises:
        TypeError: Integer input argument required
        ValueError: Error in input dimensions

    Returns:
        _type_: quadrature weights
        _type_: base points (nodes) in the interval [-1;1]
    
    Remarks:
        * 1st N.-method: see Sneeuw (1994) GJI 118, pp 707-716, eq. 19.5
        * 2nd N.-method: see uberall/GRULE
    
    Todo: 
        + TypeError is more relavant and shape error from numpy
    
    Uses:
        `grule`, `plm`

    Examples:
        >>> TO DO: write example how to use the function
    """

    try: #if input is an integer
        x, w = grule.grule(inn)
    except: #if input is an array
        if(len(inn)==1): #2nd Neumann method
            x, w = grule.grule(inn)
            if(np.not_equal(np.mod(x, 1), 0)): #Not integer
                raise TypeError("Integer input argument required")
            
            
        
        elif min(inn.shape) == 1: #1st Neumann method #Size gives 2 outputs for 2d array in matlab; for row and column
            x = inn
            theRAD = np.arccos(x) #x in radian
            l = np.array(list(range(len(x))))
            pp = plm.plm(l, theRAD)
            
            rr = list([2])
            for i in len(x-1):
                rr.append(0)
            r = np.asarray(rr)
                
            w,resid,rank,s = np.linalg.lstsq(pp,r) #Solve system of equations; Double check this operation
            if(x.shape != w.shape):
                w = w.T
            
        else:
            raise ValueError("Error in input dimensions")
            # TO DO: Write more descriptive exception messages
    
    return w, x
        

