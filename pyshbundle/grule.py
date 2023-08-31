#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Created on Fri Jun 17 15:31:15 2022
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# SHBundle original Code

# This function computes Gauss base points and weight factors
# using the algorithm given by Davis and Rabinowitz in 'Methods
# of Numerical Integration', page 365, Academic Press, 1975.

# bp, wf = grule(n)

# IN:
#    n ....... number of base points required. 

# OUT:
#    bp ...... cosine of the base points
#    wf ...... weight factors for computing integrals and such 

#@author: Amin Shakya, Interdisciplinary Center for Water Research (ICWaR), Indian Institute of Science (IISc)

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
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

import numpy as np

def grule(n: int):
    """This function computes Gauss base points and weight factors
    using the algorithm-see Reference

    Args:
        n (int): number of base points required

    Returns:
        np.array: cosine of the base points
        np.array: weight factors for computing integrals and such
    
    References:
        1. 'Methods of Numerical Integration' by Davis and Rabinowitz, page 365, Academic Press, 1975.
    
    """
    bp = np.zeros((n,1))
    wf = bp
    iter = 2
    m = np.floor((n+1)/2)
    e1 = n * (n+1)
    
    
    mm = 4*m - 1
    t = (np.pi / (4*n + 2)) * np.arange(3,mm+4,4)
    nn = (1 - (1 - 1/n)/(8*n*n))
    x0 = nn * np.cos(t)
    
    
    for i in range(iter):
        pkm1 = 1
        pk = x0
        
        for kk in range(n-1):
            k = kk + 2
            t1 = x0 * pk
            pkp1 = t1 - pkm1 - (t1-pkm1)/k  + t1
            pkm1=pk
            pk=pkp1
            
        den = 1 - x0*x0
        d1 = n * (pkm1 - x0*pk)
        dpn = d1/den
        
        
        d2pn = (2*x0*dpn - e1*pk) / den
        d3pn = (4*x0*d2pn + (2-e1)*dpn)/den
        d4pn = (6*x0*d3pn + (6-e1)*d2pn)/den
        u = pk/dpn
        v = d2pn/dpn
        h = -u * (1+(.5*u)*(v+u*(v*v - u*d3pn/(3*dpn))))
        p = pk + h*(dpn+(0.5*h)*(d2pn+(h/3)*(d3pn + 0.25*h*d4pn)))
        dp = dpn + h*(d2pn+(0.5*h)*(d3pn+h*d4pn/3))
        h = h-p/dp
        x0 = x0+h
    
    bp = -x0-h
    fx = d1 - h*e1*(pk+(h/2)*(dpn+(h/3)*(d2pn+(h/4)*(d3pn+(0.2*h)*d4pn))))
    wf = [2 * (1 - np.power(bp,2))]/(fx*fx)
                 
                
    for i in range(len(bp),n):
        bp = np.append(bp,[0])
        wf = np.append(wf,[0])
    
    if ((m)+(m)) != (n):
        m = m-1
    
    for i in range(1,int(m+1)):
        bp[-i] = -bp[i-1]
        wf[-i] = wf[i-1] 
    return bp, wf