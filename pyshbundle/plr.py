#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 14:55:09 2022

%
% Input: L [t x s]
%        S [t x n]
%        r  number of selected modes

        aw = area vector for grid cells belonging to the catchment 
        nG = no of GRACE products
        
% Output: H [s x n] model parameters
%


This file is part of PySHbundle. 
    PySHbundle is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

Acknowledgement Statement:
    Please note that PySHbundle has adapted the following code packages, 
    both licensed under GNU General Public License
    1. SHbundle: https://www.gis.uni-stuttgart.de/en/research/downloads/shbundle/

    2. Downscaling GRACE Total Water Storage Change using 
    Partial Least Squares Regression
    https://springernature.figshare.com/collections/Downscaling_GRACE_Total_Water_Storage_Change_using_Partial_Least_Squares_Regression/5054564 
    
Key Papers Referred:
    1. Vishwakarma, B. D., Horwath, M., Devaraju, B., Groh, A., & Sneeuw, N. (2017). 
    A data‐driven approach for repairing the hydrological catchment signal damage 
    due to filtering of GRACE products. Water Resources Research, 
    53(11), 9824-9844. https://doi.org/10.1002/2017WR021150

    2. Vishwakarma, B. D., Zhang, J., & Sneeuw, N. (2021). 
    Downscaling GRACE total water storage change using 
    partial least squares regression. Scientific data, 8(1), 95.
    https://doi.org/10.1038/s41597-021-00862-6 

@author: Amin Shakya
"""

def plr(L, S, r, aw, nG):
    import numpy as np
    from scipy import linalg
    
    t = L.shape[0]
    
    #Covariance C [s x n]
    C = np.dot(L.T, S/t)
    
    U, s, Vh = linalg.svd(C) #svd is a inbuilt function in Matlab
    V = Vh.T
    sum_S = np.sum(s)
    
    
    U = U[:, :r]
    sigma = np.diag(s)
    sigma = sigma[:r,:r]
    V = V[:, :r]

    ratio_Sigma = np.sum(s) / sum_S
    
    Ul = np.dot(L,U)
    K = np.linalg.lstsq(Ul, S)[0]
    
    Ut = np.vstack([Ul, L@U]) #Doublecheck this line
    P = np.array(S, (L[: -nG] @ aw.T) / (aw.T @ aw)) #Double-check aw definition
    Kt = np.linalg.lstsq(Ut, P)
    
    #Ht [s x n]
    
    H = np.dot(U, Kt)
    
    #Error covariance #Check order of operation
    Css = np.dot(S.T, S / t)
    Cll = np.dot(L.T, L / t)
    Csl = np.dot(S.T, L / t)
    Cee = Css - Csl @ U @ np.invert(U.T @ Cll @ U) @ U.T @ Csl.T
    
    return H, ratio_Sigma, Cee