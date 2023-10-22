#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Created on Fri Dec  9 11:43:15 2022

#@author: Vivek Yadav, Interdisciplinary Center for Water Research (ICWaR), Indian Institute of Science (IISc)

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
from pyshbundle import sc2cs, clm2sc

# CLM to C|S format
def clm2cs(data):
    """_summary_

    Args:
        data (_type_): _description_
    
    Suggestion: 
        Instead of printing "conversion complete" let's show a progress bar
    """
    
    # Load data 
    # data = np.load(path, allow_pickle=True)
    # data needs to be loaded in numpy format
    
    # Read variables 
    no_of_years = len(data[0])
    degree = data[0]
    clm = data[2]
    slm = data[3]
    
    # Count no of months of data 
    month_count =0
    for i in range(0,len(data[0]),1):
        month_count= month_count+len(data[0][i])/4750
    # clm >>> cs 
    month = 0
    Lmax = degree[0][-1]
    cs_mat = np.zeros([int(month_count),Lmax+1,Lmax+1])
    for year in range(0,no_of_years,1):
        for tile in range(0,int(len(clm[year])/4752),1):
            i,j = 0,0 
            for index1 in range(2,Lmax+1,1):
                for index2 in range(0,index1+1,1):
                    cs_mat[month,index1,index2] = clm[year][i + tile*4752]
                    i = i + 1
            for index3 in range(2,Lmax+1,1):
                for index4 in range(0,index3,1):
                    cs_mat[month,index4,index3] = slm[year][j+1 +tile*4752]
                    j = j + 1
                j = j + 1
            month = month + 1
    print('Conversion into clm format complete')        
    #np.save('/path/SH_coeff_cs.npy', cs_mat)
    return cs_mat

def clm2cs_new(data_mat: np.ndarray, lmax: int, sigma_flag=False):
    """Converts the format from CLM to |C\S|
    Under the hood uses the `clm2sc` and `sc2cs` function

    Args:
        data_mat (numpy.ndarray): list containing [degree;  order; clm; slm; delta clm; delta slm; start data; end date]
        lmax (int): Max Degree of the spherical harmonic expansion
        sigma_flag (boolean): Flag to return the standard deviation data in |C\S| format or not. Defaults to False

    Returns:
        numpy.ndarray: Spherical Harmonic Coefficients in |C\S| format
        
    """
    if sigma_flag:
        sc_mat, dev_sc = clm2sc.clm2sc(data_mat=data_mat, lmax=lmax, sigma_flag=True)
        return sc2cs.sc2cs(sc_mat), sc2cs.sc2cs(dev_sc)
    else:
        sc_mat = clm2sc.clm2sc(data_mat=data_mat, lmax=lmax, sigma_flag=False)
        return sc2cs.sc2cs(sc_mat)

    
