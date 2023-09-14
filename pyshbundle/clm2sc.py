#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Created on Fri Dec  9 10:14:32 2022

# author: Vivek Yadav, Interdisciplinary Center for Water Research (ICWaR), Indian Institute of Science (IISc)
# Debug, documentation and new feature: Abhishek Mhamane, MS-Research Geoinformatics, Indian Institute of Technology (IIT) Kanpur

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
from tqdm import tqdm

# CLM to /S|C\
def clm2sc(data):
    """Converts the spherical harmonic coefficients from clm format to /S|C\ format

    Args:
        data (list): list containing [degree;  order; clm; slm; delta clm; delta slm; start data; end date]

    Returns:
        np.ndarray: Spherical Harmonic Coefficients in /S|C\ format [[files or months]; [2-D matrix of /S|C\ format]]
        np.ndarray: Standard Deviations of correcponding Spherical Harmonic Coefficients in /S|C\ format [[files or months]; [2-D matrix of /S|C\ format]]
    
    References:
        Refer to the SHBundle or PySHBundle docs for the different data storage and retrival formats.
    """
    # import pickle
    # data=np.load("/path/saved_as-num.npy",allow_pickle=True)
    # with open("/path/saved_as_num", "rb") as pk:
    #     data=pickle.load(pk)

    # Read variables 
    no_of_years = len(data[0])
    degree = data[0]
    clm = data[2]
    slm = data[3]
    clm_std_dev = data[4]
    slm_std_dev = data[5]

    Lmax = degree[0][-1]
    degree_order = int((Lmax+1) * (Lmax+2)/2)
    # Count no of months of data 
    month_count = 0
    for i in range(0, len(data[0]), 1):
        month_count = month_count+round(len(data[0][i])/degree_order)
    
    # clm >>> sc 
    month = 0
    sc_mat = np.zeros([month_count, Lmax+1, 2*Lmax + 2])
    dev_sc_mat = np.zeros([month_count, Lmax+1, 2*Lmax + 2])
    for year in range(0, no_of_years, 1):
        for tile in range(0,int(len(degree[year])/degree_order), 1):
            i = 0
            for index1 in range(0,Lmax+1, 1):
                for index2 in range(0,index1+1, 1):
                    
                    sc_mat[month, index1, Lmax-index2] = slm[year][i + tile*degree_order]
                    sc_mat[month, index1, Lmax+index2+1] = clm[year][i + tile*degree_order]

                    dev_sc_mat[month, index1, Lmax-index2] = slm_std_dev[year][i + tile*degree_order]
                    dev_sc_mat[month, index1, Lmax+index2+1] = clm_std_dev[year][i + tile*degree_order]

                    i = i + 1
            month = month + 1
            
    # delete order 0 column
    sc_mat = np.delete(sc_mat, Lmax, 2)
    dev_sc_mat = np.delete(dev_sc_mat, Lmax, 2)
    print('Conversion into clm format complete')
    return sc_mat, dev_sc_mat

# mean=np.mean(sc_mat[18:102], axis=0)


def clm2sc_new(data_mat: np.ndarray, lmax: int, sigma_flag=False):
    """Converts the spherical harmonic coefficients from clm format to /S|C\ format

    Args:
        data_mat (numpy.ndarray): list containing [degree;  order; clm; slm; delta clm; delta slm; start data; end date]
        lmax (int): Max Degree of the spherical harmonic expansion
        sigma_flag (boolean): Flag to return the standard deviation data in /S|C\ format or not. Defaults to False

    Returns:
        numpy.ndarray: Spherical Harmonic Coefficients in /S|C\ format
    
    References:
        Refer to the SHBundle or PySHBundle docs for the different data storage and retrival formats.
    
    """

    

    sc_mat = np.zeros((lmax+1, 2*lmax + 2))
    dev_sc_mat = np.zeros((lmax+1, 2*lmax + 2))
    
    # as per the convention
    clm = data_mat[:, 2]
    slm = data_mat[:, 3]
    clm_std_dev = data_mat[:, 4]
    slm_std_dev = data_mat[:, 5]

    i = 0
    for index1 in range(0,lmax+1, 1):
        for index2 in range(0,index1+1, 1):
            
            sc_mat[index1, lmax-index2] = slm[i]
            sc_mat[index1, lmax+index2+1] = clm[i]

            dev_sc_mat[index1, lmax-index2] = slm_std_dev[i]
            dev_sc_mat[index1, lmax+index2+1] = clm_std_dev[i]

            i = i + 1
    
    sc_mat = np.delete(sc_mat, lmax, 1)
    dev_sc_mat = np.delete(dev_sc_mat, lmax, 1)

    if sigma_flag:
        return sc_mat, dev_sc_mat
    else:
        return sc_mat