#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Created on Tue Jan 24 10:14:27 2023

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

# author: Vivek Yadav, Interdisciplinary Center for Water Research (ICWaR), Indian Institute of Science (IISc)
# document and new feature - Abhishek Mhamane, MS-R Geoinformatics, IIT Kanpur
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

import numpy as np

def klm2sc(data):
    """Converts the spherical harmonic coefficients from klm format to /S|C\ format

    Args:
        data (list): list containing [degree;  order; clm; slm; delta clm; delta slm; start data; end date]

    Returns:
        np.ndarray: Spherical Harmonic Coefficients in /S|C\ format [[files or months]; [2-D matrix of /S|C\ format]]
        np.ndarray: Standard Deviations of correcponding Spherical Harmonic Coefficients in /S|C\ format [[files or months]; [2-D matrix of /S|C\ format]]
    """
    # import pickle
    # with open("/path/saved_as_num", "rb") as pk:
    #     data=pickle.load(pk)
        
    ''' Read variables '''
    no_of_years = len(data[0])
    degree = data[0]
    clm = data[2]
    slm = data[3]
    clm_std_dev = data[4]
    slm_std_dev = data[5]
    
    lmax=degree[0][-1]
    degree_order=int((lmax+1)*(lmax+2)/2)
    
    ''' Count no of months of data '''
    month_count =0
    for i in range(0,len(data[0]),1):
        month_count= month_count+round(len(data[0][i])/degree_order)
        
    ''' klm >>> sc '''
    month = 0
    sc_mat = np.zeros([month_count,lmax+1,2*lmax+2])
    dev_sc_mat = np.zeros((month_count, lmax+1, 2*lmax + 2))

    for year in range(0,no_of_years,1):
        index2 =0
        for tile in range(0,int(len(clm[year])/degree_order),1):  
            for index1 in range(0,lmax+1,1):
                sc_mat[month,index1:,lmax-index1] = slm[year][(index2):(index2+lmax-index1+1)]
                sc_mat[month,index1:,index1+lmax] = clm[year][(index2):(index2+lmax-index1+1)]

                dev_sc_mat[month,index1:,lmax-index1] = slm_std_dev[year][(index2):(index2+lmax-index1+1)]
                dev_sc_mat[month,index1:,index1+lmax] = clm_std_dev[year][(index2):(index2+lmax-index1+1)]

        
                
                #print(month,'\t',index1,'\t',Lmax-index1,'\t',year,'\t',index2,'\t',index2+Lmax-index1+1)
                index2 = index2+lmax-index1+1
            month=month+1
    
    sc_mat=np.delete(sc_mat,lmax,axis=2)
    dev_sc_mat=np.delete(dev_sc_mat,lmax,axis=2)

    print('Conversion into clm format complete')
    return sc_mat, dev_sc_mat

def klm2sc_new(data_mat: np.ndarray, lmax: int, sigma_flag=False):
    """Converts the spherical harmonic coefficients from klm format to /S|C\ format

    Args:
        data_mat (np.ndarray): A 2-D matrix(numpy ndarray)
        lmax (int): maximum degree of spherical harmonic expansion
        sigma_flag (bool, optional): Flag to return the associated standard deviation values. Defaults to False.

    Returns:
        np.ndarray: Spherical Harmonic Coefficients or/and associated standard deviations in /S|C\ format
    """
    sc_mat = np.zeros((lmax+1, 2*lmax + 2))
    dev_sc_mat = np.zeros((lmax+1, 2*lmax + 2))
    clm = data_mat[:, 2]
    slm = data_mat[:, 3]
    clm_std_dev = data_mat[:, 4]
    slm_std_dev = data_mat[:, 5]
    
    # first place the slm and then clm
    index2 =0
    for index1 in range(0,lmax+1,1):
        sc_mat[index1:, lmax-index1] = slm[(index2):(index2 + lmax-index1+1)]
        sc_mat[index1:, index1+lmax] = clm[(index2):(index2 + lmax-index1+1)]

        dev_sc_mat[index1:, lmax-index1] = slm_std_dev[(index2):(index2 + lmax-index1+1)]
        dev_sc_mat[index1:, index1+lmax] = clm_std_dev[(index2):(index2 + lmax-index1+1)]
        
        index2 = index2 + lmax-index1+1

    sc_mat=np.delete(sc_mat,lmax,axis=1)
    dev_sc_mat=np.delete(dev_sc_mat,lmax,axis=1)

    if sigma_flag:
        return sc_mat, dev_sc_mat
    else: 
        return sc_mat