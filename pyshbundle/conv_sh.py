#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Created on Fri Dec  9 11:43:15 2022

#@author1: Abhishek Mhamane, MS-R Geoinformatics, IIT Kanpur, India
#@author2: Vivek Yadav, Interdisciplinary Center for Water Research (ICWaR), Indian Institute of Science (IISc)

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


def clm2cs(data_mat: np.ndarray, lmax: int, sigma_flag=False):
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
        sc_mat, dev_sc = clm2sc(data_mat=data_mat, lmax=lmax, sigma_flag=True)
        return sc2cs(sc_mat), sc2cs.sc2cs(dev_sc)
    else:
        sc_mat = clm2sc(data_mat=data_mat, lmax=lmax, sigma_flag=False)
        return sc2cs(sc_mat)

    

def clm2sc(data_mat: np.ndarray, lmax: int, sigma_flag=False):
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


def cs2sc(field):
    """converts the square (L+1)x(L+1) matrix 'field', containing
    spherical harmonics coefficients in |C\S| storage format into a 
    rectangular (L+1)x(2L+1) matrix in  /S|C\ format.

    Args:
        field (np.ndarray): the square (L+1)x(L+1) np matrix field , containing
                   spherical harmonics coefficients in |C\S| storage format
    
    Returns:
        np.ndarray: Rectangular (L+1)x(2L+1) np matrix in  /S|C\ format

    Raises:
        TypeError: Input neither in cs nor in sc format
    
    Todo:
        + Rather use TypeError instead of base Exception
    
    Examples:
        >>> sc_shcoeff = cs2sc(cs_shcoeff)
        TO DO: write the output
    """
    rows = len(field)
    cols = len(field[0])

    if (rows != cols) and (cols != 2*rows - 1):
        raise TypeError("Input neither in cs nor in sc format")
    elif cols == 2*rows - 1:
        sc = field
    else:
        c    = np.tril(field)
        ut   = np.triu(field)
        i = np.identity(rows)
        i = 1-i
        s    = np.fliplr(np.transpose(np.multiply(ut, i, )))
        sc   = np.concatenate((s[:,1:rows], c), axis=1)
        
    return(sc)


def sc2cs(field):
    """converts the rectangular $(L+1) \times (2L+1)$ matrix FIELD, containing
    spherical harmonics coefficients in /S|C\ storage format into a 
    square (L+1)x(L+1) matrix in |C\S| format.

    Parameters:
        field (np.ndarray()):
            the rectangular (L+1)x(2L+1) matrix FIELD, containing
            spherical harmonics coefficients in /S|C\ storage format
        
    Returns: 
        cs (np.ndarray): 
            square (L+1)x(L+1) matrix in |C\S| format
    
    References:
        See the SHBundle docs or PySHBundle docs for more info about SH coeff. storage and retrival formats being implementd.

    Examples:
        >>> cs_fmt = sc2cs(field)
        TO DO: show suitable output
    """

    rows = len(field)
    cols = len(field[0])

    if (rows!=cols) and (cols!=2*rows - 1):
        sc2cs.exit("Input neither in cs nor in sc format")
    elif cols == rows:
        cs = field
    else:
        c    = field[:, rows-1:cols]
        st   = np.transpose(np.fliplr(field[:, 0:rows-1]))
        z    = np.zeros([1,rows])
        s    = np.concatenate((st, z), axis=0)
        cs   = np.add(c, s)
        
    return(cs)


def klm2sc(data_mat: np.ndarray, lmax: int, sigma_flag=False):
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