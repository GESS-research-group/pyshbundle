#Created on Fri Dec  9 11:43:15 2022
# Script contains various methods to reshape the extracted SH coefficients from different formats to the desired format

#@author1: Abhishek Mhamane, MS-R Geoinformatics, IIT Kanpur, India
#@author2: Vivek Yadav, Interdisciplinary Center for Water Research (ICWaR), Indian Institute of Science (IISc)
# Updated, Vivek, 2024-10-04

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

def sc2cs(field):
    """
    Converts SH coefficients from SC to CS format.

    Converts the rectangular (L+1) x (2L+1) matrix `field`, containing
    spherical harmonics coefficients in SC storage format, into a 
    square (L+1) x (L+1) matrix in CS format.

    Args:
        field (numpy.ndarray): The rectangular (L+1) x (2L+1) matrix, containing the
            spherical harmonics coefficients in SC storage format.
        
    Returns:
        (numpy.ndarray): Square (L+1) x (L+1) matrix in CS format.
    
    References:
        See the SHBundle docs or PySHBundle docs for more info about SH coefficient storage and retrieval formats being implemented.

    Examples:
        sc2cs(field)
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



def clm2cs(data_mat: np.ndarray, lmax: int, sigma_flag=False):
    """
    Converts the format from CLM to CS.

    Under the hood uses the `clm2sc` and `sc2cs` functions.

    Args:
        data_mat (numpy.ndarray): List containing [degree, order, clm, slm, delta clm, delta slm, start date, end date].
        lmax (int): Max Degree of the spherical harmonic expansion.
        sigma_flag (bool, optional): Flag to return the standard deviation data in CS format or not. Defaults to False.

    Returns:
        (numpy.ndarray): Spherical Harmonic Coefficients in CS format.
    """
    if sigma_flag:
        sc_mat, dev_sc = clm2sc(data_mat=data_mat, lmax=lmax, sigma_flag=True)
        return sc2cs(sc_mat), sc2cs.sc2cs(dev_sc)
    else:
        sc_mat = clm2sc(data_mat=data_mat, lmax=lmax, sigma_flag=False)
        return sc2cs(sc_mat)

    

def clm2sc(data_mat: np.ndarray, lmax: int, sigma_flag=False):
    """
    Converts the spherical harmonic coefficients from CLM format to SC format.

    Args:
        data_mat (numpy.ndarray): List containing [degree, order, clm, slm, delta clm, delta slm, start date, end date].
        lmax (int): Max Degree of the spherical harmonic expansion.
        sigma_flag (bool, optional): Flag to return the standard deviation data in CS format or not. Defaults to False.

    Returns:
        (numpy.ndarray): Spherical Harmonic Coefficients in SC format.

    References:
        Refer to the SHBundle or PySHBundle docs for the different data storage and retrieval formats.
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
    """
    Converts SH coefficients from CS to SC format.

    Converts the square (L+1)x(L+1) matrix `field`, containing
    spherical harmonics coefficients in CS storage format, into a 
    rectangular (L+1)x(2L+1) matrix in SC format.

    Args:
        field (numpy.ndarray): The square (L+1)x(L+1) matrix, containing
            spherical harmonics coefficients in CS storage format.
    
    Returns:
        (numpy.ndarray): Rectangular (L+1)x(2L+1) matrix in SC format.

    Raises:
        TypeError: If the input is neither in CS nor in SC format.
    
    Examples:
        cs2sc(field)
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


def klm2sc(data_mat: np.ndarray, lmax: int, sigma_flag=False):
    """
    Converts the spherical harmonic coefficients from klm format to SC format.

    Args:
        data_mat (numpy.ndarray): List containing [degree, order, clm, slm, delta clm, delta slm, start date, end date].
        lmax (int): Max Degree of the spherical harmonic expansion.
        sigma_flag (bool, optional): Flag to return the standard deviation data in CS format or not. Defaults to False.

    Returns:
        (numpy.ndarray): Spherical Harmonic Coefficients or/and associated standard deviations in SC format.
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

def cklm2sc_new(clm_mat, lmax: int):
    """
    Transforms the spherical harmonics coefficients data in clm or klm format into a SC matrix.

    Args:
        clm_mat (numpy.ndarray): The input matrix containing spherical harmonics coefficients.
        lmax (int): The maximum degree of the spherical harmonic expansion.

    Returns:
        tuple: A tuple containing:
            - scmat (numpy.ndarray): The SC matrix.
            - dev_scmat (numpy.ndarray): The deviation SC matrix.
    """

    # initialise an empty sc matrix
    sc_mat = np.zeros([lmax+1, 2*lmax + 1])
    dev_sc_mat = np.zeros([lmax+1, 2*lmax + 1])

    # navigating the maze of indices
    
    # Use logical indices

    # sc mat requires padding - Taken care of by the earlier initialisation
    # 
    # filling the value at appropriate locaation is the key
    # 
    # Approach-1
        # run through rows(degree) and fill the cols(order) respectively

    # Approach -2
        # create a row_func s.t. [....., C22, C21, C20, S21, S22, .......]
        # then stack the rows
    
    # First flatten the SC matrix - column wise aka Fortran style
    # get the flattented idx to be raplaced using sub2ind 
    # replace the indices at those locations using 
    # unflatten the matrix

    shape_sc = sc_mat.shape

    # following the approach similar to Octave implementation
    # using matrix operations to improve the time efficiency as compared to looping
    idx_s = sub2ind(sc_mat.shape, clm_mat[:, 0].astype('i'), (lmax - clm_mat[:, 1]).astype('i')).astype('i')
    idx_c = sub2ind(sc_mat.shape, clm_mat[:, 0].astype('i'), (lmax + clm_mat[:, 1]).astype('i')).astype('i')

    
    flat_sc = sc_mat.flatten("F")
    # Attention first place the slm coeff. or else it will relace zonal clm coeff.
    flat_sc[idx_s] = clm_mat[:, 3]
    flat_sc[idx_c] = clm_mat[:, 2]

    flat_sc2 = dev_sc_mat.flatten("F")
    flat_sc2[idx_s] = clm_mat[:, 5]
    flat_sc2[idx_c] = clm_mat[:, 4]

    dev_scmat = flat_sc2.reshape(shape_sc)

    scmat = flat_sc.reshape(shape_sc)

    # with one flag include for 

    return scmat, dev_scmat