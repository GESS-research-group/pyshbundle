#!/usr/bin/env python3
#-*- coding: utf-8 -*-

# Created on Thu Jun 30 14:34:55 2022
# NORMALKLM returns an ellipsoidal normal field
# consisting of normalized -Jn, n=0,2,4,6,8

# IN:
#    lmax ....... maximum degree
#    typ ........ either 'wgs84' (equipotential ellipsoid), default,
#                        'grs80',
#                 or     'he' (hydrostatic equilibrium ellipsoid)
# OUT:
#    nklm ....... normal field in CS-format (sparse)

# REMARKS:
#    .J2,J4 values for hydrostatic equilibrium ellipsoid from Lambeck (1988)
#    "Geophysical Geodesy", p.18
   
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
from scipy import sparse


def normalklm(lmax: int, typ: str = 'wgs84'):
    """ NORMALKLM returns an ellipsoidal normal field
    consisting of normalized -Jn, n=0,2,4,6,8

    Args:
        lmax (int): maximum degree
        typ (str): Ellipsoids can be either 
                    'wgs84' - World Geodetic System 84, 
                    'grs80' - , 
                    'he' - hydrostatic equilibrium ellipsoid
    
    Returns:
        nklm (np.array): normal field in CS-format (sparse array - [1, -J2, -J4, -J6, -J8])
    
    TODO: 
        Find type of nklm; I think raising TypeError, VlueError or NameError instad of general Exception

    Raises:
        TypeError: lmax should be an integer
        ValueError: lmax should be positive
        ValueError: Unknown type of ellipsoid, supports 'wgs84', `GRS80` and 'he'
    
    References:
        1. J2,J4 values for hydrostatic equilibrium ellipsoid from Lambeck (1988)
        "Geophysical Geodesy", p.18    
    """
    
    if type(lmax) != int:
        raise TypeError("lmax should be integer")
        
    if lmax < 0:
        raise ValueError("lmax should be positive")
        
    
    typ_ = typ.lower()
    if (typ_ == 'wgs84'):
        J2     =  1.08262982131e-3     #% earth's dyn. form factor (= -C20 unnormalized)
        J4     = -2.37091120053e-6    #% -C40 unnormalized
        J6     =  6.08346498882e-9     #% -C60 unnormalized
        J8     = -1.42681087920e-11    #% -C80 unnormalized
        jcoefs = np.array([1, -J2, -J4, -J6, -J8]).T.reshape(5,1)
        # as lmax + 2 is requires 
        l      = np.arange(0,min(lmax + 2,8 + 2), 2).T
        l.reshape(l.shape[0],1)
        
    elif (typ_ == 'grs80'):
        J2     =  1.08263e-3         # % earth's dyn. form factor (= -C20 unnormalized)
        J4     = -2.37091222e-6     #% -C40 unnormalized
        J6     =  6.08347e-9        #% -C60 unnormalized
        J8     = -1.427e-11         #% -C80 unnormalized
        jcoefs = np.array([1, -J2, -J4, -J6, -J8]).reshape(5,1)
        l      = np.arange(0,min(lmax + 2,8 + 2), 2).T
        l.reshape(l.shape[0],1)
        
    elif ((typ_ == 'he') or (typ_ == 'hydro')):
        J2     = 1.072618e-3		#% earth's dyn. form factor (= -C20 unnormalized)
        J4     = 0.2992e-5     	#% -C40 unnormalized
        jcoefs = np.array([1, -J2, -J4]).T.reshape(5,1)
        # adding (2) beacuse of arange function is only uptp last integer and not including last
        l      = np.arange(0,min(lmax + 2,4 + 2), 2).T
        l.reshape(l.shape[0],1)
        
    else:
        raise ValueError("Unknown type of ellipsoid:   ", typ)
    
    coefs = jcoefs[:len(l)].T / np.sqrt(2*l + 1)
#    coefs.reshape(coefs.shape[0],1)
    
    
    data = np.array(coefs)[0]
    row = np.array(l)
    col = np.zeros(len(l))
    # lmax = 96 then shape=(97, 97) -> consisitent with everything else
    nklm = sparse.coo_matrix((data,(row,col)),shape=(lmax+1,lmax+1)).toarray()
    return nklm