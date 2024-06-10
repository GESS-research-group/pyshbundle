#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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

import numpy
from pyshbundle import GRACEconstants as GC


def upwcon(degree: int, height):
    """Returns the upward continuation $(R/r)^l$

    Args:
        degree (int): Spherical harmonic degree
        height (int): Height above mean Earth radius [m] [scalar/vector]
    
    Returns:
        uc (_type_): Upward continuation terms
    
    Uses:
        `GRACEconstants.GC`
    
    Todo:
        + Add input checking functionality and raise exceptions
        + Add reference to formula
    """
    # Created on Sat May  9 18:49:45 2022
    rr = numpy.divide(GC.ae, numpy.add(GC.ae,height))
    uc = numpy.power(rr, degree)

    return(uc)    

def lovenr(lmax: int):
    """
    Created on Mon May 11 11:09:28 2022

    Todo:
        + Add type and input checking functionality
    
    _author_: Dr. Bramha Dutt Vishwakarma, Interdisciplinary Center for Water Research (ICWaR), Indian Institute of Science (IISc)
    """
    l  = [0,  1,    2,    3,    4,    5,   6,   7,   8,   9,  10,  12,  15,  20,  30,  40,  50,  70, 100, 150, 200]
    kl = numpy.divide([0,270,-3030,-1940,-1320,-1040,-890,-810,-760,-720,-690,-640,-580,-510,-400,-330,-270,-200,-140,-100, -700],1e4)
    n = range(0, lmax+1, 1)
    kn = numpy.interp(n,l,kl)
    return(kn)

def lovenrPREM(lmax:int, frame):
    """
    Created on Mon May 11 11:51:29 2022
    
    @author:  Dr. Bramha Dutt Vishwakarma, Interdisciplinary Center for Water Research (ICWaR), Indian Institute of Science (IISc)
    """
    
    data = numpy.array([[ 1,  -0.28476,   0.00000,   0.10462],
    [2,  -0.99297,  -0.61274,   0.04661] ,
    [3,  -1.05142,  -0.58897 ,  0.21048] ,
    [4,  -1.05378,  -0.53513 ,  0.23564] ,
    [5,  -1.08658,  -0.52382 ,  0.23186] ,
    [6,  -1.14404,  -0.54222 ,  0.23263] ,
    [7,  -1.21254,  -0.57464 ,  0.24058] ,
    [8,  -1.28403,  -0.61256 ,  0.25308] ,
    [9,  -1.35479,  -0.65203 ,  0.26799] ,
    [10,  -1.42330,  -0.69140,   0.28419] ,
    [11,  -1.48909,  -0.72998,   0.30121] ,
    [12,  -1.55204,  -0.76749,   0.31880] ,
    [13,  -1.61221,  -0.80381,   0.33684] ,
    [14,  -1.66968,  -0.83886,   0.35522] ,
    [15,  -1.72454,  -0.87260,   0.37382] ,
    [16,  -1.77684,  -0.90499,   0.39251] ,
    [17,  -1.82668,  -0.93599,   0.41119] ,
    [18,  -1.87414,  -0.96560,   0.42973] ,    
    [19,  -1.91928,  -0.99382,   0.44804] ,
    [20,  -1.96220,  -1.02066,   0.46603] ,
    [21,  -2.00297,  -1.04614,   0.48363] ,
    [22,  -2.04169,  -1.07029,   0.50078] ,
    [23,  -2.07844,  -1.09313,   0.51742] ,
    [24,  -2.11332,  -1.11472,   0.53355] ,
    [25,  -2.14642,  -1.13511,   0.54912] ,
    [30,  -2.28839,  -1.22067,   0.61848] ,
    [40,  -2.48641,  -1.33024,   0.71925] ,
    [50,  -2.61710,  -1.39016,   0.78410] ,
    [60,  -2.71254,  -1.42377,   0.82683] ,
    [70,  -2.78865,  -1.44313,   0.85550] ,
    [80,  -2.85368,  -1.45474,   0.87479] ,
    [90,  -2.91216,  -1.46226,   0.88764] ,
    [100,  -2.96672,  -1.46787,   0.89598] ,
    [120,  -3.06983,  -1.47811,   0.90421] ,
    [140,  -3.16950,  -1.49082,   0.90634] ,
    [160,  -3.26809,  -1.50771,   0.90603] ,
    [180,  -3.36633,  -1.52909,   0.90532] ,
    [200,  -3.48436,  -1.55473,   0.90547] ,
    [250,  -3.70773,  -1.63448,   0.91388] ,
    [300,  -3.94607,  -1.73053,   0.93714] ,
    [350,  -4.17591,  -1.83593,   0.97495] ,
    [400,  -4.39433,  -1.94515,   1.02467] ,
    [500,  -4.78872,  -2.15940,   1.14615] ,
    [600,  -5.12008,  -2.35243,   1.27714] ,
    [800,  -5.59959,  -2.64798,   1.50995] ,
    [1000,  -5.88447,  -2.83157,  1.67325] ,
    [1500,  -6.15106,  -3.00957,   1.84797] ,
    [2000,  -6.20058,  -3.04408,   1.88423] ,
    [3000,  -6.21044,  -3.05176,   1.89114] ,
    [5000,  -6.21155,  -3.05324,   1.89118] ,
    [10000,  -6.21226,  -3.05427,   1.89110]])
    
    l  =  data[:,0]
    hl =  data[:,1]
    kl =  numpy.divide(data[:,2], l)
    ll =  numpy.divide(data[:,3], l)
    
    if frame == 'CM':
        hl[0] = hl[0] - 1
        ll[0] = ll[0] - 1
        kl[0] = kl[0] - 1
        print('Love numbers are in center of mass frame')
    elif frame == 'CF':
        hlo = hl[0]
        llo = ll[0]
        hl[0] = (hlo - llo) * 2/3
        ll[0] = (hlo - llo) * (-1/3)
        kl[0] = ((-2/3)*llo) - ((-1/3)*hlo)
        print('Love numbers are in center of figure frame')
    elif frame == 'CE':
        print('Love numbers are in center of solid Earth frame')
    else:
        lovenrPREM.exit('Please choose a compatible frame of reference: one of CM, CF, or CE')
    
    
    n = range(0, lmax+1, 1)
    kn = numpy.interp(n,l,kl)
    hn = numpy.interp(n,l,hl)
    ln = numpy.interp(n,l,ll)
    kn[0] = 0 
    hn[0] = 0
    ln[0] = 0
    return(kn,hn,ln)
    
