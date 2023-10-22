#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Created on Mon Aug 29 09:47:38 2022

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
   
# @author: Vivek Yadav, Interdisciplinary Center for Water Research (ICWaR), Indian Institute of Science (IISc)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

import os
#os.chdir(path_functions)
from pyshbundle import gaussian
from pyshbundle import gshs
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

def tws_cal(data, lmax: int, gs: float, r, m):
    """_summary_

    Args:
        data (np.ndarray): SC coefficients
        lmax (int): maximum degree
        gs (float): grid size
        r (_type_): _description_
        m (_type_): _description_
    """
    SC = data
    
    gfilter = gaussian.gaussian(lmax,r)
    grid_y = int(180/gs)
    grid_x = int(360/gs)
    tws_f = np.zeros([m,grid_y,grid_x], dtype ='longdouble')
    for i in tqdm(range(0,m,1)):
        field = SC[i,0:lmax+1,96-lmax:96+lmax+1]
        shfil = np.zeros([lmax+1,2*lmax+1])

        for j in range(0,2*lmax+1,1):
            shfil[:,j] = gfilter[:,0] * field[:,j]
        
        quant = 'water' 
        grd = 'cell'
        n = int(180/gs) 
        h = 0 
        jflag = 0
        
        ff = gshs.gshs(shfil, quant, grd, n, h, jflag)[0]
        
        ff = ff*1000
        tws_f[i,:,0:int(grid_x/2)] = ff[:,int(grid_x/2):]
        tws_f[i,:,int(grid_x/2):] = ff[:,0:int(grid_x/2)]   
    
    plt.imshow(tws_f[0,:,:])
    return(tws_f)

def apply_gaussian(sc_coeff, gaussian_coeff, lmax):
    
    # filtered SH Coeff
    shfil = np.zeros([lmax+1, 2 * lmax+1])

    # applying filter on substracted coeff
    for j in range(0, 2*lmax+1, 1):
        shfil[:,j] = gaussian_coeff[:,0] * sc_coeff[:,j]
    
    return shfil