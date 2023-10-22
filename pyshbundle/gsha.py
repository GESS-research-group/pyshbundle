#!/usr/bin/env python3
#-*- coding: utf-8 -*-

#Created on Wed Aug 24 09:26:32 2022
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# SHBundle original Docs
# GSHA global spherical harmonic analysis

# IN:
#    f ....... global field of size (lmax+1)*2*lmax or lmax*2*lmax
#    method .. string argument, defining the analysis method:
#              - 'ls' ..... least squares
#              - 'wls' .... weighted least squares 
#              - 'aq' ..... approximate quadrature 
#              - 'fnm' .... first neumann method
#              - 'snm' .... second neumann method
#              - 'mean' ... block mean values (use of integrated Plm)
#    grid .... optional string argument, defining the grid:
#              - 'pole', 'mesh' ...... (default if lmax+1), equi-angular (lmax+1)*2*lmax, 
#                                      including poles and Greenwich meridian.
#              - 'block', 'cell' ..... (default if lmax), equi-angular block midpoints lmax*2lmax
#              - 'neumann', 'gauss' .. Gauss-Neumann grid (lmax+1)*2*lmax
#    lmax .... maximum degree of development

# OUT:
#    cs ...... Clm, Slm in |C\S| format

# USES:
#    plm, iplm, neumann, sc2cs 

# SEE ALSO:
#    GSHS

# REMARKS:
#    TBD - Zlm-functions option
#        - eigengrav, GRS80
#        - When 'pole' grid, m = 1 yields singular Plm-matrix!
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

# @author: Amin Shakya, Interdisciplinary Center for Water Research (ICWaR), Indian Institute of Science (IISc)

import numpy as np
from pyshbundle import neumann
from pyshbundle import plm
from scipy import sparse
from scipy import linalg
from pyshbundle import iplm
from pyshbundle import sc2cs

def gsha(f, method: str, grid: str = None, lmax: int = -9999):
    """ GSHA - Global Spherical Harmonic Analysis

    Args:
        f (np.ndarray): global field of size $(l_{max} + 1) * 2 * l_{max}$ or $l_{max} * 2 * l_{max}$
        method (str): method to be used
        grid (str, optional): choose between 'block' or 'cell'. Defaults to None.
        lmax (int, optional): maximum degree of development. Defaults to -9999.

    Returns:
        np.ndarray: Clm, Slm in |C\S| format

    Raises:
        ValueError: grid argument can only be 'block' or 'cell'
        ValueError: Grid type entered is not right
        TypeError: Invalid size of matrix F
        TypeError: GRID and METHOD must be strings
        ValueError: 2nd Neumann method ONLY on a ''neumann''/''gauss'' GRID'
        ValueError: Block mean method ONLY on a ''block''/''cell'' GRID
        ValueError: Maximum degree of development is higher than number of rows of input.
    
    Uses:
        `plm`, `neumann`, `iplm`, `sc2cs`
    """
    rows, cols = f.shape
    
    if cols == 2 * rows: #Check conditions
        if lmax == -9999:
            lmax = rows
        
        if grid == None:
            grid = 'block'
        
        if (grid != 'block') and (grid != 'cell'):
            raise ValueError("Your GRID variable should be either block or cell")
        
        n = rows
        dt = 180 / n
        theta = np.arange(dt/2, 180+(dt/4), dt)
        lam = np.arange(dt/2, 360+(dt/4), dt)
    
    elif cols == 2 * rows - 2:
        if lmax == -9999:
            lmax = rows - 1
        if grid == None:
            grid = 'pole'
        
        n = rows - 1
        dt = 180 / n
        
        if (grid == 'pole') or (grid == 'mesh'):                   
            theta = np.arange(0, 180+(dt/4), dt)
            lam = np.arange(0, 360+(dt/4) - dt, dt)
        elif (grid == 'neumann') or (grid == 'gauss'): 
        # gw, gx = neumann(n+1) #For some reason, grule does not work for even values
            gw, gx = neumann.neumann(n)
            theta = np.arccos(np.flipud(gx)) * 180 / np.pi
            lam = np.arange(0, 360+(dt/4)-dt, dt)
            
            if len(gw.shape) == 1:
                gw = gw.reshape(gw.shape[0],1)
    
            if len(gx.shape) == 1:
                gx = gx.reshape(gx.shape[0],1)
        else:
            raise ValueError("Grid type entered is not right")
    else:
        raise TypeError("Invalid size of matrix F")
    
    theRAD = theta * np.pi / 180
    # if len(theRAD.shape) == 1:
    # theRAD = theRAD.reshape(theRAD.shape[0],1)
    

    # further diagnostics

    if (type(grid) != str) or (type(method) != str):
        raise TypeError("GRID and METHOD must be strings.")
    
    if (method == 'snm') and ((grid != 'neumann') and (grid != 'gauss')):
        raise ValueError('2nd Neumann method ONLY on a ''neumann''/''gauss'' GRID')
    
    if (method == 'mean') and ((grid != 'block') and (grid != 'cell')):
        raise ValueError('Block mean method ONLY on a ''block''/''cell'' GRID')
        
    if lmax > n:
        raise ValueError('Maximum degree of development is higher than number of rows of input.')
        
    # Reshape variables as required
        
    if len(lam.shape) == 1:
        lam = lam.reshape(1,lam.shape[0])
        
    # Init
    
    L = n
    clm = np.zeros((L+1, L+1), dtype='longdouble')
    slm = np.zeros((L+1, L+1), dtype='longdouble')
    
    
    # First step of analysis

    m = np.arange(L+1).reshape(1,L+1)
    c = np.cos((lam.T @ m) * np.pi/180)
    s = np.sin((lam.T @ m) * np.pi/180)
    
    
    # % preserving the orthogonality (except for 'mean' case)
    # % we distinguish between 'block' and 'pole' type grids (in lambda)
    
    if (grid == 'block') or (grid == 'cell'):
        if method == 'mean':
            dl = dt
            c[:,0] = dl / 360
            m = np.arange(1, L+1)
            ms = 2 / m * np.sin(m * dl/2 * np.pi/180) / np.pi
            c[:,1:(L+1)+1] = c[:,1:(L+1)+1] * ms  
            s[:,1:(L+1)+1] = s[:,1:(L+1)+1] * ms

        else:
            c = c/L
            s = s/L
            c[:,0] = c[:,1]/2
            s[:,L] = s[:,L]/2
            c[:,L] = np.zeros(2*n)
            s[:,0] = np.zeros(2*n)
    else:
        c = c/L
        s = s/L
        c[:,[0, L]] = c[:,[0, L]]/2	
        s[:,[0, L]] = np.zeros((2*n,2))	  
    
    
    a = f @ c
    b = f @ s    
    
    # Second step of analysis: Clm and Slm

    if method == 'ls':
        for m in range(L+1):
#            l = np.arange(m,L+1)
            l = np.arange(m,L+1).reshape(L+1-m, 1)
            l = l.T
            
            p = plm(l,m,theRAD, 3, 1)
            p = p[:,:,0]
            ai = a[:, m]
            bi = b[:, m]
            
            clm[m+1:L+2, m+1] = linalg.lstsq(p, ai)
            slm[m+1:L+2, m+1] = linalg.lstsq(p, bi)
            
            
    elif method == 'aq': #Approximate Quadrature
        si = np.sin(theRAD)
        si = 2 * si / np.sum(si)
        
        for m in range(L+1):
            l = np.arange(m, L+1).reshape(L+1-m, 1)
            l = l.T
            
            p = plm(l,m,theRAD, 3, 1)
            
            ai = a[:, m]
            bi = b[:, m]
                        
            clm[m:L+1, m] = (1 + (m == 0))/ 4 * p.T @ (si * ai)
            slm[m:L+1, m] = (1 + (m == 0))/ 4 * p.T @ (si * bi)
    
    elif method == 'fnm': #1st Neumann method (exact upto L/2)
        w = neumann.neumann(np.cos(theRAD))
        
        for m in range(L+1):
            l = np.arange(m, L+1).reshape(L+1-m, 1)
            l = l.T
            
            p = plm(l,m,theRAD, 3, 1)
            
            ai = a[:, m]
            bi = b[:, m]
                        
            clm[m:L+1, m] = (1 + (m == 0))/ 4 * p.T @ (w * ai)
            slm[m:L+1, m] = (1 + (m == 0))/ 4 * p.T @ (w * bi)
        
    elif method == 'snm': #2nd Neumann method (exact)
        for m in range(L+1):
            l = np.arange(m, L+1).reshape(L+1-m, 1)
            l = l.T
            
            p = plm(l,m,theRAD, 3, 1)
            
            ai = a[:, m]
            bi = b[:, m]
            
            clm[m:L+1, m] = (1 + (m == 0))/ 4 * p.T @ (gw * ai)
            slm[m:L+1, m] = (1 + (m == 0))/ 4 * p.T @ (gw * bi)
            
    elif method == 'mean':
        for m in range(L+1):
            print(m)
            #l = np.arange(m,L+1).reshape(L+1-m,1)
            #l = l.T
            
            
            l = np.array([np.arange(m,L+1, 1)])
        # l = np.array([[m]])
            
            p = iplm(l,m,theRAD)
        # p = p[:,-1]
            ai = a[:, m]
            bi = b[:, m]
            
            clm[m:L+1, m] = (1 + (m == 0))/ 4 * p.T @ ai
            slm[m:L+1, m] = (1 + (m == 0))/ 4 * p.T @ bi
                
    # Write the coefficients Clm & Slm in |C\S| format

    slm = np.fliplr(slm)
    cs = sc2cs(np.concatenate((slm[:, np.arange(L)], clm), axis = 1))
    cs = cs[:int(lmax+1), :int(lmax+1)]
    
    
    # np.save('/path/csRb.npy',cs)
    
    return cs
    