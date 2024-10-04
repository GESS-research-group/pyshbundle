# Core functionalities of PySHBundle
# Curator: Abhishek Mhamane, Vivek Kumar Yadav

# - - - - - - - - - - - - - - 
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
import numpy.matlib as npm
import scipy as sc
from scipy import signal
from scipy import linalg

from os import chdir, getcwd

from pyshbundle.reshape_SH_coefficients import cs2sc, sc2cs
from pyshbundle.shutils import normalklm, plm, iplm, eigengrav, ispec, Gaussian, naninterp, neumann


def gshs(field, quant = 'none', grd = 'mesh', n = -9999, h = 0, jflag = 1):
    """
    Global Spherical Harmonic Synthesis

    Args:
        field (numpy.ndarray): Matrix of SH coefficients, either in SC-triangle or CS-square format.
        quant (str, optional): Defining the field quantity. Defaults to 'none'.
        grd (str, optional): Defining the grid. Defaults to 'mesh'.
        n (int, optional): Defaults to -9999.
        h (int, optional): Defaults to 0.
        jflag (int, optional): Defaults to 1.

    Returns:
        (numpy.ndarray): The global Spherical Harmonics field.
        (numpy.ndarray): Vector of co-latitudes in radians.
        (numpy.ndarray): Vector of longitudes in radians.

    Raises:
        Exception: If the format of the field is incorrect.
        Exception: If n is not scalar.
        Exception: If n is not an integer.
        Exception: If the grid argument is not a string.

    Uses:
        `cs2sc`, `normalklm`, `plm`, `eigengrav`, `ispec`

    Todo:
        * Change general exceptions to specific and descriptive built-in ones.
        * Using the not and then check is not always advisable.
        * Check how to document valid options.

    Author:
        Amin Shakya, Interdisciplinary Center for Water Research (ICWaR), Indian Institute of Science (IISc)
        Vivek Kumar Yadav, Interdisciplinary Center for Water Research (ICWaR), Indian Institute of Science (IISc)
    """

    wd = getcwd()
    chdir(wd)
    
    rows, cols = field.shape
    
    if rows == cols:                    #field in CS-format 
        lmax = rows - 1
        field = cs2sc(field)
    elif cols - 2 * rows == -1:         #field in SC-format already
        lmax = rows - 1
    else:
        raise Exception("Check format of the field")
    
    if n == -9999:                      #(no value of n is input) -> set n = lmax
        n = lmax
    
    if not np.isscalar(n):
        raise Exception("n must be scalar")
    
    if not np.issubdtype(type(n), np.integer):
        raise Exception("n must be integer")
    
    if not type(grd) == str:
        raise Exception("Grid argument must be string")
        
    grd = grd.lower()
    
    
    
    #Grid Definition
    dt = np.pi/n
    
    if grd == 'pole' or grd == 'mesh':
        theRAD = np.arange(0, np.pi+dt*0.5, dt, dtype='float')
        lamRAD = np.arange(0, 2*np.pi, dt, dtype='float')
    elif grd == 'block' or grd == 'cell':
        theRAD = np.arange(dt/2, np.pi + dt*0.5, dt, dtype='float')
        lamRAD = np.arange(dt/2, 2*np.pi + dt*0.5, dt, dtype='float')
    else:
        raise Exception("Incorrect grid type input")
    
    nlat = len(theRAD)
    nlon = len(lamRAD)
    
    

#   -------------------------------------------------------------------------
#   Preprocessing on the coefficients:
        # - subtract reference field (if jflag is set)
        # - specific transfer
        # - upward continuation
#   -------------------------------------------------------------------------
    
    if jflag:
        field = field - cs2sc(normalklm(lmax+1))
        
    l = np.arange(0, lmax+1)
    transf = np.array([eigengrav(lmax, quant, h)])[0, :, :].T
    
    field = field * np.matmul(transf, np.ones((1, 2*lmax+1)), dtype='float')
    
    

# -------------------------------------------------------------------------
        # Size declarations and start the waitbar:
        # Note that the definition of dlam causes straight zero-padding in case N > L.
        # When N < L, there will be zero-padding up to the smallest integer multiple
        # of N larger than L. After the Fourier transformation (see below), the
        # proper samples have to be picked out, with stepsize dlam.
# -------------------------------------------------------------------------

    
    dlam = int(np.ceil(lmax/n))             #longitude step size
    abcols = dlam*n + 1                     #columns required in A and B
    a = np.zeros((nlat, int(abcols)), dtype='float')
    b = np.zeros((nlat, int(abcols)), dtype='float')
     


    m = 0
    c = field[m:lmax+1, lmax+m] 
    l = np.array([np.arange(m,lmax+1)])
    p = plm(l, m, theRAD, nargin = 3, nargout = 1)[:,:,0]
    a[:, m] = np.dot(p,c) 
    b[:, m] = np.zeros(nlat) 
    
    
    
    for m in range(1,lmax+1,1):
        c = field[m:lmax+1,lmax+m]
        s = field[m:lmax+1,lmax-m]
        
        l = np.array([np.arange(m,lmax+1)])
        p = plm(l, m, theRAD, nargin = 3, nargout = 1)[:,:,0]
        a[:, m] = np.dot(p,c)
        b[:, m] = np.dot(p,s)
        
    del field


#   -------------------------------------------------------------------------
        # The second synthesis step consists of an inverse Fourier transformation
        # over the rows of a and b. 
        # In case of 'block', the spectrum has to be shifted first.
        # When no zero-padding has been applied, the last b-coefficients must be set to
        # zero. Otherwise the number of longitude samples will be 2N+1 instead of 2N.
        # For N=L this corresponds to setting SLL=0!
#  -------------------------------------------------------------------------


    if grd =='block' or grd == 'cell': 
      m      = np.arange(0,abcols,1)
      cshift = np.array([np.ones(nlat)], dtype='float').T * np.array([np.cos(m*np.pi/2/n)], dtype='float');	# cshift/sshift describe the 
      sshift = np.array([np.ones(nlat)], dtype='float').T * np.array([np.sin(m*np.pi/2/n)], dtype='float');	# half-blocksize lambda shift.
      atemp  =  cshift*a + sshift*b
      b      = -sshift*a + cshift*b
      a      = atemp
    
    
    
    if np.remainder(n,lmax) == 0:               #Case without zero-padding
        b[:,abcols-1] = np.zeros(nlat)
    
    #Code for ispec
    
    f = ispec(a.T, b.T).T
    if dlam > 1: 
        f = f[:,np.arange(1,dlam*nlon+1,dlam)]

    return f, theRAD, lamRAD
    


def gsha(f, method: str, grid: str = None, lmax: int = -9999):
    """
    Global Spherical Harmonic Analysis, inverse of GSHS.

    Args:
        f (numpy.ndarray): Global field of size $(l_{max} + 1) * 2 * l_{max}$ or $l_{max} * 2 * l_{max}$.
        method (str): Method to be used.
        grid (str, optional): Choose between 'block' or 'cell'. Defaults to None.
        lmax (int, optional): Maximum degree of development. Defaults to -9999.

    Returns:
        (numpy.ndarray): Spherical harmonics coefficients Clm, Slm in |C\S| format.

    Raises:
        ValueError: If grid argument can only be 'block' or 'cell'.
        ValueError: If grid type entered is not right.
        TypeError: If invalid size of matrix F.
        TypeError: If GRID and METHOD are not strings.
        ValueError: If 2nd Neumann method ONLY on a 'neumann'/'gauss' GRID.
        ValueError: If Block mean method ONLY on a 'block'/'cell' GRID.
        ValueError: If maximum degree of development is higher than number of rows of input.

    Remarks:
        TBD - Zlm-functions option
            - eigengrav, GRS80
            - When 'pole' grid, m = 1 yields singular Plm-matrix!

    Uses:
        `plm`, `neumann`, `iplm`, `sc2cs`

    Author:
        Amin Shakya, Interdisciplinary Center for Water Research (ICWaR), Indian Institute of Science (IISc)
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
            gw, gx = neumann(n)
            theta = np.arccos(np.flipud(gx)) * 180 / np.pi
            lam = np.arange(0, 360+(dt/4)-dt, dt)
            
            if len(gw.shape) == 1:
                gw = gw.reshape(gw.shape[0],1)
    
            if len(gx.shape) == 1:
                gx = gx.reshape(gx.shape[0],1)
        else:
            raise ValueError("Grid type entered is incorrect")
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
    clm = np.zeros((L+1, L+1), dtype='float')
    slm = np.zeros((L+1, L+1), dtype='float')
    
    
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
        w = neumann(np.cos(theRAD))
        
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
    
    
    return cs


def PhaseCalc(fts, ffts):
    """
    Calculates the phase difference between two time series based on the
    Hilbert transform method explained by Phillip et al.

    Args:
        fts (numpy.ndarray): Time-series 1.
        ffts (numpy.ndarray): Time-series 2.

    Returns:
        (numpy.ndarray): Phase difference between the two time series.
    
    References:
        Phillips, T., R. S. Nerem, B. Fox-Kemper, J. S. Famiglietti, and B. Rajagopalan (2012),
        The influence of ENSO on global terrestrial water storage using GRACE, Geophysical
        Research Letters, 39 (16), L16,705, doi:10.1029/2012GL052495.
    
    Author:
        Amin Shakya, Interdisciplinary Center for Water Research (ICWaR), Indian Institute of Science (IISc)
    """
    c = fts.shape[1]
    
    ps = np.zeros((1, c))
    
    filter_ = ~np.isnan(fts)
    filter__ = ~np.isnan(ffts)
    
    fts_ = fts[filter_] #Extract values and leave Nan
    ffts_ = ffts[filter__] #Extract values and leave Nan
    
    fts = fts_.reshape(int(fts_.shape[0]/c),c)
    ffts = ffts_.reshape(int(ffts_.shape[0]/c),c)
    
    rn = fts.shape[0]
    
    for i in range(c):
        # A = np.concatenate(np.ones((rn,1)), np.real(signal.hilbert(ffts[:, i])), np.imag(signal.hilbert(ffts[:, i]))) #design matrix
        
        A = np.array((np.ones((rn)), np.real(signal.hilbert(ffts[:, i])), np.imag(signal.hilbert(ffts[:, i])))).T
        
        A = A.astype('double')
        B = fts[:,i]
        B = B.astype('double')
        abc = np.linalg.lstsq(A.T @ A, A.T @ B)[0]
        
        ps[0,i] = np.arctan2(abc[3-1],abc[2-1])*(180/np.pi) #check indices and degree/radian
    return ps


# ----------------------------- GDDC --------------------------------------------------#

def deg_to_rad(deg: float):
    """Converts angle from degree to radian.

    Args:
        deg (float): Angle in degree.

    Returns:
        float: Angle in radian.

    Todo:
        - Inbuilt function available in numpy module.
    """
    return deg * np.pi/180

def GRACE_Data_Driven_Correction_Vishwakarma(F, cf, GaussianR, basins):
    """
    Signal leakage correction using data-driven methods.

    When GRACE data is applied for hydrological studies, the signal leakage is a common
    problem. This function uses data-driven methods to correct signal leakage in GRACE data.
    Please refer to paper (1) above for more details.

    Args:
        F (numpy.ndarray): A cell matrix with one column containing SH coefficients.
        cf (int): The column in F that contains SH coefficients from GRACE.
        GaussianR (float): Radius of the Gaussian filter (recommended = 400).
        basins (numpy.ndarray): Mask functions of basin, a cell data format with one
            column and each entry is a 360 x 720 matrix with 1 inside the
            catchment and 0 outside.

    Returns:
        tuple: A tuple containing:
            - RecoveredTWS (numpy.ndarray): Corrected data-driven time-series (Least Squares fit method).
            - RecoveredTWS2 (numpy.ndarray): Corrected data-driven time-series (shift and amplify method).
            - FilteredTS (numpy.ndarray): Gaussian filtered GRACE TWS time-series for all the basins.

    Raises:
        Exception: If there is an error in the corrected data-driven time-series (Least Squares fit method).
        Exception: If there is an error in the corrected data-driven time-series (shift and amplify method).
        Exception: If there is an error in the Gaussian filtered GRACE TWS time-series for all the basins.

    Author:
        Amin Shakya, Interdisciplinary Center for Water Research (ICWaR), Indian Institute of Science (IISc)
    """
    deg = 0.5
    x = np.linspace(0, 360-deg, int(360/deg))
    y = np.linspace(0, 180-deg, int(180/deg))
    x1 = np.linspace(deg, 360, int(360/deg))
    y1 = np.linspace(deg, 180, int(180/deg))
    lambdd,theta = np.meshgrid(x,y)  
    lambdd1,theta1 = np.meshgrid(x1,y1)
    
    theta_rad = deg_to_rad(theta)
    theta1_rad = deg_to_rad(theta1)
    
    #Areahalfdeg = (6378.137**2)*np.power(10,6)*np.pi/180*(np.multiply(a,b)) #Area matrix
    Areahalfdeg = (6378.137**2)*(((np.pi/180)*lambdd1) - ((np.pi/180)*lambdd))*(np.sin((np.pi/2) - theta_rad) - np.sin((np.pi/2) - theta1_rad))
    
    qty = 'water'
    
    if type(F) != np.ndarray:
        raise Exception("input GRACE field should be in Numpy Ndarray format, please check guidelines")
        
    
    if type(basins) != np.ndarray:
        raise Exception("input basin field should be in Numpy NdArray format, please check guidelines")
        

    r = F.shape[0] #No of entries in F numpy ndarrray
    
    cid = 1 #number of river catchments
    
    f = F[:,cf-1:cf]
    l = f[0][0].shape[0]
    cfield = f[0][0].shape[1]
    if cfield == l:
        flag_cs = 0
    else:
        flag_cs = 1

    Weights = Gaussian(l-1, GaussianR) 
    #gaussian returns weights as a list #gaussian is np.array()
    
    try: #Broadcase Weights into dimensions
        filter_ = np.ones([1,(2*(l-1))+1]) * Weights
    except:
        w0 = Weights.shape[0]
        Weights = Weights.reshape(w0,1)
        filter_ = np.ones([1,(2*(l-1))+1]) * Weights
    
    
    #SH Synthesis
    if l == cfield:
        for m in range(r):
            if flag_cs == 0:
                Ft = cs2sc(f[m][0]).astype('double') 
            else:
                Ft = f[m][0].astype('double') 
                
           
            fFld__, _, _ = gshs(Ft * filter_, qty, 'cell', int(180/deg), 0, 0) 
            ffFld__, _, _ = gshs((Ft * filter_ * filter_), qty, 'cell', int(180/deg), 0, 0)
            
            if m == 0:
                fFld = np.zeros((r,fFld__.shape[0],fFld__.shape[1]), dtype='double') 
                ffFld = np.zeros((r, ffFld__.shape[0], ffFld__.shape[1]), dtype='double')
                
            fFld[m] = fFld__
            ffFld[m] = ffFld__
            
        long = 360/deg
        Area = Areahalfdeg
    else:
        raise Exception("enter CS coefficients")
        

        
    #Declaration of size of the vectors:
    cid = len(basins) #Here basins is a dictionary with each element storing nd array
    tsleaktotalf = np.zeros([r, cid], dtype='double')
    tsleaktotalff = np.zeros([r, cid], dtype='double')
    
    ftsleaktotal = np.zeros([r, cid], dtype='double')
    fftsleaktotal = np.zeros([r, cid], dtype='double')
    
    lhat = np.zeros([r, cid], dtype='double')
    
    bfDevRegAv = np.zeros([r, cid], dtype='double')
    bbfDevRegAv = np.zeros([r, cid], dtype='double')

    FilteredTS = np.zeros([r, cid], dtype='double')
    filfilts = np.zeros([r, cid], dtype='double')

    leakage = np.zeros([r, cid], dtype='double')
    leakager = np.zeros([r, cid], dtype='double')   
    

    
    for rbasin in range(0, cid):
        #Get the basin functions ready
       
        #Basin functions, filtered basin function and transfer function Kappa
        Rb = basins[rbasin][0] 
        csRb = gsha(Rb, 'mean', 'block', long/2) 
        csF = cs2sc(csRb[0:l, 0:l]) 
        filRb_ = gshs(csF * filter_, 'none', 'cell', int(long/2), 0, 0) 
        filRb = filRb_[0]
        kappa = (1-Rb) * filRb
         
        
    
        fF = np.zeros((fFld__.shape[0],fFld__.shape[1]), dtype='double')
        ffF = np.zeros((fFld__.shape[0],fFld__.shape[1]), dtype='double')
        for m in range(0,r):
            
            
            fF = np.concatenate((fFld[m,:,int(fF.shape[1]/2):], fFld[m,:,:int(fF.shape[1]/2)]), axis=1)
            ffF = np.concatenate((ffFld[m,:,int(ffF.shape[1]/2):], ffFld[m,:,:int(ffF.shape[1]/2)]), axis=1)
            #if False:    
            if np.isnan(fF[:20,:20]).any(): #if there is a gap in time series, fill it with NaNs
                                                

                tsleaktotalf[m][rbasin] = np.nan
                tsleaktotalff[m][rbasin] = np.nan
                FilteredTS[m][rbasin] = np.nan
                filfilts[m][0:rbasin] = np.nan
                bfDevRegAv[m][rbasin] = np.nan
                bbfDevRegAv[m][0:rbasin] = np.nan
            
            else:
                #leakage time series from filtered and twice filtered fields
                tsleaktotalf[m][rbasin] = np.sum(fF * kappa * Area) / np.sum(Rb * Area)
                tsleaktotalff[m][rbasin] = np.sum(ffF * kappa * Area) / np.sum(Rb * Area)
                
                #time series from filtered fields
                FilteredTS[m][rbasin] = np.sum(fF * Rb * Area) / np.sum(Rb * Area)
                filfilts[m][rbasin] = np.sum(ffF * Rb * Area) / np.sum(Rb * Area)
                
                #Deviation integral timeseries
                bfDevRegAv[m][rbasin] = np.sum((fF * Rb - FilteredTS[m][rbasin]) * filRb * Area) / np.sum(Rb * Area) #working 2022-10-20
                bbfDevRegAv[m][rbasin] = np.sum((ffF * Rb - filfilts[m][rbasin]) * filRb * Area) / np.sum(Rb * Area)
                print(m)
                
                
       
    
    
    b = list()
    bl = list()
    for i in range(0, cid):
                
        A = np.ones([60,2])
        A[:,1] = naninterp(bbfDevRegAv[:, i]) #Pchip interpolate should contain atleast two elements
        
        lssol_ = sc.linalg.lstsq(A, naninterp(bfDevRegAv[:, i])) #returns a tuple of solution "x", residue and rank of matrix A; for A x = B
        lssol = lssol_[0] 
        
        b.append(lssol[2-1])
        
        
        A = np.ones([60,2])
        A[:,1] = naninterp(tsleaktotalff[:, i])
        lssol_ = sc.linalg.lstsq(A, naninterp(tsleaktotalf[:, i])) #returns a tuple of solution "x", residue and rank of matrix A; for A x = B
        lssol = lssol_[0]
        bl.append(lssol[2-1])
        #Working till here 2022-10-21 1530pm
    
    multp = npm.repmat(b, r, 1) 
    devint = bfDevRegAv * multp
    multp = npm.repmat(bl, r, 1)
    leakLS = tsleaktotalf * multp
    
    
    ps = PhaseCalc(tsleaktotalf,tsleaktotalff)
    
    
    
    #Compute the near true leakage
    
    for i in range(0, cid):   
        ftsleaktotal[:,i] = naninterp(tsleaktotalf[:,i]) #Replaces gaps (NaN values) with an itnerpolated value in the leakage time series from once filtered fields
        fftsleaktotal[:,i] = naninterp(tsleaktotalff[:,i]) #replace the gaps (NaN values) with an interpolated value in leakage time series from twice filtered fields
        
        X = sc.fft.fft(ftsleaktotal[:,i]) #take fast Fourier transform #check shape of X 2022-10-21
        p = -ps[0,i] / r #compute the fraction of the time period by which the time series is to be shiftes
        Y = np.exp(1j * np.pi * p * ((np.arange(r)) - r/2) / r) #compute the Conjugate-Symmetric shift 
        Z = X * Y #Apply the shift
        
        a = sc.fft.ifft(Z) #apply inverse fft
        
        con = np.conj(a)
        
        s = a + con
        
        z = s/2
        
        leakage[:,i] = z #shifted timeseries
        
        
        #Shift timeseriecs from once filtered fields in the direction of the time series from twice filtered fields, to later compute the amplitude ratio
        p = ps[0,i] / r #Fraction of a time period to shift data
        Y = np.exp(1j * np.pi * p * ((np.arange(r)) - r/2) / r) #compute the Conjugate-Symmetric shift
        Z = X * Y
        
        a = sc.fft.ifft(Z) #apply inverse fft
        
        con = np.conj(a)
        
        s = a + con
        
        z = s/2
        
        leakager[:,i] = z #shifted timeseries
        
    
    
    #compute the ratio between the amplitude of the shifted leakage from once filtered fields and leakage from twice filtered fields
    rfn = leakage/fftsleaktotal
    rfn[(rfn) >= 2] = 1
    rfn[(rfn) <= -2] = -1
    rfn = np.sum(np.abs(rfn), axis = 0)
    rfn=rfn/r # amplitude ratio
    
    
    lhat = leakager * rfn #apply the amplitude ratio to the shifted leakage timeseries from the once filtered fields to get the near true leakage
    lhat[np.isnan(FilteredTS)] = np.nan #reintroduce nan for data gaps
    leakLS[np.isnan(FilteredTS)] = np.nan
    RecoveredTWS = FilteredTS - leakLS - devint
    RecoveredTWS2 = FilteredTS - lhat - devint
    
    return RecoveredTWS, RecoveredTWS2, FilteredTS
