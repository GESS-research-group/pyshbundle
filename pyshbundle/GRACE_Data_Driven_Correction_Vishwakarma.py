#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Created on Tue Jun 28 14:17:50 2022
#@author: Amin Shakya, Interdisciplinary Center for Water Research (ICWaR), Indian Institute of Science (IISc)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Original Codes

#Data-driven approach for leakage and attenuation control
#based on Vishwakarma et. al., (2021)

#   Detailed explanation goes here

# Input: F, a cell matrix with one column containing SH coefficients
#      : cf, the column in F that contains SH coefficients from GRACE
#      : GaussianR, radius of the Gaussian filter (recommened = 400)
#      : basins, mask functions of basin, a cell data format with one
#      column and each entry is a 360 x 720 matrix with 1 inside the
#      catchment and 0 outside

# Output : every output has a size (number of months x basins)
#        : RecoveredTWS, corrected data-driven time-series (Least Squares fit method)
#        : RecoveredTWS2, corrected data-driven time-series (shift and amplify method)
#        : FilteredTS, gaussian filtered GRACE TWS time-series for all the basins. 
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

import numpy as np
import numpy.matlib as npm
import scipy as sc
import scipy.io
#CS2SC, gsha, gshs, gaussian

from pyshbundle import gaussian
from pyshbundle import cs2sc
from pyshbundle import gshs
from pyshbundle import gsha
from pyshbundle import naninterp
from pyshbundle import Phase_calc

def deg_to_rad(deg: float):
    """Converts angle from degree to radian

    Args:
        deg (float): Angle in degree

    Returns:
        float: Angle in Radian
    
    Todo:
        + Inbuilt function available in numpy module
    """
    return deg * np.pi/180

def GRACE_Data_Driven_Correction_Vishwakarma(F, cf, GaussianR, basins):
    """_summary_

    Args:
        F (_type_): a cell matrix with one column containing SH coefficients
        cf (_type_): the column in F that contains SH coefficients from GRACE
        GaussianR (_type_): radius of the Gaussian filter (recommened = 400)
        basins (_type_): mask functions of basin, a cell data format with one
                        column and each entry is a 360 x 720 matrix with 1 inside the
                        catchment and 0 outside

    Raises:
        Exception: corrected data-driven time-series (Least Squares fit method)
        Exception: corrected data-driven time-series (shift and amplify method)
        Exception: gaussian filtered GRACE TWS time-series for all the basins.

    Returns:
        _type_: _description_
    
    Todo:
        + TypeError
    """
    deg = 0.5
    deg_rad = deg_to_rad(deg)
    
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

    Weights = gaussian.gaussian(l-1, GaussianR) 
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
                Ft = cs2sc.cs2sc(f[m][0]).astype('longdouble') 
            else:
                Ft = f[m][0].astype('longdouble') 
                
           
            fFld__, _, _ = gshs.gshs(Ft * filter_, qty, 'cell', int(180/deg), 0, 0) 
            ffFld__, _, _ = gshs.gshs((Ft * filter_ * filter_), qty, 'cell', int(180/deg), 0, 0)
            
            if m == 0:
                fFld = np.zeros((r,fFld__.shape[0],fFld__.shape[1]), dtype = 'longdouble') 
                ffFld = np.zeros((r, ffFld__.shape[0], ffFld__.shape[1]), dtype = 'longdouble')
                
            fFld[m] = fFld__
            ffFld[m] = ffFld__
            
        long = 360/deg
        Area = Areahalfdeg
    else:
        raise Exception("enter CS coefficients")
        

        
    #Declaration of size of the vectors:
    cid = len(basins) #Here basins is a dictionary with each element storing nd array
    tsleaktotalf = np.zeros([r, cid], dtype = 'longdouble')
    tsleaktotalff = np.zeros([r, cid], dtype = 'longdouble')
    
    ftsleaktotal = np.zeros([r, cid], dtype = 'longdouble')
    fftsleaktotal = np.zeros([r, cid], dtype = 'longdouble')
    
    lhat = np.zeros([r, cid], dtype = 'longdouble')
    
    bfDevRegAv = np.zeros([r, cid], dtype = 'longdouble')
    bbfDevRegAv = np.zeros([r, cid], dtype = 'longdouble')

    FilteredTS = np.zeros([r, cid], dtype = 'longdouble')
    filfilts = np.zeros([r, cid], dtype = 'longdouble')

    leakage = np.zeros([r, cid], dtype = 'longdouble')
    leakager = np.zeros([r, cid], dtype = 'longdouble')   
    

    
    for rbasin in range(0, cid):
        #Get the basin functions ready
       
        #Basin functions, filtered basin function and transfer function Kappa
        Rb = basins[rbasin][0] 
        csRb = gsha.gsha(Rb, 'mean', 'block', long/2) 
        csF = cs2sc.cs2sc(csRb[0:l, 0:l]) 
        filRb_ = gshs(csF * filter_, 'none', 'cell', int(long/2), 0, 0) 
        filRb = filRb_[0]
        kappa = (1-Rb) * filRb
         
        
    
        fF = np.zeros((fFld__.shape[0],fFld__.shape[1]), dtype = 'longdouble')
        ffF = np.zeros((fFld__.shape[0],fFld__.shape[1]), dtype = 'longdouble')
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
        A[:,1] = naninterp.naninterp(tsleaktotalff[:, i])
        lssol_ = sc.linalg.lstsq(A, naninterp.naninterp(tsleaktotalf[:, i])) #returns a tuple of solution "x", residue and rank of matrix A; for A x = B
        lssol = lssol_[0]
        bl.append(lssol[2-1])
        #Working till here 2022-10-21 1530pm
    
    multp = npm.repmat(b, r, 1) 
    devint = bfDevRegAv * multp
    multp = npm.repmat(bl, r, 1)
    leakLS = tsleaktotalf * multp
    
    
    ps = Phase_calc.Phase_calc(tsleaktotalf,tsleaktotalff)
    
    
    
    #Compute the near true leakage
    
    for i in range(0, cid):   
        ftsleaktotal[:,i] = naninterp.naninterp(tsleaktotalf[:,i]) #Replaces gaps (NaN values) with an itnerpolated value in the leakage time series from once filtered fields
        fftsleaktotal[:,i] = naninterp.naninterp(tsleaktotalff[:,i]) #replace the gaps (NaN values) with an interpolated value in leakage time series from twice filtered fields
        
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
        
        
        
    
    
            