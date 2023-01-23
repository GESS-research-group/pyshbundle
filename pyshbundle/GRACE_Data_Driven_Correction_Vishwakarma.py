#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 14:17:50 2022
%% data-driven approach for leakage and attenuation control
% by Bramha Dutt Vishwakarma, University of Bristol, UK
% date: 18/11/2019
% the function uses data-driven approach published and explained in
% Vishwakarma et al. 2017, AGU WRR paper. doi:10.1002/2017WR021150.
%   Detailed explanation goes here

% Input: F, a cell matrix with one column containing SH coefficients
%      : cf, the column in F that contains SH coefficients from GRACE
%      : GaussianR, radius of the Gaussian filter (recommened = 400)
%      : basins, mask functions of basin, a cell data format with one
%      column and each entry is a 360 x 720 matrix with 1 inside the
%      catchment and 0 outside

% Output : every output has a size (number of months x basins)
%        : RecoveredTWS, corrected data-driven time-series (Least Squares fit method)
%        : RecoveredTWS2, corrected data-driven time-series (shift and amplify method)
%        : FilteredTS, gaussian filtered GRACE TWS time-series for all the basins. 

%% Further details are below:

%% create area matrix
@author: Amin Shakya
"""
import numpy as np
import numpy.matlib as npm
import scipy as sc
import scipy.io
#CS2SC, gsha, gshs, gaussian

from . import gaussian
from . import cs2sc
from . import gshs
from . import gsha
from . import naninterp
from . import Phase_calc


'''
GRACE_sh = scipy.io.loadmat('/home/bramha/Desktop/5hk/Papers/01_IISc/Bramha/Data_Nature/GRACE_SH_product.mat')
field1 = GRACE_sh["GRACE_SH"]
field = field1.T[2][0]
F = field1
#F = GRACE_sh
cf = 3
GaussianR = 400

example_data = scipy.io.loadmat('/home/bramha/Desktop/5hk/Papers/01_IISc/Bramha/Data_Nature/example_data.mat')
basins = example_data["mask_10_catchments"]
'''

    
'''
Debug on 2022-09-06
All data loaded till Matlab line 95 (fF and ffF)

GDC_run_line95_20220902 = scipy.io.loadmat('/home/bramha/Desktop/5hk/Papers/01_IISc/Bramha/Data_Nature/Downscaling_implementation_scripts/GDDC_run_line95_20220902.mat')

fF = GDC_run_line95_20220902["fF"]
ffF = GDC_run_line95_20220902["ffF"]
tsleaktotalf = GDC_run_line95_20220902["tsleaktotalf"]
tsleaktotalff = GDC_run_line95_20220902["tsleaktotalff"]
FilteredTS = GDC_run_line95_20220902["FilteredTS"]
filfilts = GDC_run_line95_20220902["filfilts"]
bfDevRegAv = GDC_run_line95_20220902["bfDevRegAv"]
bbfDevRegAv = GDC_run_line95_20220902["bbfDevRegAv"]
'''
    #example_data = scipy.io.loadmat('/home/bramha/Desktop/5hk/Papers/01_IISc/Bramha/Data_Nature/example_data.mat')
#    basins = example_data["mask_10_catchments"]
    
#    GaussianR = 400
#   cf = 3

#csRb is obtained from gsha; however gsha is not completed as of 2022-08-23; gsha completed as of 2022-08-26    
#    csRb_GDDC_gsha_20220823 = scipy.io.loadmat('/mnt/Data/5hk/Project_STC/Mat2Py/sample_data/GDDC_gsha/csRb_GDDC_gsha_20220823.mat')
#    csRb = csRb_GDDC_gsha_20220823["csRb"]

def GRACE_Data_Driven_Correction_Vishwakarma(F,cf,GaussianR, basins):
    def deg_to_rad(deg):
        return deg * np.pi/180
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

    Weights = gaussian(l-1, GaussianR) 
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
                Ft = cs2sc(f[m][0]).astype('longdouble') #In matlab, cs2sc is defined with 2 params; the second optional parameter being backval; not considered in the python script
            else:
                Ft = f[m][0].astype('longdouble') #In matlab F[m][cf] is the definition; double check conversion in case of error
                
        #Code checked till here on July 15; looks ok till here; need gshs to be converted before we can proceed
        #Code checked till here on Aug 26, 2022; for l == cfield and flag_cs == 0, the code is working well.
            
            fFld__, _, _ = gshs(Ft * filter_, qty, 'cell', int(180/deg), 0, 0) 
            ffFld__, _, _ = gshs((Ft * filter_ * filter_), qty, 'cell', int(180/deg), 0, 0)
            
            if m == 0:
                fFld = np.zeros((r,fFld__.shape[0],fFld__.shape[1]), dtype = 'longdouble') #Value not matching with matlab; datatype changed to longdouble 2022-10-17
                ffFld = np.zeros((r, ffFld__.shape[0], ffFld__.shape[1]), dtype = 'longdouble')
                # fFld = np.zeros((r), dtype = 'longdouble') #Value not matching with matlab; datatype changed to longdouble 2022-10-17
                # ffFld = np.zeros((r), dtype = 'longdouble')
            fFld[m] = fFld__
            ffFld[m] = ffFld__
            
        long = 360/deg
        Area = Areahalfdeg
    else:
        raise Exception("enter CS coefficients")
        
#        Checked till here on Aug 22, 2022 12.15pm
#Working till here 2022-10-19 16.00pm
        
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
    
    
    
    
    #gshs not working in this loop 2022-09-22
    # fF = np.zeros((r,fFld__.shape[0],fFld__.shape[1]), dtype = 'longdouble')
    # ffF = np.zeros((r,fFld__.shape[0],fFld__.shape[1]), dtype = 'longdouble')
    
    # fF = np.concatenate((fFld[:,:,int(fF.shape[2]/2):], fFld[:,:,:int(fF.shape[2]/2)]), axis=2)
    # ffF = np.concatenate((ffFld[:,:,int(ffF.shape[2]/2):], ffFld[:,:,:int(ffF.shape[2]/2)]), axis=2)
    
    for rbasin in range(0, cid):
        #Get the basin functions ready
    #Start timer to be added later, if required
    
    #Basin functions, filtered basin function and transfer function Kappa
        Rb = basins[rbasin][0] #looks good; double check later 2022-08-23
        csRb = gsha(Rb, 'mean', 'block', long/2) #csRb values not matching with csRb in matlab 2022-10-17
        csF = cs2sc(csRb[0:l, 0:l]) #seems to be working 2022-08-23
        filRb_ = gshs(csF * filter_, 'none', 'cell', int(long/2), 0, 0) #This line not working 2022-09-22 2022-10-19
        filRb = filRb_[0]
        kappa = (1-Rb) * filRb
    
    
        
        
    #Checked till here on 2022-08-23; everything above seems to work fine; csRb requires gsha; gsha not yet coded
    #Checked till here on 2022-08-29; looks to be working good
        #Please double check how the nan definitions are being enforced ~5hk 2022-08-30
        fF = np.zeros((fFld__.shape[0],fFld__.shape[1]), dtype = 'longdouble')
        ffF = np.zeros((fFld__.shape[0],fFld__.shape[1]), dtype = 'longdouble')
        for m in range(0,r):
#            fF  = fFld[m][0][:, 360:], fFld[m][0][:, 0:360]
#            ffF = ffFld[m][0][:, 360:], ffFld[m][0][:, 0:360]
            
            
            fF = np.concatenate((fFld[m,:,int(fF.shape[1]/2):], fFld[m,:,:int(fF.shape[1]/2)]), axis=1)
            ffF = np.concatenate((ffFld[m,:,int(ffF.shape[1]/2):], ffFld[m,:,:int(ffF.shape[1]/2)]), axis=1)
            #if False:    
            if np.isnan(fF[:20,:20]).any(): #if there is a gap in time series, fill it with NaNs; logic looks to be working
                                                #maybe the logic is not working afterall; should items be selectively nan? Pls doublecheck 2022-08-29
#                m, rbasin = np.nan
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
                
                
       #bfDevRegAv and bbfDevRegAv working as required 2022-10-20 12pm
       #Code working till here as required 2022-10-20 12pm
    
    
    b = list()
    bl = list()
    for i in range(0, cid):
        # A = np.ones([r,1]), naninterp(bbfDevRegAv[:, i]) #This line needs to be reconfigured to generate a 60*2 array 2022-09-06
        
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
    #Checked till here 2022-10-21 looks to be working well
    
    ps = Phase_calc(tsleaktotalf,tsleaktotalff)
    #Checked till here on 2022-10-20; phase_calc gives nearly same values between matlab and python
    #Checked till here on 2022-10-21; phase_calc gives nearly same values between matlab and python
    
    
    #Compute the near true leakage
    #Check here tomorrow 2022-10-20
    for i in range(0, cid):   
        ftsleaktotal[:,i] = naninterp(tsleaktotalf[:,i]) #Replaces gaps (NaN values) with an itnerpolated value in the leakage time series from once filtered fields
        fftsleaktotal[:,i] = naninterp(tsleaktotalff[:,i]) #replace the gaps (NaN values) with an interpolated value in leakage time series from twice filtered fields
        
        X = sc.fft.fft(ftsleaktotal[:,i]) #take fast Fourier transform #check shape of X 2022-10-21
        p = -ps[0,i] / r #compute the fraction of the time period by which the time series is to be shiftes
        Y = np.exp(1j * np.pi * p * ((np.arange(r)) - r/2) / r) #compute the Conjugate-Symmetric shift #Working 2022-10-21
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
        #Working till here 2022-10-21
    
    
    #compute the ratio between the amplitude of the shifted leakage from once filtered fields and leakage from twice filtered fields
    rfn=leakage/fftsleaktotal
    rfn[(rfn)>=2] = 1
    rfn[(rfn)<=-2] = -1
    rfn = np.sum(np.abs(rfn), axis = 0)
    rfn=rfn/r # amplitude ratio
    
    
    lhat = leakager * rfn #apply the amplitude ratio to the shifted leakage timeseries from the once filtered fields to get the near true leakage
    lhat[np.isnan(FilteredTS)] = np.nan #reintroduce nan for data gaps
    leakLS[np.isnan(FilteredTS)] = np.nan
    RecoveredTWS = FilteredTS - leakLS - devint
    RecoveredTWS2 = FilteredTS - lhat - devint
    
    return RecoveredTWS, RecoveredTWS2, FilteredTS
        
        
        
    
    
            