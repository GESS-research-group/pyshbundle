#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# @author: Vivek Yadav, Interdisciplinary Center for Water Research (ICWaR), Indian Institute of Science (IISc)
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
# Original Codes

# IPLM Integrals of the fully normalized associated Legendre functions
# over blocks for a selected order M. 

#    HOW: p = iplm(l,m,theRAD)		- assumes dt = theRAD(2)-theRAD(1)
#        p = iplm(l,m,theRAD,dt)

#    IN:
#       l ........ degree (vector). Integer, but not necessarily monotonic.
#                For l < m a vector of zeros will be returned.
#       m ........ order (scalar)
#       theRAD ... co-latitude [rad] (vector)
#       dt ....... integration block-size [rad] (scalar). Default: dt = theRAD(2)-theRAD(1)

#    OUT: 
#       p ........ Matrix with integrated Legendre functions.
#                Functions are integrated from theRAD(i)-dt/2 till theRAD(i)+dt/2.
#                The matrix has length(TH) rows and length(L) columns, unless L 
#                or TH is scalar. Then the output vector follows the shape of 
#                respectively L or TH. 

#    USES:
#       plm

#    REMARKS:
#       The blocks at the pole might become too large under circumstances.
#       This is not treated separately, i.e. unwanted output may appear.
#       In case TH is scalar, dt will be 1 (arbitrarily).
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

import numpy as np
import sys
from pyshbundle import plm

def iplm(l, m:int, theRAD, dt=-9999):
    """IPLM Integrals of the fully normalized associated Legendre functions
        over blocks for a selected order M. 

    Args:
        l (_type_): degree (vector). Integer, but not necessarily monotonic.
                For l < m a vector of zeros will be returned.
        m (int): order (scalar)
        theRAD (np.array): co-latitude [rad] (vector)
        dt (int, optional): integration block-size [rad] (scalar). Defaults to -9999.
    
    Returns:
        np.ndarray: Matrix with integrated Legendre functions.
                Functions are integrated from theRAD(i)-dt/2 till theRAD(i)+dt/2.
                The matrix has length(TH) rows and length(L) columns, unless L 
                or TH is scalar. Then the output vector follows the shape of 
                respectively L or TH.
    
    Notes:
        The blocks at the pole might become too large under circumstances.
        This is not treated separately, i.e. unwanted output may appear.
        In case TH is scalar, dt will be 1 (arbitrarily).
    
    TO DO:
        Instead of using sys.exit() we could raise exceptions - that would be a better way of error handling
    
    Uses:
        `plm`
    """
    if dt == -9999:
        if len(theRAD) == 1:
            dt = np.pi/180
        else:
            dt = theRAD[1] - theRAD[0]
    if  min(l.shape) != 1:
        print('Degree l must be a vector (or scalar)') 
        sys.exit([])
    if  np.remainder(l,1).all() != 0:
        print('Vector l contains non-integers !!') 
        sys.exit([])
    # if 
    #     print('Order m must be a scalar') 
    #     sys.exit([])
    if  np.remainder(m,1) != 0:
        print('Order must be integer')
        sys.exit([])
    # if min(dt.shape) !=1:
    #     print('Block size DT must be scalar.') 
        sys.exit([])
    if dt == 0: 
        print('DT cannot be zero') 
        sys.exit([])
        
    # init
    lcol = len(l)
    trow = len(theRAD)
    n = len(theRAD)
    theRAD.T
    if min(theRAD) < 0 or max(theRAD) > np.pi:
        print('Is the co-latitude ''theta'' given in radian?')
        sys.exit([])
    lmax = max(l[0])
    mfix = m
    lvec = np.transpose(l) 
    l = np.arange(mfix,lmax+1,1)
    
    # Initialization of cosine, sine and Plm functions
    stplus  = np.sin(theRAD+dt/2)
    stmin   = np.sin(theRAD-dt/2)
    ctplus  = np.cos(theRAD+dt/2)
    ctmin   = np.cos(theRAD-dt/2)
    plmplus = np.ones([n,lmax+1])
    plmmin = np.ones([n,lmax + 1])
    plmplus[:,l] = plm.plm(np.array([l]),mfix,(theRAD + dt/2),3,1)[:,:,0]                  # Tesserals
    plmmin[:,l] = plm.plm(np.array([l]),mfix,(theRAD - dt/2),3,1)[:,:,0] 
    if mfix > 0:
        m = np.arange(1,mfix + 1,1)
        mm = 2*m
        fac = np.sqrt(2*np.cumprod((mm+1)/mm))
        mgr, stp = np.meshgrid(m, stplus)
        fgr, stm = np.meshgrid(fac, stmin)
        plmplus[:, m] = fgr * np.power(stp, mgr)
        plmmin[:, m] = fgr * np.power(stm, mgr)
    ptmp = np.zeros([n, lmax +2 ])
    ptmp00 = np.cos(theRAD - dt/2) - ctplus
    ptmp11 = np.sqrt(3)/2 * (dt - ctplus* stplus + ctmin* stmin)
    ptmp10 = np.sqrt(3)/2 * (np.power(stplus,2) - np.power(stmin,2))
    ptmp[:,0] = ptmp00
    
    # Compute first the integrals of order m == 0
    if mfix == 0:
        ptmp[:,1] = ptmp10
        for l in range(2,lmax+1,1):              #loop over the degree l 
            rootnm = np.sqrt( (2*l+1)*(2*l-1)/np.power(l,2))
            root1nm = np.sqrt( (2*l-1)*(2*l-3)/np.power(l-1,2))
            ptmp[:,l] = rootnm/(l+1)*(((l-2)*ptmp[:,l-2]/root1nm).T + np.power(stplus,2)*plmplus[:,l-1].T - np.power(stmin,2)*plmmin[:,l-1].T )
    else:
        # Compute the integrals of order m > 0

        # First we compute the diagonal element IPmm (lmax == mfix)

        ptmp[:,1] = ptmp11
        for l in range(2,mfix+1,1):
            # print(l)
            rootmm = np.sqrt( (2*l+1)/(2*l) )
            root1mm = np.sqrt( (2*l-1)/(2*l-2))
            if l == 2:
                root1mm = np.sqrt(3)
            
            ptmp[:,l] = rootmm/(l+1)*( l*root1mm*ptmp[:,l-2].T - (ctplus*plmplus[:,l].T -ctmin*plmmin[:,l].T)/rootmm )
    #the arbitrary element IPlm ( computed only when lmax > mfix)        
        if lmax > mfix:
            l = mfix + 1
        #first we do the element IPlm, for which l - m = 1    
            rootnm = np.sqrt( (2*l+1)*(2*l-1)/(l+mfix)/(l-mfix))
            ptmp[:,l] = rootnm/(l+1)*(np.power(stplus,2)*plmplus[:,l-1].T - np.power(stmin,2)*plmmin[:,l-1].T)
        #now we do the rest
            for l in range(mfix+2,lmax+1,1):              #loop over the degree l
                rootnm = np.sqrt( (2*l+1) * (2*l-1) / (l+mfix) / (l-mfix) )
                root1nm = np.sqrt( (2*l-1) * (2*l-3) / (l-1+mfix) / (l-1-mfix) )
                ptmp[:,l] = rootnm/(l+1)*( (l-2)*ptmp[:,l-2].T/root1nm + np.power(stplus,2)*plmplus[:,l-1].T -np.power(stmin,2)*plmmin[:,l-1].T)


# The integrated functions have been computed. What remains to be done, is to
# extract the proper columns from ptmp, corresponding to the vector lvec. 
# If l or theta is scalar the output matrix p reduces to a vector. It should
# have the shape of respectively theta or l in that case.

# p     = zeros(n, length(lvec))
    lind = np.argwhere(lvec<mfix)[:,0]      #index into l < m
    pcol = lvec + 1                         #index into columns of ptmp
    pcol[lind] = (lmax + 2)*np.ones([len(lind),1])   #Now l < m points to last col
    p = ptmp[:,pcol[:,0]-1]                 #proper column extraction 
    
    if max(lvec.shape) == 1 and min(np.array([theRAD]).shape) == 1 and trow == 1:
        p = p.T
    if max(np.array([theRAD]).shape) == 1 and min(lvec.shape) == 1 and lcol == 1:
        p = p.T
    return p