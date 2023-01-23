#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 09:26:32 2022

% GSHA global spherical harmonic analysis
%
% IN:
%    f ....... global field of size (lmax+1)*2*lmax or lmax*2*lmax
%    method .. string argument, defining the analysis method:
%              - 'ls' ..... least squares
%              - 'wls' .... weighted least squares 
%              - 'aq' ..... approximate quadrature 
%              - 'fnm' .... first neumann method
%              - 'snm' .... second neumann method
%              - 'mean' ... block mean values (use of integrated Plm)
%    grid .... optional string argument, defining the grid:
%              - 'pole', 'mesh' ...... (default if lmax+1), equi-angular (lmax+1)*2*lmax, 
%                                      including poles and Greenwich meridian.
%              - 'block', 'cell' ..... (default if lmax), equi-angular block midpoints lmax*2lmax
%              - 'neumann', 'gauss' .. Gauss-Neumann grid (lmax+1)*2*lmax
%    lmax .... maximum degree of development
%
% OUT:
%    cs ...... Clm, Slm in |C\S| format
%
% USES:
%    plm, iplm, neumann, sc2cs 
%
% SEE ALSO:
%    GSHS
%
% REMARKS:
%    TBD - Zlm-functions option
%        - eigengrav, GRS80
%        - When 'pole' grid, m = 1 yields singular Plm-matrix!

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Dimitris TSOULIS (DT), IAPG, TU-Munich
%    Nico SNEEUW (NS), IAGP, TU Munich
%    Markus ANTONI (MA), GI, Uni Stuttgart 
%    Matthias ROTH (MR), GI, Uni Stuttgart 
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    2014-01-15: MR, revise help text, beautify code
%    2012-01-23: MA, input of plm in radian
%    1999-02-01: NS, brush up (help text, layout, removal of unused commands, ...)
%                    restructuring of 'mean' method over 1st and 2nd analysis
%                    'mean' method as quadrature instead of LS in 2nd step 
%    1998-11-??: DT, initial version
% -------------------------------------------------------------------------
% license:
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the  Free  Software  Foundation; either version 3 of the License, or
%    (at your option) any later version.
%  
%    This  program is distributed in the hope that it will be useful, but 
%    WITHOUT   ANY   WARRANTY;  without  even  the  implied  warranty  of 
%    MERCHANTABILITY  or  FITNESS  FOR  A  PARTICULAR  PURPOSE.  See  the
%    GNU General Public License for more details.
%  
%    You  should  have  received a copy of the GNU General Public License
%    along with Octave; see the file COPYING.  
%    If not, see <http://www.gnu.org/licenses/>.
% -------------------------------------------------------------------------

%--------------------------------------------------------------------------
% diagnostics and preliminaries
%--------------------------------------------------------------------------
%narginchk(2, 4) % error(nargchk(2, 4, nargin))

%--------------------------------------------------------------------------
% Grid definition
%--------------------------------------------------------------------------

@author: Amin Shakya
"""
import numpy as np
from . import neumann
from . import plm
from scipy import sparse
from scipy import linalg
from . import iplm
from . import sc2cs
'''
Use the other code for debug for now
##SC = np.load("SC_numpy.npy")
#f1 = scipy.io.loadmat('/home/bramha/Desktop/5hk/Papers/01_IISc/Bramha/Data_Nature/f_gsha_20220825.mat')
##f = SC[0]
#f = f1['f'][0]
'''

'''
import scipy.io
f1 = scipy.io.loadmat('/home/bramha/Desktop/5hk/Papers/01_IISc/Bramha/Data_Nature/Downscaling_implementation_scripts/Rb_GDDC_gsha_20220826.mat', mat_dtype= 1)
f = f1['Rb']
method = 'mean'
grid = 'block'
lmax = 720/2
'''

def gsha(f, method, grid = None, lmax = -9999):
    rows, cols = f.shape
    
    if cols == 2 * rows: #Check conditions
        if lmax == -9999:
            lmax = rows
        
        if grid == None:
            grid = 'block'
        
        if (grid != 'block') and (grid != 'cell'):
            raise Exception("Your GRID variable should be either block or cell")
        
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
            lam = np.arange(0, 360+(dt/4)-dt, dt)
        elif (grid == 'neumann') or (grid == 'gauss'): 
#            gw, gx = neumann(n+1) #For some reason, grule does not work for even values
            gw, gx = neumann(n)
            theta = np.arccos(np.flipud(gx)) * 180 / np.pi
            lam = np.arange(0, 360+(dt/4)-dt, dt)
            
            if len(gw.shape) == 1:
                gw = gw.reshape(gw.shape[0],1)
    
            if len(gx.shape) == 1:
                gx = gx.reshape(gx.shape[0],1)
        else:
            raise Exception("Grid type entered is not right")
    else:
        raise Exception("Invalid size of matrix F")
    
    theRAD = theta * np.pi / 180
#    if len(theRAD.shape) == 1:
#        theRAD = theRAD.reshape(theRAD.shape[0],1)
    
    '''
    %--------------------------------------------------------------------------
% further diagnostics
%--------------------------------------------------------------------------
    '''
    
    if (type(grid) != str) or (type(method) != str):
        raise Exception("GRID and METHOD must be strings.")
    
    if (method == 'snm') and ((grid != 'neumann') and (grid != 'gauss')):
        raise Exception('2nd Neumann method ONLY on a ''neumann''/''gauss'' GRID')
    
    if (method == 'mean') and ((grid != 'block') and (grid != 'cell')):
        raise Exception('Block mean method ONLY on a ''block''/''cell'' GRID')
        
    if lmax > n:
        raise Exception('Maximum degree of development is higher than number of rows of input.')
        
    '''
    Reshape variables as required
    '''
    
    
    
    
    
#    if len(theta.shape) == 1:
#        theta = theta.reshape(theta.shape[0],1)
        
    if len(lam.shape) == 1:
        lam = lam.reshape(1,lam.shape[0])
        
    '''
    Init
    '''
    
    L = n
    clm = np.zeros((L+1, L+1), dtype='longdouble')
    slm = np.zeros((L+1, L+1), dtype='longdouble')
    
    
    '''
    First step of analysis
    '''
    m = np.arange(L+1).reshape(1,L+1)
    c = np.cos((lam.T @ m) * np.pi/180)
    s = np.sin((lam.T @ m) * np.pi/180)
    
    
    '''
    % preserving the orthogonality (except for 'mean' case)
% we distinguish between 'block' and 'pole' type grids (in lambda)
    '''
    
    if (grid == 'block') or (grid == 'cell'):
        if method == 'mean':
            dl = dt
            c[:,0] = dl / 360
            m = np.arange(1, L+1)
            ms = 2 / m * np.sin(m * dl/2 * np.pi/180) / np.pi
            c[:,1:(L+1)+1] = c[:,1:(L+1)+1] * ms  #Double check how to run these lines 2022-08-26
            s[:,1:(L+1)+1] = s[:,1:(L+1)+1] * ms
#            raise Exception("This use case is not yet fully developed")
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
        c[:,[0, L]] = c[:,[0, L]]/2	#% / (1 + dm0 + dmL)
        s[:,[0, L]] = np.zeros((2*n,2))	#% / (1 - dm0 - dmL), Sl0 & SLL unestimable  
    
    
    a = f @ c
    b = f @ s
    
    
    
    
    '''
    Second step of analysis: Clm and Slm
    '''
    if method == 'ls':
        for m in range(L+1):
#            l = np.arange(m,L+1)
            l = np.arange(m,L+1).reshape(L+1-m,1)
            l = l.T
            
            p = plm(l,m,theRAD, 3,1)
            p = p[:,:,0]
            ai = a[:,m]
            bi = b[:,m]
            
            clm[m+1:L+2, m+1] = linalg.lstsq(p, ai)
            slm[m+1:L+2, m+1] = linalg.lstsq(p, bi)
            
#    elif method == 'wls': #wls method not yet checked
#        si = np.sin(theRAD)
#        si = 2 * si / np.sum(si)
#        
#        for m in range(L+1):
#            l = np.arange(m,L+1).reshape(L+1-m,1)
#            l = l.T
#            
#            p = plm(l,m,theRAD, 3,1)
#            ai = a[:,m+1]
#            bi = b[:,m+1]
#            d = np.arange(len(theRAD))
#            pts = p.T @ sparse.coo_matrix((si,(d,d))
#                
##            How to write a sparse:     sparse.coo_matrix((data,(row,col)),shape=(lmax,lmax)).toarray()
#            clm[m+1:L+2, m+1] = linalg.lstsq(pts @ p, pts @ ai)
#            slm[m+1:L+2, m+1] = linalg.lstsq(pts @ p, pts @ bi)
            
    elif method == 'aq': #Approximate Quadrature
        si = np.sin(theRAD)
        si = 2 * si / np.sum(si)
        
        for m in range(L+1):
            l = np.arange(m,L+1).reshape(L+1-m,1)
            l = l.T
            
            p = plm(l,m,theRAD, 3,1)
            
            ai = a[:,m]
            bi = b[:,m]
                        
            clm[m:L+1, m] = (1 + (m==0))/ 4 * p.T @ (si * ai)
            slm[m:L+1, m] = (1 + (m==0))/ 4 * p.T @ (si * bi)
    
    elif method == 'fnm': #1st Neumann method (exact upto L/2)
        w = neumann(np.cos(theRAD))
        
        for m in range(L+1):
            l = np.arange(m,L+1).reshape(L+1-m,1)
            l = l.T
            
            p = plm(l,m,theRAD, 3,1)
            
            ai = a[:,m]
            bi = b[:,m]
                        
            clm[m:L+1, m] = (1 + (m==0))/ 4 * p.T @ (w * ai)
            slm[m:L+1, m] = (1 + (m==0))/ 4 * p.T @ (w * bi)
        
    elif method == 'snm': #2nd Neumann method (exact)
        for m in range(L+1):
            l = np.arange(m,L+1).reshape(L+1-m,1)
            l = l.T
            
            p = plm(l,m,theRAD, 3,1)
            
            ai = a[:,m]
            bi = b[:,m]
            
            clm[m:L+1, m] = (1 + (m==0))/ 4 * p.T @ (gw * ai)
            slm[m:L+1, m] = (1 + (m==0))/ 4 * p.T @ (gw * bi)
            
    elif method == 'mean':
        for m in range(L+1):
            print(m)
            #l = np.arange(m,L+1).reshape(L+1-m,1)
            #l = l.T
            
            
            l = np.array([np.arange(m,L+1,1)])
#            l = np.array([[m]])
            
            p = iplm(l,m,theRAD)
#            p = p[:,-1]
            ai = a[:,m]
            bi = b[:,m]
            
            clm[m:L+1, m] = (1 + (m==0))/ 4 * p.T @ ai
            slm[m:L+1, m] = (1 + (m==0))/ 4 * p.T @ bi
            
            
            
    '''
    %--------------------------------------------------------------------------
% Write the coefficients Clm & Slm in |C\S| format
%--------------------------------------------------------------------------
    '''
        
    #Double check steps here
    slm = np.fliplr(slm)
    cs = sc2cs(np.concatenate((slm[:, np.arange(L)], clm), axis = 1))
    cs = cs[:int(lmax+1), :int(lmax+1)]
    
    '''
    np.save('/mnt/Data/5hk/Project_STC/Mat2Py/mat2py/test/gddc_csRb_r0/csRb_20220922a.npy',cs)
    '''
    return cs
    