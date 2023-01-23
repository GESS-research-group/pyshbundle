#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 15:55:22 2022

% GSHS global spherical harmonic synthesis 
% f = gshs(field)
%
% IN:
%    field ... set of SH coefficients, either in SC-triangle or CS-square format 
%    quant ... optional string argument, defining the field quantity:
%              - 'none' ............. (default), coefficients define the output
%              - 'geoid' ............ geoid height [m],
%              - 'potential' ........ potential [m^2/s^2],
%              - 'dg', 'gravity' .... gravity anomaly [mGal], 
%              - 'tr' ............... grav. disturbance, 1st rad. derivative [mGal],
%              - 'trr' .............. 2nd rad. derivative [1/s^2],
%              - 'water' ............ equiv. water height [m],
%              - 'smd' .............. surface mass density [kg/m^2]. 
%     grd .... optional string argument, defining the grid:
%              - 'pole', 'mesh' ..... (default), equi-angular (n+1)*2n, including 
%                                     poles and Greenwich meridian.
%              - 'block', 'cell' .... equi-angular block midpoints. n*2n
%              - 'neumann', 'gauss' . Gauss-grid (n+1)*2n
%     n ...... grid size parameter n. (default: n = lmax, determined from field)
%              #longitude samples: 2*n
%              #latitude samples n ('blocks') or n+1.
%     h ...... (default: 0), height above Earth mean radius [m].
%     jflag .. (default: true), determines whether to subtract GRS80.
%
% OUT:
%    f ....... the global field
%    theRAD .. vector of co-latitudes [rad] 
%    lamRAD .. vector of longitudes [rad]
%
% EXAMPLE: see SHBUNDLE/examples/example_gshs.m
%
% USES:
%    vec2cs, cs2sc, eigengrav, plm, normalklm, 
%    uberall/grule, uberall/standing, uberall/isint, uberall/ispec
% 
% SEE ALSO:
%    GSHSAG, RSHS, GSHA

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Nico SNEEUW (NS), IAPG, TU-Munich
%    Matthias WEIGELT (MW), DoGE, UofC  
%    Markus ANTONI (MA), GI, Uni Stuttgart 
%    Matthias ROTH (MR), GI, Uni Stuttgart
%    Balaji DEVARAJU (BD), GI, Uni Stuttgart
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    2014-03-09: BD, changed the variable 'grid' to 'grd' as 'grid' conflicts
%                    with the function 'grid'
%    2014-01-14: MR, brush up help text, beautify code, exchange deprecated
%                    'isstr' by 'ischar'
%    2013-02-13: MR, change function names, brush up comments
%    2013-01-30: MA, comments/removing of smoothing option
%    2013-01-23: MA, output in radian
%    2013-01-18: MA, replacement: isscal -> isscalar
%    1998-08-28: NS, brushing up and redefinition of checks
%    1994-??-??: NS, initial version
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

% -------------------------------------------------------------------------
% Some checking, default settings and further preliminaries.
% A lot of checking is done by EIGENGRAV as well.
% -------------------------------------------------------------------------

@author: Amin Shakya, ICWaR, IISc Bangalore
"""
def gshs(field, quant = 'none', grd = 'mesh', n = -9999, h = 0, jflag = 1):
    import numpy as np
    from os import chdir, getcwd
    
    # wd = getcwd()
    # chdir(wd)

    from . import cs2sc
    from . import normalklm
    from . import plm
    from . import eigengrav
    from . import ispec
    
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
        theRAD = np.arange(0, np.pi+dt*0.5, dt, dtype='longdouble')
        lamRAD = np.arange(0, 2*np.pi, dt, dtype='longdouble')
    elif grd == 'block' or grd == 'cell':
        theRAD = np.arange(dt/2, np.pi + dt*0.5, dt, dtype='longdouble')
        lamRAD = np.arange(dt/2, 2*np.pi + dt*0.5, dt, dtype='longdouble')
    #The elif below does not work
#    elif grd == 'neumann' or grd == 'gauss':
#        gx, _ = grule(n+1)
#        theRAD = np.flipud(np.cosh(np.standing(gx)))
#        lamRAD = np.arange(0, 2*np.pi-0.5*dt, dt)
    else:
        raise Exception("Incorrect grid type input")
    
    nlat = len(theRAD)
    nlon = len(lamRAD)
    
    
    '''
    % -------------------------------------------------------------------------
% Preprocessing on the coefficients: 
%    - subtract reference field (if jflag is set)
%    - specific transfer
%    - upward continuation
% -------------------------------------------------------------------------
    '''
    
    if jflag:
        field = field - cs2sc(normalklm(lmax+1))
        
    l = np.arange(0, lmax+1)
    transf = np.array([eigengrav(lmax, quant, h)], dtype='longdouble').T
    
    field = field * np.matmul(transf, np.ones((1, 2*lmax+1), dtype = 'longdouble'), dtype='longdouble') #removed transf[0] 20221019
    
    
    '''
    % -------------------------------------------------------------------------
% Size declarations and start the waitbar:
% Note that the definition of dlam causes straight zero-padding in case N > L.
% When N < L, there will be zero-padding up to the smallest integer multiple
% of N larger than L. After the Fourier transformation (see below), the
% proper samples have to be picked out, with stepsize dlam.
% -------------------------------------------------------------------------
    '''
    
    dlam = np.ceil(lmax/n) #longitude step size
    abcols = dlam*n + 1 #columns required in A and B
    a = np.zeros((nlat, int(abcols)), dtype='longdouble')
    b = np.zeros((nlat, int(abcols)), dtype='longdouble')
    
    '''
    %hwb    = waitbar(0, 'Percentage of orders m ready ...');
%set(hwb, 'NumberTitle', 'off', 'Name', 'GSHS')

% -------------------------------------------------------------------------
% Treat m = 0 separately
% -------------------------------------------------------------------------
    '''
    
    m = 0
    if len(field.shape) == 3: #handle cases where field is a tuple
        c = field[0,m:lmax+1, lmax+m] 
    else:
        c = field[m:lmax+1, lmax+m] 
    l = np.array([np.arange(m,lmax+1)])
    p = plm(l, m, theRAD, nargin = 3, nargout = 1)[:,:,0]
    a[:, m] = np.dot(p,c) 
    b[:, m] = np.zeros(nlat) 
    
    '''
    %waitbar((m+1) / (lmax+1)) 			% Update the waitbar

% -------------------------------------------------------------------------
% Do loop m = 1 to lmax
% -------------------------------------------------------------------------
    '''
    
    for m in range(1,lmax+1,1):
        if len(field.shape) == 3:
            c = field[0,m:lmax+1,lmax+m]
            s = field[0,m:lmax+1,lmax-m]
        else:
            c = field[m:lmax+1,lmax+m]
            s = field[m:lmax+1,lmax-m]
        
        l = np.array([np.arange(m,lmax+1)])
        p = plm(l, m, theRAD, nargin = 3, nargout = 1)[:,:,0]
        a[:, m] = np.dot(p,c)
        b[:, m] = np.dot(p,s)
        
    del field
    '''
    % -------------------------------------------------------------------------
% The second synthesis step consists of an inverse Fourier transformation
% over the rows of a and b. 
% In case of 'block', the spectrum has to be shifted first.
% When no zero-padding has been applied, the last b-coefficients must be set to
% zero. Otherwise the number of longitude samples will be 2N+1 instead of 2N.
% For N=L this corresponds to setting SLL=0!
% -------------------------------------------------------------------------
    '''
    
    
    

    if grd =='block' or grd == 'cell': 
      m      = np.arange(0,abcols,1)
      cshift = np.array([np.ones(nlat)], dtype='longdouble').T * np.array([np.cos(m*np.pi/2/n)], dtype='longdouble');	# cshift/sshift describe the 
      sshift = np.array([np.ones(nlat)], dtype='longdouble').T * np.array([np.sin(m*np.pi/2/n)], dtype='longdouble');	# half-blocksize lambda shift.
      atemp  =  cshift*a + sshift*b;
      b      = -sshift*a + cshift*b;
      a      = atemp;
    
    
    
    if np.remainder(n,lmax) == 0:               #Case without zero-padding
        b[:,int(abcols-1)] = np.zeros(nlat)
    
    #Code for ispec
#    '''
    f = ispec(a.T,b.T).T
    if dlam > 1: 
        f = f[:,np.arange(1,dlam*nlon+1,dlam)]
                                              #This line needs to be worked on
#'''
    return f, theRAD, lamRAD
    

