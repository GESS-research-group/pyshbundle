#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 18 14:38:16 2022
% NEUMANN returns the weights and nodes for Neumann's numerical integration
%
% NB1 In 1st N-method, length(x) should not become too high, 
% since a linear system of this size is solved. Eg: 500.
% NB2 No use is made of potential symmetries of nodes.
%
% HOW:   w     = neumann(x)   -- 1st Neumann method 
%        [w,x] = neumann(n)   -- 2nd Neumann method (Gauss quadrature)
%
% IN:
%    x ...... base points (nodes) in the interval [-1;1]
%    n ...... number of weights and number of base points
%
% OUT:
%    w ...... quadrature weights
%    x ...... base points (nodes) in the interval [-1;1]
%
% USES:
%    plm, uberall/grule
%
% REMARKS:
%    1st N.-method: see Sneeuw (1994) GJI 118, pp 707-716, eq. 19.5
%    2nd N.-method: see uberall/GRULE

% ----------------------------------------------------------------------------
% project: SHBundle 
% ----------------------------------------------------------------------------
% authors:
%    Nico SNEEUW (NS), IAGP, TU Munich
%    Markus ANTONI (MA), GI, Uni Stuttgart 
%    Matthias ROTH (MR), GI, Uni Stuttgart
%    <bundle@gis.uni-stuttgart.de>
% ----------------------------------------------------------------------------
% revision history:
%    2013-02-13: MR, change function names, brush up comments
%    2013-01-23: MA, angular input of plm in radian
%    2000-02-14: NS, remove dependencies on isinteger, isvector, isscalar
%    1999-06-30: NS, NB1 and NB2 added
%    1998-11-25: NS, interchange names 1st and 2nd
%    1995-08-23: NS, initial version
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
@author: 5hk
"""

import numpy as np
from . import grule
from . import plm
def neumann(inn):
    try: #if input is an integer
        x, w = grule(inn)
    except: #if input is an array
        if(len(inn)==1): #2nd Neumann method
            x,w = grule(inn)
            if(np.not_equal(np.mod(x, 1), 0)): #Not integer
                raise Exception("Integer input argument required")
            
            
        
        elif min(inn.shape) == 1: #1st Neumann method #Size gives 2 outputs for 2d array in matlab; for row and column
            x = inn
            theRAD = np.arccos(x) #x in radian
            l = np.array(list(range(len(x))))
            pp = plm(l, therad)
            
            rr = list([2])
            for i in len(x-1):
                rr.append(0)
            r = np.asarray(rr)
                
            w,resid,rank,s = np.linalg.lstsq(pp,r) #Solve system of equations; Double check this operation
            if(x.shape != w.shape):
                w=w.T
            
        else:
            raise Exception("Error in input dimensions")
    
    return w, x
        

