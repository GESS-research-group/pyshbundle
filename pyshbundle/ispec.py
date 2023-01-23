#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 10:48:53 2022
% ISPEC(A,B) returns the function F from the spectra A and B.
%
%IN:
%    a ...... cosine coefficients 
%    b ...... sine coefficients          
%             
%             a and b are defined by:
%             f(t) = a_0 + SUM_(i=1)^n2 a_i*cos(iwt) + b_i*sin(iwt)
%   
%             with w = ground-frequency and n2 half the number of samples (+1).
%             Note that no factor 2 appears in front of the sum.
% 
% OUT:
%    F = ISPEC(A,B) considers A the cosine- and B the sine-spectrum.
%    F = ISPEC(S) assumes S = [A B].
%    If  A and B are matrices, Fourier operations are columnwise.
% 
% USES: 
%    spec
%
% SEE ALSO:
%    SPEC, FFT

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Nico SNEEUW (NS), IAPG, TU-Munich
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    1994-06-29: NS, initial version
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
@author: Amin Shakya, ICWAR, Indian Institute of Science
"""
import numpy as np
import scipy
import scipy.fft

def ispec(a,b = -9999):
    '''
    if (b == -9999).all():                          #Only one input
        if min(a.shape) == 2 and a.shape[1] == 2:
            a = a.T
        
        m = a.shape[1]
        
        if(np.remainder(m,2) != 0):
            raise Exception("If one input argument, number of columns must be even")
        
        b = a[:, m/2 : m + 1]     
        a = a[:, 0:m/2]     
    else:
        if (a.shape[0] != b.shape[0]) or (a.shape[1] != b.shape[1]):
            raise Exception("Size of a and b do not match")
        
        if min(a.shape) == 1:
            a = a.T                              #Put a and b upright
            b = b.T
    '''
    
    n2 = a.shape[0]
    a[0,:] = a[0, :]*2

    
    if (np.absolute(b[n2-1,:]) < 1e-10).all():
        n = 2 * n2 - 2     
        a[n2-1,:] = a[n2-1,:] * 2               #Simulate 100% aliasing
        fs = (a - 1j * b)/2
        fs  = (np.concatenate((fs,np.conj(fs[np.arange(n2-2,0,-1),:])), axis = 0))*max(n,1)

    else:
        n = 2 * n2 - 1                          #Double check this
        fs = (a - 1j * b)/2
        fs = (np.concatenate((fs,np.conj(fs[np.arange(n2-1,0,-1),:])), axis = 0))*n
    
        
    
    #fs = np.csingle(fs)
    f = np.real(scipy.fft.ifft(fs.T).T)
    return f

# import scipy.io
# err_dict = {"fsf": np.imag(fs)}
# scipy.io.savemat('fsf.mat', err_dict)




        
    