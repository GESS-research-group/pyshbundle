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


This file is part of PySHbundle. 
    PySHbundle is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

Acknowledgement Statement:
    Please note that PySHbundle has adapted the following code packages, 
    both licensed under GNU General Public License
    1. SHbundle: https://www.gis.uni-stuttgart.de/en/research/downloads/shbundle/

    2. Downscaling GRACE Total Water Storage Change using 
    Partial Least Squares Regression
    https://springernature.figshare.com/collections/Downscaling_GRACE_Total_Water_Storage_Change_using_Partial_Least_Squares_Regression/5054564 
    
Key Papers Referred:
    1. Vishwakarma, B. D., Horwath, M., Devaraju, B., Groh, A., & Sneeuw, N. (2017). 
    A data‐driven approach for repairing the hydrological catchment signal damage 
    due to filtering of GRACE products. Water Resources Research, 
    53(11), 9824-9844. https://doi.org/10.1002/2017WR021150

    2. Vishwakarma, B. D., Zhang, J., & Sneeuw, N. (2021). 
    Downscaling GRACE total water storage change using 
    partial least squares regression. Scientific data, 8(1), 95.
    https://doi.org/10.1038/s41597-021-00862-6 

@author: Amin Shakya, ICWAR, Indian Institute of Science
"""

def ispec(a,b = -9999):
    import numpy as np
    import scipy
    import scipy.fft
    
    n2 = a.shape[0]
    a[0,:] = a[0, :]*2

    
    if (np.absolute(b[n2-1,:]) < 1e-10).all():
        n = 2 * n2 - 2     
        a[n2-1,:] = a[n2-1,:] * 2            
        fs = (a - 1j * b)/2
        fs  = (np.concatenate((fs,np.conj(fs[np.arange(n2-2,0,-1),:])), axis = 0))*max(n,1)

    else:
        n = 2 * n2 - 1                        
        fs = (a - 1j * b)/2
        fs = (np.concatenate((fs,np.conj(fs[np.arange(n2-1,0,-1),:])), axis = 0))*n

    f = np.real(scipy.fft.ifft(fs.T).T)
    return f


        
    