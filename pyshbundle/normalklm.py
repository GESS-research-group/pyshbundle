#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 14:34:55 2022
% NORMALKLM returns an ellipsoidal normal field
% consisting of normalized -Jn, n=0,2,4,6,8
%
% IN:
%    lmax ....... maximum degree
%    typ ........ either 'wgs84' (equipotential ellipsoid), default,
%                        'grs80',
%                 or     'he' (hydrostatic equilibrium ellipsoid)
% OUT:
%    nklm ....... normal field in CS-format (sparse)
%
% REMARKS:
%    .J2,J4 values for hydrostatic equilibrium ellipsoid from Lambeck (1988)
%    "Geophysical Geodesy", p.18

% -------------------------------------------------------------------------
% project: SHBundle 
@author: Amin Shakya
"""

def normalklm(lmax, typ = 'wgs84'):
    import numpy as np 
    from scipy import sparse
    if type(lmax) != int:
        raise Exception("lmax should be integer")
        
    if lmax < 0:
        raise Exception("lmax should be positive")
        
    
    typ_ = typ.lower()
    if (typ_ == 'wgs84'):
        J2     =  1.08262982131e-3     #% earth's dyn. form factor (= -C20 unnormalized)
        J4     = -2.37091120053e-6    #% -C40 unnormalized
        J6     =  6.08346498882e-9     #% -C60 unnormalized
        J8     = -1.42681087920e-11    #% -C80 unnormalized
        jcoefs = np.array([1, -J2, -J4, -J6, -J8]).T.reshape(5,1)
        l      = np.arange(0,min(lmax + 2,8 + 2), 2).T
        l.reshape(l.shape[0],1)
        
    elif (typ_ == 'grs80'):
        J2     =  1.08263e-3         # % earth's dyn. form factor (= -C20 unnormalized)
        J4     = -2.37091222e-6     #% -C40 unnormalized
        J6     =  6.08347e-9        #% -C60 unnormalized
        J8     = -1.427e-11         #% -C80 unnormalized
        jcoefs = np.array([1, -J2, -J4, -J6, -J8]).reshape(5,1)
        l      = np.arange(0,min(lmax + 2,8 + 2), 2).T
        l.reshape(l.shape[0],1)
        
    elif ((typ_ == 'he') or (typ_ == 'hydro')):
        J2     = 1.072618e-3		#% earth's dyn. form factor (= -C20 unnormalized)
        J4     = 0.2992e-5     	#% -C40 unnormalized
        jcoefs = np.array([1, -J2, -J4]).T.reshape(5,1)
        l      = np.arange(0,min(lmax + 2,4 + 2), 2).T
        l.reshape(l.shape[0],1)
        
    else:
        raise Exception("Unknown type of ellipsoid:   ", typ)
    
    coefs = jcoefs[:len(l)].T / np.sqrt(2*l + 1)
#    coefs.reshape(coefs.shape[0],1)
    
    
    data = np.array(coefs)[0]
    row = np.array(l)
    col = np.zeros(len(l))
    nklm = sparse.coo_matrix((data,(row,col)),shape=(lmax,lmax)).toarray()
    return nklm