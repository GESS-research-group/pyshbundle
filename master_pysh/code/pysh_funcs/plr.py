#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 14:55:09 2022

%
% Input: L [t x s]
%        S [t x n]
%        r  number of selected modes

        aw = area vector for grid cells belonging to the catchment 
        nG = no of GRACE products
        
% Output: H [s x n] model parameters
%

@author: Amin Shakya
"""

def plr(L, S, r, aw, nG):
    import numpy as np
    from scipy import linalg
    
    t = L.shape[0]
    
    #Covariance C [s x n]
    C = np.dot(L.T, S/t)
    
    U, s, Vh = linalg.svd(C) #svd is a inbuilt function in Matlab
    V = Vh.T
    sum_S = np.sum(s)
    
    
    U = U[:, :r]
    sigma = np.diag(s)
    sigma = sigma[:r,:r]
    V = V[:, :r]

    ratio_Sigma = np.sum(s) / sum_S
    
    Ul = np.dot(L,U)
    K = np.linalg.lstsq(Ul, S)[0]
    
    Ut = np.vstack([Ul, L@U]) #Doublecheck this line
    P = np.array(S, (L[: -nG] @ aw.T) / (aw.T @ aw)) #Double-check aw definition
    Kt = np.linalg.lstsq(Ut, P)
    
    #Ht [s x n]
    
    H = np.dot(U, Kt)
    
    #Error covariance #Check order of operation
    Css = np.dot(S.T, S / t)
    Cll = np.dot(L.T, L / t)
    Csl = np.dot(S.T, L / t)
    Cee = Css - Csl @ U @ np.invert(U.T @ Cll @ U) @ U.T @ Csl.T
    
    return H, ratio_Sigma, Cee