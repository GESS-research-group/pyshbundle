#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 15:31:15 2022

% This function computes Gauss base points and weight factors
% using the algorithm given by Davis and Rabinowitz in 'Methods
% of Numerical Integration', page 365, Academic Press, 1975.

% [bp, wf] = grule(n)


% IN:
%    n ....... number of base points required. 
%
% OUT:
%    bp ...... cosine of the base points
%    wf ...... weight factors for computing integrals and such 


@author: 5hk
"""

import numpy as np

def grule(n):
    bp = np.zeros((n,1))
    wf = bp
    iter=2
    m = np.floor((n+1)/2)
    e1 = n * (n+1)
    
    
    mm = 4*m - 1
    t = (np.pi / (4*n + 2)) * np.arange(3,mm+4,4)
    nn = (1 - (1-1/n)/(8*n*n))
    x0 = nn * np.cos(t)
    
    
    for i in range(iter):
        pkm1 = 1
        pk = x0
        
        for kk in range(n-1):
            k = kk+2
            t1 = x0 * pk
            pkp1 = t1-pkm1-(t1-pkm1)/k +t1
            pkm1=pk
            pk=pkp1
            
        den = 1 - x0*x0
        d1 = n * (pkm1 - x0*pk)
        dpn = d1/den
        
        
        d2pn = (2*x0*dpn - e1*pk) / den
        d3pn = (4*x0*d2pn + (2-e1)*dpn)/den
        d4pn = (6*x0*d3pn + (6-e1)*d2pn)/den
        u = pk/dpn
        v = d2pn/dpn
        h = -u * (1+(.5*u)*(v+u*(v*v - u*d3pn/(3*dpn))))
        p = pk + h*(dpn+(0.5*h)*(d2pn+(h/3)*(d3pn + 0.25*h*d4pn)))
        dp = dpn + h*(d2pn+(0.5*h)*(d3pn+h*d4pn/3))
        h=h-p/dp
        x0=x0+h
    
    bp = -x0-h
    fx = d1 - h*e1*(pk+(h/2)*(dpn+(h/3)*(d2pn+(h/4)*(d3pn+(0.2*h)*d4pn))))
    wf = [2 * (1 - np.power(bp,2))]/(fx*fx)
                 
                
    for i in range(len(bp),n):
        bp = np.append(bp,[0])
        wf = np.append(wf,[0])
    
    if ((m)+(m)) != (n):
        m = m-1
    
    #jj = np.arange(1,m+1)
    #nlj = n+1 - jj
    #bp[nlj-1] = -bp[jj-1]
    
    for i in range(1,int(m+1)):
        bp[-i]=-bp[i-1]
        wf[-i] = wf[i-1] 
    return bp, wf