# -*- coding: utf-8 -*-
"""
% This function uses cubic interpolation to replace NaNs
% See INTERP1 for more info
"""

def naninterp(X):
    from scipy.interpolate import PchipInterpolator
    import numpy as np
    nan = np.nan
    
    ok = ~np.isnan(X)
    xp = ok.ravel().nonzero()[0] #Indices of xs with values
    fp = X[~np.isnan(X)]
    
    x  = np.isnan(X).ravel().nonzero()[0] #Indices of xs without values
    
    pchip = PchipInterpolator(xp,fp) #Initialize scipy PHCIP cubic interpolation
    X[np.isnan(X)] = pchip(x) #Interpolate Nan values in X
    
    return X
    

    