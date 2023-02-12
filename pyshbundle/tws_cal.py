#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 09:47:38 2022

@author: wslvivek, Amin Shakya
"""

def tws_cal(data,lmax,gs,r,m):
    from gaussian import gaussian
    from gshs import gshs
    import numpy as np

    SC = data
    
    gfilter = gaussian(lmax,r)
    grid_y = int(180/gs)
    grid_x = int(360/gs)
    tws_f = np.zeros([m,grid_y,grid_x], dtype ='longdouble')
    for i in range(0,m,1):
        field = SC[i,0:lmax+1,96-lmax:96+lmax+1]
        shfil = np.zeros([lmax+1,2*lmax+1])
        for j in range(0,2*lmax+1,1):
            shfil[:,j] = gfilter[:,0] * field[:,j]
        quant = 'water' 
        grd = 'cell'
        n = int(180/gs) 
        h = 0 
        jflag = 0
        
        
        ff = gshs(shfil,quant,grd,n,h,jflag)[0]
        
        ff = ff*1000
        tws_f[i,:,0:int(grid_x/2)] = ff[:,int(grid_x/2):]
        tws_f[i,:,int(grid_x/2):] = ff[:,0:int(grid_x/2)]   
    
    import matplotlib.pyplot as plt
    plt.imshow(tws_f[0,:,:])
    return(tws_f)
