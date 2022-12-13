#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  9 10:14:32 2022

@author: wslvivek
"""
# CLM to /S|C\
def clm2sc(data):
    import numpy as np
    
    ''' Load data '''
    #data = np.load(path+str('saved_as_num.npy'), allow_pickle=True)
    #data = np.load('/home/wslvivek/Desktop/level2/preprocess/output/saved_as_num.npy', allow_pickle=True)
    
    ''' Read variables '''
    no_of_years = len(data[0])
    degree = data[0]
    clm = data[2]
    slm = data[3]
    
    ''' Count no of months of data '''
    month_count =0
    for i in range(0,len(data[0]),1):
        month_count= month_count+len(data[0][i])/4750
    
    ''' clm >>> sc '''
    month = 0
    Lmax = degree[0][-1]
    sc_mat = np.zeros([int(month_count),Lmax+1,2*Lmax+2])
    for year in range(0,no_of_years,1):
        for tile in range(0,int(len(clm[year])/4752),1):
            i = 0
            for index1 in range(1,Lmax+1,1):
                for index2 in range(0,index1+1,1):
                    sc_mat[month,index1,Lmax+index2+1] = clm[year][i+tile*4752]
                    sc_mat[month,index1,Lmax-index2] = slm[year][i+tile*4752]
                    i = i + 1
            month = month + 1
            
    "delete order 0 column"
    sc_mat = np.delete(sc_mat, Lmax, 2)
    ''' Calculation of long-term mean '''
    #SH_long_mean = sc_mat[26:98].mean(0)
    #np.save('/home/wslvivek/Desktop/level2/preprocess/output/SH_long_mean',SH_long_mean)            
    return sc_mat