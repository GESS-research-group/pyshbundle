#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  9 11:43:15 2022

@author: wslvivek
"""
# CLM to C|S format
def clm2cs(data):
    import numpy as np
    
    
    ''' Load data '''
    #data = np.load(path, allow_pickle=True)
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
    ''' clm >>> cs '''
    month = 0
    Lmax = degree[0][-1]
    cs_mat = np.zeros([int(month_count),Lmax+1,Lmax+1])
    for year in range(0,no_of_years,1):
        for tile in range(0,int(len(clm[year])/4752),1):
            i,j = 0,0 
            for index1 in range(2,Lmax+1,1):
                for index2 in range(0,index1+1,1):
                    cs_mat[month,index1,index2] = clm[year][i + tile*4752]
                    i = i + 1
            for index3 in range(2,Lmax+1,1):
                for index4 in range(0,index3,1):
                    cs_mat[month,index4,index3] = slm[year][j+1 +tile*4752]
                    j = j + 1
                j = j + 1
            month = month + 1
    print('Conversion into clm format complete')        
    #np.save('/home/wslvivek/Desktop/level2/preprocess/output/SH_coeff_cs.npy', cs_mat)
