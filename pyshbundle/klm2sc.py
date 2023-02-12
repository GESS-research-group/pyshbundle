#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 10:14:27 2023

@author: wslvivek
"""

def klm2sc(data):
    import numpy as np
    
    ''' Read variables '''
    no_of_years = len(data[0])
    degree = data[0]
    clm = data[2]
    slm = data[3]
    
    ''' Count no of months of data '''
    month_count =0
    for i in range(0,len(data[0]),1):
        month_count= month_count+len(data[0][i])/4753
        
    ''' klm >>> sc '''
    month = 0
    Lmax = degree[0][-1]
    sc_mat = np.zeros([int(month_count),Lmax+1,2*Lmax+1])
    for year in range(0,no_of_years,1):
        index2 =0
        for tile in range(0,int(len(clm[year])/4753),1):  
            for index1 in range(0,Lmax+1,1):
                sc_mat[month,index1:,(Lmax-index1)] = slm[year][(index2):(index2+Lmax-index1+1)]
                sc_mat[month,index1:,(index1+Lmax)] = clm[year][(index2):(index2+Lmax-index1+1)]
                index2 = index2+Lmax-index1+1
            month=month+1
            
    return sc_mat