#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  9 10:14:32 2022

@author: wslvivek
"""
# CLM to /S|C\
def clm2sc(data):
    import numpy as np
    import pickle
    # data=np.load("/home/wslvivek/Desktop/level2/pysh_v2/output/saved_as-num.npy",allow_pickle=True)
    # with open("/home/wslvivek/Desktop/level2/pysh_v2/output/saved_as_num", "rb") as pk:
    #     data=pickle.load(pk)
    ''' Read variables '''
    no_of_years = len(data[0])
    degree = data[0]
    clm = data[2]
    slm = data[3]

    Lmax=degree[0][-1]
    degree_order=int((Lmax+1)*(Lmax+2)/2)
    ''' Count no of months of data '''
    month_count =0
    for i in range(0,len(data[0]),1):
        month_count= month_count+round(len(data[0][i])/degree_order)
    
    ''' clm >>> sc '''
    month = 0
    sc_mat = np.zeros([month_count,Lmax+1,2*Lmax+2])
    for year in range(0,no_of_years,1):
        for tile in range(0,int(len(degree[year])/degree_order),1):
            i = 0
            for index1 in range(0,Lmax+1,1):
                for index2 in range(0,index1+1,1):
                    sc_mat[month,index1,Lmax+index2+1] = clm[year][i+tile*degree_order]
                    sc_mat[month,index1,Lmax-index2] = slm[year][i+tile*degree_order]
                    i = i + 1
            month = month + 1
            
    "delete order 0 column"
    sc_mat = np.delete(sc_mat, Lmax, 2)
    print('Conversion into clm format complete')
    return sc_mat

# mean=np.mean(sc_mat[18:102], axis=0)
