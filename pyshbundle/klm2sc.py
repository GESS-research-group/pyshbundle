#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 10:14:27 2023

@author: wslvivek
"""

def klm2sc(data):
    import numpy as np
    # import pickle
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
        
    ''' klm >>> sc '''
    month = 0
    sc_mat = np.zeros([month_count,Lmax+1,2*Lmax+2])
    for year in range(0,no_of_years,1):
        index2 =0
        for tile in range(0,int(len(clm[year])/degree_order),1):  
            for index1 in range(0,Lmax+1,1):
                sc_mat[month,index1:,Lmax-index1] = slm[year][(index2):(index2+Lmax-index1+1)]
                sc_mat[month,index1:,index1+Lmax] = clm[year][(index2):(index2+Lmax-index1+1)]
                print(month,'\t',index1,'\t',Lmax-index1,'\t',year,'\t',index2,'\t',index2+Lmax-index1+1)
                index2 = index2+Lmax-index1+1
            month=month+1
    sc_mat=np.delete(sc_mat,193,axis=2)
    print('Conversion into clm format complete')
    return sc_mat