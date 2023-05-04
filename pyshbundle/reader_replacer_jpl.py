#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Created on Fri Dec  9 10:08:55 2022
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# Acknowledgement Statement:
#    Please note that PySHbundle has adapted the following code packages, 
#    both licensed under GNU General Public License
#       1. SHbundle: https://www.gis.uni-stuttgart.de/en/research/downloads/shbundle/

#       2. Downscaling GRACE Total Water Storage Change using Partial Least Squares Regression
#          https://springernature.figshare.com/collections/Downscaling_GRACE_Total_Water_Storage_Change_using_Partial_Least_Squares_Regression/5054564 
    
# Key Papers Referred:
#    1. Vishwakarma, B. D., Horwath, M., Devaraju, B., Groh, A., & Sneeuw, N. (2017). 
#       A data‐driven approach for repairing the hydrological catchment signal damage 
#       due to filtering of GRACE products. Water Resources Research, 
#       53(11), 9824-9844. https://doi.org/10.1002/2017WR021150

#    2. Vishwakarma, B. D., Zhang, J., & Sneeuw, N. (2021). 
#       Downscaling GRACE total water storage change using 
#       partial least squares regression. Scientific data, 8(1), 95.
#       https://doi.org/10.1038/s41597-021-00862-6 

# @author: Vivek Yadav, Interdisciplinary Center for Water Research (ICWaR), Indian Institute of Science (IISc)

import gzip
import re
import os
import numpy as np
import julian
import math
# Function to read files & extract data
def reader(file_name,line_num,degree,order,clm,slm,slm_std,clm_std,start_date,end_date,year_start,time_axes):
    """[Helper function]

    Args:
        file_name (str): _description_
        line_num (int): _description_
        degree (int): _description_
        order (int): _description_
        clm (_type_): _description_
        slm (_type_): _description_
        slm_std (_type_): _description_
        clm_std (_type_): _description_
        start_date (_type_): _description_
        end_date (_type_): _description_
        year_start (_type_): _description_
        time_axes (_type_): _description_
    """
    with gzip.open(file_name,"r") as file:
        stuff = file.readlines()
        stuff
        for i in range (0,len(stuff),1):
            if str(stuff[i]) == str( b'# End of YAML header\n',):
                line_num = i+1                                                 #Line number of starting line of data
                break
        pattern = '\s+'                                                        #Delimiter for splitting                                                          
        while(line_num<len(stuff)):
            split = re.split(  pattern, str(stuff[line_num]))                  #split wordwise with Regex
            degree[time_axes].append(  int(split[1])  )
            order[time_axes].append( int(split [2])    )
            clm[time_axes].append( float(split [3])    )
            slm[time_axes].append( float(split [4])    )
            slm_std[time_axes].append( float(split [5])    )
            clm_std[time_axes].append(  float(split[6])    )
            start_date[time_axes].append(  split[7][0:4]+'-'+split[7][4:6]+'-'+split[7][6:8]  )
            end_date[time_axes].append(  split[8][0:4]+'-'+split[8][4:6]+'-'+split[8][6:8]  )
            line_num = line_num + 1
           
# Function for yearwise            
def TIME(year_start,file_name,time_axes):
    """[Helper function]

    Args:
        year_start (_type_): _description_
        file_name (_type_): _description_
        time_axes (_type_): _description_

    Returns:
        _type_: _description_
    """
    if  year_start == file_name[-39:-35]:
        time_axes = time_axes
        year_start = year_start
    else:
        time_axes = time_axes + 1  
        year_start = file_name[-39:-35]
    return year_start, time_axes
# Main code
def last_4chars(x):
        #print(x[-39:-32])
        return(x[-39:-32])

def reader_replacer_jpl(path, path_tn14, path_tn13):
    """Replaces the certain coefficients with other given coeff.

    Args:
        path (str): _description_
        path_tn14 (str): _description_
        path_tn13 (str): _description_
    
    Returns:
        _type_: _description_
    """
    file_list = os.listdir(path)

    filenames = os.listdir(path)       #Names of files in folder    
    # Identify the data product source
    if 'GFZ' in str(filenames[0]):
        source = str('GFZ')
    elif 'CSR' in str(filenames[0]):
        source = str('CSR')
    elif 'JPL' in str(filenames[0]):
        source = str('JPL')
    elif 'ITSG' in str(filenames[0]):
        source = str('ITSG')
    print(source)
    aa = sorted(file_list, key = last_4chars) #when does your data start at ?
    year_start = aa[0][6:10]
    time_axes = 0

    # Counts the number of files
    no_of_files = 0

    # Empty numpy arrays to store variables
    # Empty numpy arrays to store variables
    line_num = 0
    oi = 22
    degree = [[] for x in range(oi)]
    order = [[] for x in range(oi)]
    clm = [[] for x in range(oi)]
    slm = [[] for x in range(oi)]
    slm_std = [[] for x in range(oi)]
    clm_std = [[] for x in range(oi)]
    start_date = [[] for x in range(oi)]
    end_date = [[] for x in range(oi)]
    
    # Iterate through all files
    for file in sorted(file_list, key = last_4chars):
        #print(file)
        if file.endswith(".gz"):
            file_name = f"{path}/{file}"
            #print(file_name[-39:-32])
            
            # Save yearwise
            year_start, time_axes = TIME(year_start,file_name,time_axes)
            
            # Call the function 'reader'
            reader(file_name,line_num,degree,order,clm,slm,slm_std,clm_std,start_date,end_date,year_start,time_axes)
            no_of_files = no_of_files + 1
    print('Reading into clm format complete!')
    print("Number of files read:", no_of_files)
    
    Lmax=degree[0][-1]
    degree_order=int((Lmax+1)*(Lmax+2)/2)
    ''' Replacement '''
    print('Starting replacement')
    ''' Replace deg 2,3 '''
    new_file_TN14 = path_tn14
    
    rep_start_date, rep_end_date, c20, c30, c20_sigma, c30_sigma = [], [], [], [], [], []
    with open(new_file_TN14,"r") as file:
        stuff = file.readlines()
        for i in range(0,len(stuff),1):
            line_num = str('not found')
            if  stuff[i] == str('Product:\n'):
                print("found:",i+1)
                line_num = i + 1
                break
        if type(line_num) is str:
            print('Replacement data not found')
        pattern = '\s+'
        count = 0
        while (line_num<len(stuff)):
            split = re.split(pattern, str(stuff[line_num]))
            c20.append( float(split [2]) )
            c30.append( float(split[5]) )
            c20_sigma.append(float(split[4])*1e-10)
            c30_sigma.append(float(split[7])*1e-10)
            rep_start_date.append(str(julian.from_jd(float(split[0]), fmt='mjd').date()))
            rep_end_date.append(str(julian.from_jd(float(split[8]), fmt='mjd').date()))
            line_num = line_num + 1
            count = count + 1
            
    # Actual replacement    
    
    index = 0
    for year in range(0,time_axes+1,1):
        for y in range(0,round(len(degree[year])/int(degree_order-3)),1):    
            while(index<len(c20)):
                
                curr_date = start_date[year][y*degree_order]
                rep_date_cmp = rep_start_date[index]
                
                if curr_date == rep_date_cmp:
                    print(curr_date, rep_date_cmp, index)
                    clm[year][y*int(degree_order-3)] = c20[index]
                    clm_std[year][y*int(degree_order-3)]=c20_sigma[index]
                    if math.isnan(c30[index]) == False:
                        clm[year][y*int(degree_order-3)+3] = c30[index]    
                        clm_std[year][y*int(degree_order-3)+3]=c30_sigma[index]
                    index= index + 1
                    break
                else:   
                    index = index +1
    print('Degree 2,3 replacement complete!')
                    
    ''' Replace deg 1 '''
    new_file_TN13 = path_tn13
    rep_start_date_deg1, rep_end_date_deg1, c1m, s1m, c1m_sigma, s1m_sigma = [], [], [], [], [], []
    with open(new_file_TN13,"r") as file:
        stuff = file.readlines()
        for i in range(0,len(stuff),1):
            if  stuff[i] == str('end of header ===============================================================================\n'):
                print("found:",i+1)
                line_num = i + 1
                break
            else:
                line_num = str('not found')
        pattern = '\s+'
        count = 0
        while (line_num<len(stuff)):
            split = re.split(pattern, str(stuff[line_num]))
            c1m.append( float(split [3]) )
            s1m.append( float(split[4]) )
            c1m_sigma.append( float(split[5])*1e-10)
            s1m_sigma.append( float(split[6])*1e-10)
            rep_start_date_deg1.append(split[7][0:4] +'-'+split[7][4:6] +'-'+split[7][6:8])
            rep_end_date_deg1.append(split[8][0:4] +'-'+split[8][4:6] +'-'+split[8][6:8])
            line_num = line_num + 1
            count = count + 1 
      
    # replace deg 1   
    index = 0
    for year in range(0,time_axes+1,1):
        for y in range(0,int(len(degree[year])/int(degree_order-3)),1):    
            while(index<len(c1m)):
                
                curr_date_deg1 = start_date[year][y*int(degree_order-1)]
                rep_date_cmp_deg1 = rep_start_date_deg1[index]
                
                if curr_date_deg1 == rep_date_cmp_deg1:
                    print(curr_date_deg1, rep_date_cmp_deg1, index)
                    degree[year].insert(y*int(degree_order-1),int(1))
                    degree[year].insert(y*int(degree_order-1)+1,int(1))
                    order[year].insert(y*int(degree_order-1),int(0))
                    order[year].insert(y*int(degree_order-1)+1,int(1))
                    clm[year].insert(y*int(degree_order-1),float(c1m[index]))
                    clm[year].insert(y*int(degree_order-1)+1,float(c1m[index+1]))
                    slm[year].insert(y*int(degree_order-1),float(s1m[index]))
                    slm[year].insert(y*int(degree_order-1)+1,float(s1m[index+1]))
                    clm_std[year].insert(y*int(degree_order-1),float(c1m_sigma[index]))
                    clm_std[year].insert(y*int(degree_order-1)+1,float(c1m_sigma[index+1]))
                    slm_std[year].insert(y*int(degree_order-1),float(s1m_sigma[index]))
                    slm_std[year].insert(y*int(degree_order-1)+1,float(s1m_sigma[index+1]))
                    start_date[year].insert(y*int(degree_order-1),rep_start_date_deg1[index])
                    start_date[year].insert(y*int(degree_order-1)+1,rep_start_date_deg1[index+1])
                    end_date[year].insert(y*int(degree_order-1),rep_end_date_deg1[index])
                    end_date[year].insert(y*int(degree_order-1)+1,rep_end_date_deg1[index+1])
                    index = index + 2
                    break
                else:   
                    index = index +2                
    print('Degree 1 replacement complete!')
    
    
    ''' Add degree 0 coeffs (zeros) '''
    for year in range(0,time_axes+1,1):
        for y in range(0,int(len(degree[year])/int(degree_order-1)),1):
            degree[year].insert(y*degree_order, float(0))
            order[year].insert(y*degree_order, float(0))
            clm[year].insert(y*degree_order, float(0)) 
            slm[year].insert(y*degree_order, float(0))
            clm_std[year].insert(y*degree_order, float(0))
            slm_std[year].insert(y*degree_order, float(0))
            start_date[year].insert(y*degree_order, start_date[year][y*degree_order+1])
            end_date[year].insert(y*degree_order, end_date[year][y*degree_order+1]) 
    
    ''' Save everything in a list '''
    saved_as_num=[degree,order,clm,slm,clm_std,slm_std,start_date,end_date]
    
    ''' Number of years '''
    beta = np.zeros(len(start_date))
    sum = 0
    for x in range(0,time_axes+1,1):
        beta[x] = round(len(start_date[x])/degree_order)
        sum = sum + beta[x]
    
    ''' Finding the dates for time bounds of data '''
    dates_start, dates_end = [],[] 
    for i in range(0,time_axes+1,1):
        j = 0
        while j < beta[i]:
            dates_start.append(str(start_date[i][j*degree_order]))
            dates_end.append(str(end_date[i][j*degree_order]))
            j = j + 1
    print("Number of months of data in each year starting", dates_start[0], \
          "& ending", dates_end[-1], beta)   
        
    # import pickle
    # with open("/home/wslvivek/Desktop/level2/pysh_v2/output/saved_as_num","wb") as pk:
    #     pickle.dump(saved_as_num, pk)
        
    return saved_as_num, dates_start,dates_end, no_of_files