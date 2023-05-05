#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Created on Fri Dec  9 10:08:55 2022

# License:
#    This file is part of PySHbundle.
#    PySHbundle is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
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

"""_Module_Summary_
"""
import gzip
import os
import re
import numpy as np
import julian
import math

def reader(file_name: str,line_num, degree: int, order: int, clm,slm,delta_clm,delta_slm,start_date,end_date,year_start,time_axes):
    """Function to read files & extract data

    Args:
        file_name (_type_): _description_
        line_num (_type_): _description_
        degree (_type_): _description_
        order (_type_): _description_
        clm (_type_): _description_
        slm (_type_): _description_
        delta_clm (_type_): _description_
        delta_slm (_type_): _description_
        start_date (_type_): _description_
        end_date (_type_): _description_
        year_start (_type_): _description_
        time_axes (_type_): _description_
    """
    with gzip.open(file_name,"r") as file:
        stuff = file.readlines()
        stuff
        for i in range (0,len(stuff),1):
            #print(str(line))
            if str(stuff[i]) == str( b'# End of YAML header\n',):
                #print("found:",i+1)
                line_num = i+1                                                 #Line number of starting line of data
                break
        pattern = '\s+'                                                        #Delimiter for splitting                                                          
        while(line_num<len(stuff)):
            split = re.split(  pattern, str(stuff[line_num]))                  #split wordwise with Regex
            degree[time_axes].append(  int(split[1])  )
            order[time_axes].append( int(split [2])    )
            clm[time_axes].append( float(split [3])    )
            slm[time_axes].append( float(split [4])    )
            delta_clm[time_axes].append( float(split [5])    )
            delta_slm[time_axes].append(  float(split[6])    )
            start_date[time_axes].append(  split[7][0:4]+'-'+split[7][4:6]+'-'+split[7][6:8]  )
            end_date[time_axes].append(  split[8][0:4]+'-'+split[8][4:6]+'-'+split[8][6:8]  )
            line_num = line_num + 1


            
# Function for yearwise            
def TIME(year_start,file_name,time_axes):
        if  year_start == file_name[-39:-35]:
            time_axes = time_axes
            year_start = year_start
        else:
            time_axes = time_axes + 1  
            year_start = file_name[-39:-35]
        return year_start, time_axes

def last_4chars(x):
        #print(x[-39:-32])
        return(x[-39:-32])
    
# Main code
def reader_replacer(path, path_tn14, path_tn13):

    # Give path to Level2 data
    #path = r"/home/wslvivek/Desktop/level2/Level_2_Data/JPL_GSM_GRACE"
    file_list = os.listdir(path)    
    
    filenames = os.listdir(path)                                                   #Names of files in folder
    
    
    # Identify the data product source
    if 'GFZ' in str(filenames[1]):
        source = str('GFZ')
    
    elif 'CSR' in str(filenames[1]):
        source = str('CSR')
    
    elif 'JPL' in str(filenames[1]):
        source = str('JPL')
    aa = sorted(file_list, key = last_4chars)                                      #when does your data start at ?
    year_start = aa[1][6:10]
    time_axes = 0
    
    
    # Counts the number of files
    no_of_files = 0
    
    
    # Empty numpy arrays to store variables
    # Empty numpy arrays to store variables
    line_num = 0
    oi = 21
    degree = [[] for x in range(oi)]
    order = [[] for x in range(oi)]
    clm = [[] for x in range(oi)]
    slm = [[] for x in range(oi)]
    delta_clm = [[] for x in range(oi)]
    delta_slm = [[] for x in range(oi)]
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
            reader(file_name,line_num,degree,order,clm,slm,delta_clm,delta_slm,start_date,end_date,year_start,time_axes)
            no_of_files = no_of_files + 1
    print('Reading into clm format complete!')
    print("Number of files read:", no_of_files)
    ''' Replacement '''
    print('Starting replacement')
    ''' Replace deg 2,3 '''
    new_file_TN14 = path_tn14
    
    
    rep_start_date, rep_end_date, c20, c30 = [], [], [], []
    with open(new_file_TN14,"r") as file:
        stuff = file.readlines()
        for i in range(0,len(stuff),1):
            if  stuff[i] == str('Product:\n'):
                #print("found:",i+1)
                line_num = i + 1
                break
            else:
                line_num = str('not found')
        pattern = '\s+'
        count = 0
        while (line_num<len(stuff)):
            split = re.split(pattern, str(stuff[line_num]))
            c20.append( float(split [2]) )
            c30.append( float(split[5]) )
            rep_start_date.append(str(julian.from_jd(float(split[0]), fmt='mjd').date()))
            rep_end_date.append(str(julian.from_jd(float(split[8]), fmt='mjd').date()))
            line_num = line_num + 1
            count = count + 1
            
    # Actual replacement    
    
    index = 0
    for year in range(0,21,1):
        for y in range(0,int(len(clm[year])/4750),1):    
            while(index<205):
                if start_date[year][y*4750] == rep_start_date[index]:
                    clm[year][y*4750] = c20[index]
                    if math.isnan(c30[index]) == False:
                        clm[year][y*4750+3] = c30[index]                   
                    break
                else:   
                    index = index +1
    print('Degree 2,3 replacement complete!')
                    
    ''' Replace deg 1 '''
    

    new_file_TN13 = path_tn13
    
    rep_start_date_deg1, rep_end_date_deg1, c1m, s1m = [], [], [], []
    with open(new_file_TN13,"r") as file:
        stuff = file.readlines()
        for i in range(0,len(stuff),1):
            if  stuff[i] == str('end of header ===============================================================================\n'):
                #print("found:",i+1)
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
            rep_start_date_deg1.append(split[7][0:4] +'-'+split[7][4:6] +'-'+split[7][6:8])
            rep_end_date_deg1.append(split[8][0:4] +'-'+split[8][4:6] +'-'+split[8][6:8])
            line_num = line_num + 1
            count = count + 1 
      
    # replace deg 1   
    index = 0
    for year in range(0,21,1):
        for y in range(0,int(len(clm[year])/4750),1):    
            while(index<424):
                if start_date[year][y*4752] == rep_start_date_deg1[index]:
                    clm[year].insert(y*4752,float(c1m[index]))
                    clm[year].insert(y*4752+1,float(c1m[index+1]))
                    slm[year].insert(y*4752,float(s1m[index]))
                    slm[year].insert(y*4752+1,float(s1m[index+1]))
                    break
                else:   
                    index = index +2            
    saved_as_num = np.array([np.array(degree),np.array(order),np.array(clm),np.array(slm),np.array(delta_clm),np.array(delta_slm),np.array(start_date),np.array(end_date)])            
    print('Degree 1 replacement complete!')
    
    
    ''' Number of years '''
    beta = np.zeros(len(start_date))
    sum = 0
    for x in range(0,21,1):
        beta[x] = (len(start_date[x])/4750)
        sum = sum + len(start_date[x])/4750 
    
    ''' Finding the dates for time bounds of data '''
    dates_start, dates_end = [],[] 
    for i in range(0,21,1):
        j = 0
        while j < beta[i]:
            #print(  str(start_date[i][j*4750]) + '     '+'till' +'     ' + str(end_date[i][j*4750])  )  
            #dates.append(str(start_date[i][j*4750]) +  '     ' + 'till' +'     ' + str(end_date[i][j*4750]))
            dates_start.append(str(start_date[i][j*4750]))
            dates_end.append(str(end_date[i][j*4750]))
            j = j + 1
    print("Number of months of data in each year starting", dates_start[0], \
          "& ending", dates_end[-1], beta)       
    return saved_as_num, dates_start,dates_end, no_of_files
# Saved as numpy array
# saved_as_num = np.array([np.array(degree),np.array(order),np.array(clm),np.array(slm),np.array(delta_clm),np.array(delta_slm),np.array(start_date),np.array(end_date)])
# np.save('/home/wslvivek/Desktop/level2/preprocess/saved_as_num', saved_as_num)
