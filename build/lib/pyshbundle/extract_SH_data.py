#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Created on Fri Dec  9 10:08:55 2022
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
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

# author: Vivek Yadav, Interdisciplinary Center for Water Research (ICWaR), Indian Institute of Science Bengaluru
# Debug and restructure: Abhishek Mhamane, MS-Research Geoinformatics, IIT Kanpur
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

import numpy
import gzip
import re

# Function to read files & extract data
def extract_SH_data(file_path, source):
    # Initialize an empty dictionary to store the coefficients and dates
    data = {
        'coefficients': {},
        'time_coverage_start': None,
        'time_coverage_end': None
    }


    # Regular expression pattern to match the lines with coefficients
    coeff_pattern_csr = re.compile(r'^GRCOF2\s+(\d+)\s+(\d+)\s+([-+]?\d*\.\d+E[-+]?\d+)\s+([-+]?\d*\.\d+E[-+]?\d+)\s+([-+]?\d*\.\d+E[-+]?\d+)\s+([-+]?\d*\.\d+E[-+]?\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+')
    coeff_pattern_jpl = re.compile(r'^GRCOF2\s+(\d+)\s+(\d+)\s+([-+]?\d*\.\d+e[-+]?\d+)\s+([-+]?\d*\.\d+e[-+]?\d+)\s+([-+]?\d*\.\d+e[-+]?\d+)\s+([-+]?\d*\.\d+e[-+]?\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)')
    coeff_pattern_itsg = re.compile(r'^gfc\s+(\d+)\s+(\d+)\s+([-+]?\d*\.\d+e[-+]?\d+)\s+([-+]?\d*\.\d+e[-+]?\d+)\s+([-+]?\d*\.\d+e[-+]?\d+)\s+([-+]?\d*\.\d+e[-+]?\d+)$')


    if source=='jpl': coeff_pattern=coeff_pattern_jpl
    elif source=='csr': coeff_pattern=coeff_pattern_csr
    elif source=='itsg': coeff_pattern=coeff_pattern_itsg


    # Regular expression patterns to match the time coverage start and end lines
    start_pattern = re.compile(r'time_coverage_start\s*:\s*([\d\-T:.]+)')
    end_pattern = re.compile(r'time_coverage_end\s*:\s*([\d\-T:.]+)')
    timeindex_itsg = re.compile(r'^modelname\s+(.+)$')


    # Open and read the gzipped file to extract the time coverage start and end dates
    if source=='itsg':
        with open(file_path, 'rt') as file:
            for line in file:
                # Strip any leading/trailing whitespace characters
                line = line.strip()

                # Search for time coverage start
                start_match = timeindex_itsg.search(line)
                if start_match:
                    data['time_coverage_start'] = start_match.group(1)
                
                # Break the loop if both dates are found
                if data['time_coverage_start']:
                    break
            # File is automatically closed here due to the 'with' statement
        with open(file_path, 'rt') as file:
            for line in file:
                # Strip any leading/trailing whitespace characters
                line = line.strip()
                # print(line)

                # Search for the coefficient pattern in the line
                coeff_match = coeff_pattern.search(line)
                if coeff_match:
                    # Extract degree, order, Clm, and Slm
                    degree = int(coeff_match.group(1))
                    order = int(coeff_match.group(2))
                    clm = numpy.longdouble(coeff_match.group(3))
                    slm = numpy.longdouble(coeff_match.group(4))
                    clm_sdev = numpy.longdouble(coeff_match.group(5))
                    slm_sdev = numpy.longdouble(coeff_match.group(6))
                    
                    # Store the coefficients in the dictionary
                    data['coefficients'][(degree, order)] = {'Clm': clm, 'Slm': slm,
                                                            'Clm_sdev': clm_sdev, 'Slm_sdev': slm_sdev}



    elif source=='csr' or source=='jpl':
        with gzip.open(file_path, 'rt') as file:   # gzip.open
            for line in file:
                # Strip any leading/trailing whitespace characters
                line = line.strip()

                # Search for time coverage start
                start_match = start_pattern.search(line)
                if start_match:
                    data['time_coverage_start'] = start_match.group(1)
                
                # Search for time coverage end
                end_match = end_pattern.search(line)
                if end_match:
                    data['time_coverage_end'] = end_match.group(1)
                
                # Break the loop if both dates are found
                if data['time_coverage_start'] and data['time_coverage_end']:
                    break
            # File is automatically closed here due to the 'with' statement


        # Open and read the gzipped file again to extract the coefficients
        with gzip.open(file_path, 'rt') as file:
            for line in file:
                # Strip any leading/trailing whitespace characters
                line = line.strip()
                # print(line)

                # Search for the coefficient pattern in the line
                coeff_match = coeff_pattern.search(line)
                if coeff_match:
                    # Extract degree, order, Clm, and Slm
                    degree = int(coeff_match.group(1))
                    order = int(coeff_match.group(2))
                    clm = numpy.longdouble(coeff_match.group(3))
                    slm = numpy.longdouble(coeff_match.group(4))
                    clm_sdev = numpy.longdouble(coeff_match.group(5))
                    slm_sdev = numpy.longdouble(coeff_match.group(6))
                    
                    # Store the coefficients in the dictionary
                    data['coefficients'][(degree, order)] = {'Clm': clm, 'Slm': slm,
                                                            'Clm_sdev': clm_sdev, 'Slm_sdev': slm_sdev}
    return data