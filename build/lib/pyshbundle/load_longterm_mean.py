#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Created on Fri Feb  17 2023

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
    
# @author: Amin Shakya, Interdisciplinary Center for Water Research (ICWaR), Indian Institute of Science (IISc)


"""The purpose of this script is to,
    load longterm mean for our GRACE SH data

    For this, we need to input the GRACE data source as well as the path to the longterm mean values
    Data source may be CSR, JPL or ITSG.

    For RL06, example data have been provided within the package. In case this option is chosen, the program directly returns the longterm mean values.

    Returns the long_mean path
"""
import pkg_resources

def load_longterm_mean(source = "", use_sample_mean = 0):
    """_summary_

    Args:
        source (str, optional): _description_. Defaults to "".
        use_sample_mean (int, optional): _description_. Defaults to 0.

    Raises:
        Exception: _description_

    Returns:
        _type_: _description_
    
    Todo:
        + Not sure if using "source = ''" is all right
        + instead of base eception is can be ValueError
    """
    if use_sample_mean == 1:
        print("Loading preloaded RL06 long term mean values")
        print("Please ensure that your data is RL06 \nIf not, please manually input long term mean by setting use_sample_mean = 0")

        if str.upper(source) == 'CSR':
            long_mean = pkg_resources.resource_filename('pyshbundle', 'data/RL06_long_mean/SH_long_mean_csr.npy')
        elif str.upper(source) == 'JPL':
            long_mean = pkg_resources.resource_filename('pyshbundle', 'data/RL06_long_mean/SH_long_mean_itsg.npy')
        elif str.upper(source) == 'ITSG':
            long_mean = pkg_resources.resource_filename('pyshbundle', 'data/RL06_long_mean/SH_long_mean_jpl.npy')
        else:
            raise Exception("Incorrect selection of source")
        print("Successfully loaded preloaded longterm means")
    else:
        print("Please download and provide the longterm GRACE SH mean values")
        print("Instructions to download the longterm GRACE SH mean values may be referred to in https://github.com/mn5hk/pyshbundle/blob/main/docs/index.md#how-to-download-data")
        long_mean = str(input("Enter the longterm mean for the SH values in the numpy (.npy) format"))
        print("Successfully loaded path to long term mean:", long_mean)

    return long_mean