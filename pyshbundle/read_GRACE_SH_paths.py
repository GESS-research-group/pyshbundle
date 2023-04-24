#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Created on Fri Feb  17 2023
#@author: Amin Shakya, Interdisciplinary Center for Water Research (ICWaR), Indian Institute of Science (IISc)

#The purpose of this script is to,
#    firstly read what the data source is (JPL, CSR or ITSG)
#    read file path for GRACE L2 spherical harmonics inputs,
#    read replacement files for tn13 and tn14
#    source of the SH files (JPL, ITSG or CSR)
# The code returns path of data files, path of tn13 and path of tn14 replacement files

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
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

"""The purpose of this script is to,
    firstly read what the data source is (JPL, CSR or ITSG)
    read file path for GRACE L2 spherical harmonics inputs,
    read replacement files for tn13 and tn14
    source of the SH files (JPL, ITSG or CSR)
"""
import pkg_resources

def read_GRACE_SH_paths(use_sample_files = 0):
    """Returns path of data files, path of tn13 and path of tn14 replacement files

    Args:
        use_sample_files (int, optional): _description_. Defaults to 0.

    Raises:
        Exception: _description_

    Returns:
        _type_: _description_
    """
    
    print("This program supports working with GRACE L2 Spherical harmonics data from the following centers: CSR, JPL and ITSG")
    print("Instructions to download data may be referred to in https://github.com/mn5hk/pyshbundle/blob/main/docs/index.md#how-to-download-data")
    source = str(input("Enter the source of L2 SH coeffs code(jpl, csr, gfz): "))

    if use_sample_files ==1:
        
        print("You have chosen to use sample replacement files.")
        print("The replacement files for the TN13 and TN14 parameters have been preloaded into the program")
        print("Due to the size of the GRACE SH files, these have not been preloaded into the program")
        print("You may download the GRACE SH L2 files from the link below. Please ensure to download the files as per your selection of source in the prior step")
        print("Download sample files from: https://github.com/mn5hk/pyshbundle/tree/main/sample_input_data")
    path_sh = str(input("Enter the path to the folder with SH L2 data"))

    
    if str.upper(source) == 'JPL':
        if use_sample_files == 1:
            path_tn13 = pkg_resources.resource_filename('pyshbundle', 'data/sample_JPL_TN_files/TN-13_GEOC_JPL_RL06.txt')
            path_tn14 = pkg_resources.resource_filename('pyshbundle', 'data/sample_JPL_TN_files/TN-14_C30_C20_GSFC_SLR.txt')
            print("Successfully loaded preloaded TN13 and TN14 replacement files for JPL")
        else:
            path_tn13 = str(input("Enter the path to the file for tn13 replacement in .txt format"))
            path_tn14 = str(input("Enter the path to the file for tn14 replacement in .txt format"))
            print("Successfully loaded TN13 and TN14 replacement files for JPL")

    elif str.upper(source) == 'CSR':
        if use_sample_files == 1:
            path_tn13 = pkg_resources.resource_filename('pyshbundle', 'data/sample_CSR_TN_files/TN-14_C30_C20_SLR_GSFC.txt')
            path_tn14 = pkg_resources.resource_filename('pyshbundle', 'data/sample_CSR_TN_files/TN-13_GEOC_CSR_RL06.1.txt')
            print("Successfully loaded preloaded TN13 and TN14 replacement files for CSR")
        else:
            path_tn13 = str(input("Enter the path to the file for tn13 replacement in .txt format"))
            path_tn14 = str(input("Enter the path to the file for tn14 replacement in .txt format"))
            print("Successfully loaded TN13 and TN14 replacement files for CSR")

    elif str.upper(source) == 'ITSG':
        if use_sample_files == 1:
            path_tn13 = pkg_resources.resource_filename('pyshbundle', 'data/sample_ITSG_TN_files/TN-13_GEOC_CSR_RL06.1.txt')
            path_tn14 = pkg_resources.resource_filename('pyshbundle', 'data/sample_ITSG_TN_files/TN-14_C30_C20_SLR_GSFC.txt')
            print("Successfully loaded preloaded TN13 and TN14 replacement files for ITSG")
        else:
            path_tn13 = str(input("Enter the path to the file for tn13 replacement in .txt format"))
            path_tn14 = str(input("Enter the path to the file for tn14 replacement in .txt format"))
            print("Successfully loaded TN13 and TN14 replacement files for ITSG")
    else:
        raise Exception("Source selection is incorrect. Please select between JPL, CSR or gfz")

    return path_sh, path_tn13, path_tn14, source