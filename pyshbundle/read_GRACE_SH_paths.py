"""
Created on Fri Feb  17 2023
@author: Amin Shakya

The purpose of this script is to,
    firstly read what the data source is (JPL, CSR or ITSG)
    read file path for GRACE L2 spherical harmonics inputs,
    read replacement files for tn13 and tn14
    source of the SH files (JPL, ITSG or CSR)
The code returns path of data files, path of tn13 and path of tn14 replacement files
"""

def read_GRACE_SH_paths(use_sample_files = 0):
    
    print("This program supports working with GRACE L2 Spherical harmonics data from the following centers: CSR, JPL and ITSG")
    print("Instructions to download data may be referred to in https://github.com/mn5hk/pyshbundle/blob/main/docs/index.md#how-to-download-data")
    source = str(input("Enter the source of L2 SH coeffs code(jpl, csr, gfz): "))

    if use_sample_files ==1:
        import pkg_resources
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