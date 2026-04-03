import pkg_resources

def read_GRACE_SH_paths(use_sample_files = 0):
    """
    Returns path of data files, path of tn13 and path of tn14 replacement files.

    Args:
        use_sample_files (int, optional): Defaults to 0.

    Raises:
        Exception: If the source selection is incorrect.

    Returns:
        tuple: 
            - path_sh (str): Path of data files.
            - path_tn13 (str): Path of tn13 replacement file.
            - path_tn14 (str): Path of tn14 replacement file.
            - source (str): Source of the SH files (JPL, ITSG, or CSR).

    Remarks:
        The purpose of this script is to:
        - Read what the data source is (JPL, CSR, or ITSG).
        - Read file path for GRACE L2 spherical harmonics inputs.
        - Read replacement files for tn13 and tn14.
        - Identify the source of the SH files (JPL, ITSG, or CSR).
    """
    #Created on Fri Feb  17 2023
    #@author: Amin Shakya, Interdisciplinary Center for Water Research (ICWaR), Indian Institute of Science (IISc)
    
    print("This program supports working with GRACE L2 Spherical harmonics data from the following centers: CSR, JPL and ITSG")
    print("Instructions to download data may be referred to in https://github.com/GESS-research-group/pyshbundle/blob/main/docs/index.md#how-to-download-data")
    source = str(input("Enter the source of L2 SH coeffs code(jpl, csr, gfz): "))

    if use_sample_files ==1:
        
        print("You have chosen to use sample replacement files.")
        print("The replacement files for the TN13 and TN14 Args have been preloaded into the program")
        print("Due to the size of the GRACE SH files, these have not been preloaded into the program")
        print("You may download the GRACE SH L2 files from the link below. Please ensure to download the files as per your selection of source in the prior step")
        print("Download sample files from: https://github.com/GESS-research-group/pyshbundle/tree/main/sample_input_data")
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
            # path_tn13 = as_file(files('pyshbundle', 'data/sample_CSR_TN_files/TN-14_C30_C20_SLR_GSFC.txt'))
            # path_tn14 = as_file(files('pyshbundle', 'data/sample_CSR_TN_files/TN-13_GEOC_CSR_RL06.1.txt'))
            print("Successfully loaded preloaded TN13 and TN14 replacement files for CSR")
        else:
            path_tn13 = str(input("Enter the path to the file for tn13 replacement in .txt format"))
            path_tn14 = str(input("Enter the path to the file for tn14 replacement in .txt format"))
            print("Successfully loaded TN13 and TN14 replacement files for CSR")

    elif str.upper(source) == 'ITSG':
        if use_sample_files == 1:
            path_tn13 = pkg_resources.resource_filename('pyshbundle', 'data/sample_ITSG_TN_files/TN-13_GEOC_CSR_RL06.1.txt')
            path_tn14 = pkg_resources.resource_filename('pyshbundle', 'data/sample_ITSG_TN_files/TN-14_C30_C20_SLR_GSFC.txt')
            # path_tn13 = as_file(files('pyshbundle', 'data/sample_ITSG_TN_files/TN-13_GEOC_CSR_RL06.1.txt'))
            # path_tn14 = as_file(files('pyshbundle', 'data/sample_ITSG_TN_files/TN-14_C30_C20_SLR_GSFC.txt'))
            print("Successfully loaded preloaded TN13 and TN14 replacement files for ITSG")
        else:
            path_tn13 = str(input("Enter the path to the file for tn13 replacement in .txt format"))
            path_tn14 = str(input("Enter the path to the file for tn14 replacement in .txt format"))
            print("Successfully loaded TN13 and TN14 replacement files for ITSG")
    else:
        raise Exception("Source selection is incorrect. Please select between JPL, CSR or gfz")

    return path_sh, path_tn13, path_tn14, source



def load_longterm_mean(source = "", use_sample_mean = 0):
    """
    Loads the long term mean values for the GRACE SH data.

    Args:
        source (str, optional): Source of data. Defaults to "".
        use_sample_mean (int, optional): Whether to use default long-mean values provided with the data. Defaults to 0.

    Raises:
        ValueError: If the source selection is incorrect.

    Returns:
        str: Path of the appropriate long term mean file.

    Todo:
        + Not sure if using "source = ''" is all right.
        + Instead of base exception, it can be ValueError.
    """
# @author: Amin Shakya, Interdisciplinary Center for Water Research (ICWaR), Indian Institute of Science (IISc)

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
        print("Instructions to download the longterm GRACE SH mean values may be referred to in https://github.com/GESS-research-group/pyshbundle/blob/main/docs/index.md#how-to-download-data")
        long_mean = str(input("Enter the longterm mean for the SH values in the numpy (.npy) format"))
        print("Successfully loaded path to long term mean:", long_mean)

    return long_mean