"""
Created on Fri Feb  17 2023
@author: Amin Shakya

The purpose of this script is to,
    load longterm mean for our GRACE SH data

    For this, we need to input the GRACE data source as well as the path to the longterm mean values
    Data source may be CSR, JPL or ITSG.

    For RL06, example data have been provided within the package. In case this option is chosen, the program directly returns the longterm mean values.

    Returns the long_mean path
"""

def load_longterm_mean(source = "", use_sample_mean = 0):
    if use_sample_mean == 1:
        import pkg_resources

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