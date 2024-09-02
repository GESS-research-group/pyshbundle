#!/usr/bin/env python
# coding: utf-8

# Unit test script file for pyshbundle package
# This scripts follows the notebook and `04_TWS_time_series.ipynb` `validation_pyshbundle.ipynb` and is used to validate the 
# pyshbundle package. Please read them for more details.
import numpy as np
import pandas as pd
import xarray as xr
import os
import scipy.io
import sys
from datetime import datetime
from collections import OrderedDict
from pyshbundle.hydro import TWSCalc
from pyshbundle.io import extract_SH_data, extract_deg1_coeff_tn13, extract_deg2_3_coeff_tn14
import pkg_resources
long_mean_file_path = pkg_resources.resource_filename('pyshbundle', 'data/long_mean/SH_long_mean_jpl.npy')
matlab_file_path = pkg_resources.resource_filename('pyshbundle', 'data/validation_data/tws_sh.mat')
# long_mean_file_path = '../pyshbundle/data/long_mean/SH_long_mean_jpl.npy'
# matlab_file_path = '../examples/validation_data/tws_sh.mat'

ignore_warnings = True

# Add the folder path to the Python path
# Get the current directory
current_dir = os.path.dirname(os.path.abspath(__file__))

# Add the parent directory to the Python path
parent_dir = os.path.dirname(current_dir)
sys.path.append(parent_dir)

# Add the necessary paths
sys.path.append(os.path.join(parent_dir, 'pyshbundle'))
sys.path.append(os.path.join(parent_dir, 'sample_input_data'))

# Rest of the code...
def validation_pyshbundle():
    source='jpl'
    path_sh = os.path.join(parent_dir, 'data', 'JPL_input')

    path_tn14 = os.path.join(parent_dir, 'pyshbundle','data', 'JPL_TN_files', 'TN-14_C30_C20_GSFC_SLR.txt')
    path_tn13 = os.path.join(parent_dir, 'pyshbundle','data', 'JPL_TN_files', 'TN-13_GEOC_JPL_RL06.txt')
    files = os.listdir(path_sh)
    file_paths = [os.path.join(path_sh, file) for file in files if os.path.splitext(file)[1] == '.gz'];
    extracted_data={} 
    for file_path in file_paths:
        # file_data = read_sh(file_path, source=source)
        file_data = extract_SH_data(file_path, source=source)
        if file_data['time_coverage_start']:
            # Convert time_coverage_start to a datetime object and then format it as yyyy-mm
            if source == 'itsg':
                start_date = datetime.strptime(file_data['time_coverage_start'][-7:], '%Y-%m').strftime('%Y-%m')
            else:
                start_date = datetime.strptime(file_data['time_coverage_start'], '%Y-%m-%dT%H:%M:%S.%f').strftime('%Y-%m')
            # Use the formatted date as the key
            extracted_data[start_date] = file_data['coefficients']


    # Time Sort the dictionary by keys (dates)
    sorted_data = OrderedDict(sorted(extracted_data.items()));
    temp_tn14 = extract_deg2_3_coeff_tn14(path_tn14)
    for date_key in temp_tn14.keys():
        if date_key in sorted_data.keys():
            sorted_data[date_key][(2,0)]['Clm'] = temp_tn14[date_key]['c20']
            if temp_tn14[date_key]['c30'] is not None:
                sorted_data[date_key][(3,0)]['Clm'] = temp_tn14[date_key]['c30'];
    temp_tn13=extract_deg1_coeff_tn13(path_tn13)
    if str.upper(source)=='JPL':
        for key in sorted_data:
            sorted_data[key][(0, 0)] = {'Clm': 0.0, 'Slm': 0.0, 'Clm_sdev': 0.0, 'Slm_sdev': 0.0}
            sorted_data[key][(1, 0)] = {'Clm': 0.0, 'Slm': 0.0, 'Clm_sdev': 0.0, 'Slm_sdev': 0.0}
            sorted_data[key][(1, 1)] = {'Clm': 0.0, 'Slm': 0.0, 'Clm_sdev': 0.0, 'Slm_sdev': 0.0};
    else: pass
    for date_key in temp_tn13.keys():
        if date_key[0] in sorted_data.keys():
            sorted_data[date_key[0]][(date_key[1], date_key[2])]['Clm'] = temp_tn13[date_key]['Clm']
            sorted_data[date_key[0]][(date_key[1], date_key[2])]['Slm'] = temp_tn13[date_key]['Slm']
            sorted_data[date_key[0]][(date_key[1], date_key[2])]['Clm_sdev'] = temp_tn13[date_key]['Clm_sdev']
            sorted_data[date_key[0]][(date_key[1], date_key[2])]['Slm_sdev'] = temp_tn13[date_key]['Slm_sdev'];
    max_degree=np.max([degree for date in sorted_data.keys() for degree, order in sorted_data[date].keys()])
    max_order=np.max([order for date in sorted_data.keys() for degree, order in sorted_data[date].keys()])
    number_of_months=len(sorted_data.keys())
    sc_mat=np.zeros([number_of_months, max_degree+1, 2*(max_degree+1)], dtype=np.double)

    for index, key in enumerate(sorted_data.keys()):
        temp=sorted_data[key]
        for l in range(0,max_degree+1):
            for m in range(0,l+1):
                sc_mat[index, l, max_degree+m+1]=temp[(l,m)]['Clm']
                sc_mat[index, l, max_degree-m]=temp[(l,m)]['Slm']
        del temp
    sc_mat=np.delete(sc_mat, max_degree, axis=2);
    SH_long_mean_jpl = np.load(long_mean_file_path)    
    delta_sc=np.ones_like(sc_mat)*np.nan
    delta_sc = sc_mat -   SH_long_mean_jpl
    lmax,gs,half_rad_gf=96, 1, 500
    tws_fields = TWSCalc(delta_sc,lmax, gs,half_rad_gf, number_of_months)
    tws_fields = np.float32(tws_fields)
    lon = np.arange(-180,180,gs)
    lat = np.arange(89,-91,-gs)
    dates = pd.to_datetime(list(sorted_data.keys()), format='%Y-%m',) \
                + pd.offsets.MonthEnd(0)

    ds = xr.Dataset(
        data_vars=dict(
            tws=(["time","lat", "lon"], tws_fields)
        ),
        coords = {
            # "time":(('time'),dates),
            "time":dates,
            "lat":lat,
            "lon":lon },);

    # Load the .mat file
    data = scipy.io.loadmat(matlab_file_path)
    # Access the variables in the .mat file
    var1 = data['tws_m']
    ds_pysh = ds.copy()

    # #### Lets convert the shbundle datasets to a netcdf format for easier calculations.
    # * Converting `shbundle` processed data into netcdf format using xarray, to `ds_pysh`
    ds_msh = xr.Dataset(
        data_vars=dict(
            tws=(["time","lat", "lon"], var1)
        ),
        coords = {
            "time":(('time'),dates),
            "lat":lat,
            "lon":lon },);
    # print(ds_msh)
    # print(ds_pysh)
    # ## 1. Gridwise RMSE calculation
    # 
    # Before finding the grid-wise RMSE values we need to ignore the data for the missing months.
    # Calculate the difference between the two datasets
    diff = ds_msh['tws'].dropna(dim='time').values - ds_pysh['tws'].dropna(dim='time').values    # dropna is used to remove nan values, the dates where the GRACE data is missing

    # Calculate the squared difference
    squared_diff = diff**2

    # Calculate the mean squared difference along the time axis
    mean_squared_diff = np.mean(squared_diff, axis=0)

    # Calculate the root mean squared error (RMSE)
    gridwise_rmse = np.sqrt(mean_squared_diff)

    # Test whether gridwise RMSE is less than 1e-3
    if np.all(gridwise_rmse < 1e-3):
        pass
    else:
        raise ValueError('Gridwise RMSE is greater than 1e-3')

    # ## 2. Gridwise NRMSE
    # Calculate the normalized root mean squared error (NRMSE)
    gridwise_nrmse = gridwise_rmse/np.std(ds_msh['tws'].dropna(dim='time').values, axis=0)
    
    # Test whether gridwise NRMSE is less than 1e-5
    if np.all(gridwise_nrmse < 1e-5):
        pass
    else:
        raise ValueError('Gridwise NRMSE is greater than 1e-4')
    return "expected_result"