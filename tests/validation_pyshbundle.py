#!/usr/bin/env python
# coding: utf-8

# Unit test script file for pyshbundle package
# This scripts follows the notebook and `04_TWS_time_series.ipynb` `validation_pyshbundle.ipynb` and is used to validate the 
# pyshbundle package. Please read them for more details.
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import os
import scipy.io

def validation_pyshbundle():
    # Load the .mat file
    data = scipy.io.loadmat('../pyshbundle/validation_data/tws_m.mat')
    # Access the variables in the .mat file
    var1 = data['tws_m']

    temp=xr.open_dataset('../pyshbundle/validation_data/tws_py.nc', engine="netcdf4")
    var2=temp.tws.values


    # #### Lets convert both datasets to a netcdf format for easier calculations.

    # * Converting `pyshbundle` processed data into netcdf format using xarray, to `ds_pysh`
    gs=1;
    lon = np.arange(-180,180,gs)
    lat = np.arange(89,-91,-gs)
    ds_pysh = xr.Dataset(
        data_vars=dict(
            tws=(["time","lat", "lon"], var2)
        ),
        coords = {
            "time":(('time'),temp.time.data),
            "lat":lat,
            "lon":lon },);
    ds_msh = xr.Dataset(
        data_vars=dict(
            tws=(["time","lat", "lon"], var1)
        ),
        coords = {
            "time":(('time'),temp.time.data),
            "lat":lat,
            "lon":lon },);

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

    # ## 2. Gridwise NRMSE
    # Calculate the normalized root mean squared error (NRMSE)
    gridwise_nrmse = gridwise_rmse/np.std(ds_msh['tws'].dropna(dim='time').values, axis=0)

    # ## 3. Global area weighted water budget closure
    # Area of grids
    from pyshbundle.hydro import area_weighting
    global_grid_area=area_weighting(1)
    global_grid_area_sum = np.sum(global_grid_area)
    print('global surface area in m\u00b2:', global_grid_area_sum)


    # * Calculate the global area weighted water budget closure error
    # # Create a copy of the datasets
    ds_msh_area_weighted, ds_pysh_area_weighted = ds_msh.copy(), ds_pysh.copy()

    # Area weight with the global grid area and calculate the sum over lat and lon
    ds_msh_area_weighted = ds_msh['tws']*global_grid_area / global_grid_area_sum
    ds_msh_area_weighted = ds_msh_area_weighted.sum(dim=['lat', 'lon'])

    # Same for the pyshbundle dataset
    ds_pysh_area_weighted = ds_pysh['tws']*global_grid_area / global_grid_area_sum
    ds_pysh_area_weighted = ds_pysh_area_weighted.sum(dim=['lat', 'lon'])

    diff_global = ds_msh_area_weighted - ds_pysh_area_weighted

    # Reinsert the NaN values where the GRACE data is missing and the time coordinate
    ds_msh_area_weighted=ds_msh_area_weighted.where(~np.isnan(ds_pysh['tws'][:,0,0]), np.nan)
    ds_pysh_area_weighted=ds_pysh_area_weighted.where(~np.isnan(ds_pysh['tws'][:,0,0]), np.nan)
    diff_global=diff_global.where(~np.isnan(ds_pysh['tws'][:,0,0]), np.nan)

    # ## 4. Difference in basin-average Time Series
    import geopandas as gpd
    path_shapefile = '/Users/vivek/Desktop/vivek_desktop/mrb_shp_zip/mrb_basins.shp'
    shp = gpd.read_file(path_shapefile)
    shp.plot(figsize=(8, 4))  
    basin_name='KRISHNA'
    shp_basin=shp[shp['RIVER_BASI']==basin_name]
    print(shp_basin.head(), '\n')
    shp_basin.plot()
    basin_area=np.float64(shp_basin['Shape_Area'].values)*1e12          # basin area already in m^2
    print('Basin area is :', basin_area, 'm\u00b2');

    from pyshbundle.hydro import Basinaverage
    _, basin_avg_tws_msh = Basinaverage(ds_msh, gs, shp_basin, basin_area)
    _, basin_avg_tws_pysh = Basinaverage(ds_pysh, gs, shp_basin, basin_area);

    new_dates=pd.date_range(start=basin_avg_tws_msh.time[0].values, 
                            end=basin_avg_tws_msh.time[-1].values, freq='ME',);

    # Empty dataset for the gapped data, msh
    basin_avg_tws_gapped_msh = xr.Dataset(
            data_vars = dict(   tws=(["time"], np.nan*np.arange(len(new_dates)))),
            coords=dict(time=new_dates),);
    # Empty dataset for the gapped data, pysh
    basin_avg_tws_gapped_pysh = xr.Dataset(
            data_vars = dict(   tws=(["time"], np.nan*np.arange(len(new_dates)))),
            coords=dict(time=new_dates),);

    #
    basin_avg_tws_gapped_msh['tws'] = basin_avg_tws_msh['tws'].where(
        basin_avg_tws_msh['time'].isin(basin_avg_tws_gapped_msh['time']),)
    #
    basin_avg_tws_gapped_pysh['tws'] = basin_avg_tws_pysh['tws'].where(
        basin_avg_tws_pysh['time'].isin(basin_avg_tws_gapped_pysh['time']),)

    diff_global_basin = basin_avg_tws_gapped_msh['tws'] - basin_avg_tws_gapped_pysh['tws']