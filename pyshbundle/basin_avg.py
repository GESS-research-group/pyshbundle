#!/usr/bin/env python3
#-*- coding: utf-8 -*-

# Created on Wed Dec 14 22:37:19 2022

# @author: Vivek Yadav, Interdisciplinary Center for Water Research (ICWaR), Indian Institute of Science (IISc)

# ySHbundle: A Python implementation of MATLAB codes SHbundle
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
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

import xarray as xr
from matplotlib import pyplot as plt
import numpy as np
import salem

def basin_avg(data, path: str, c_rs, m, gs):
    """Computes the TWSA time-series for a given basin shape file, using the SH data.

    Args:
        data (xarray.Dataset): xarray dataset with format - {coordinates: [time, lat, lon], Data variables: [tws]}
        path (str): valid path to the basin shape file with extension (.shp)
        c_rs (crs): the crs into which the dataframe must be transformed (related to salem module)
        m (int): number of files read
        gs (float): grid size 

    Returns:
        xarray.Dataset: basin averaged values of TWS
        _type_: _description_
    """
    shdf = salem.read_shapefile(path)
    shdf.crs
    shdf.plot()
    shdf_area = sum(shdf.to_crs(c_rs).area)
    print('Area of basin in km2:',shdf_area/1e6)
    if shdf_area < 63*1e9:
        print('Warning basin too small for GRACE data application')
    
    tws_val = data.tws.values
    dates = data.time
    lat,lon = data.lat, data.lon
    lat_shape, lon_shape = data.tws.shape[1],data.tws.shape[2]
    
    # Calculation of area of each corresponding to  the latitudes and longitudes
    # not sure if ';' is proper syntax may be the octave residu

    deg = gs
    x = np.linspace(0, 359+(1-deg), int(360/deg), dtype='double')
    y = np.linspace(0, 179+(1-deg), int(180/deg), dtype='double')
    x1 = np.linspace(1*deg, 360, int(360/deg), dtype='double')
    y1 = np.linspace(1*deg, 180, int(180/deg), dtype='double')
    lambd,theta = np.meshgrid(x,y)
    lambd1,theta1 = np.meshgrid(x1,y1)
    a = np.sin(np.deg2rad(90-theta))-np.sin(np.deg2rad(90-theta1))
    b = (lambd1 - lambd)*np.pi/180
    
    
    # Area of each grid (360*720)
    area = (6378.137**2)*pow(10, 6)*(np.multiply(a, b))        # units m^2
    tot_area = np.sum(np.sum(area))
    tws_m = np.zeros([m, lat_shape, lon_shape])
    for i in range(0,m,1):
        tws_m[i, :, :] = np.multiply(tws_val[i, :, :],area)
    ds_area_w = xr.Dataset(
    data_vars=dict(
        tws=(["time","lat", "lon"], tws_m)
    ),
    coords = {
        "time":dates,
        "lat":lat,
        "lon":lon },
    attrs=dict(description="TWS Anomaly corresponding to long term (2004-2010) mean \n lmax=96 and half radius of Gaussian filter = 500Km"),
    )
    
    ds_area_w_clp= ds_area_w.salem.roi(shape=shdf)
    # Time series for the whole basin(shapefile) in user defined range
    alpha = ds_area_w_clp.tws.sum(dim=('lon','lat'))/shdf_area
    fig,ax = plt.subplots(figsize=(15,5))
    alpha.plot(ax=ax, color='b')
    ax.set_box_aspect(0.33)
    ax.set_title('Time series for the basin', size=15)
    ax.set_ylabel('TWS anomaly in mm ', size=15)
    plt.tight_layout()
    
    return alpha, ds_area_w