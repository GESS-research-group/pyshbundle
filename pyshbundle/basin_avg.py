#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 22:37:19 2022

@author: wslvivek
"""
def basin_avg(data,path,c_rs,m,gs):
    import xarray as xr
    from matplotlib import pyplot as plt
    import numpy as np
    import salem
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
    deg = gs;
    x = np.linspace(0, 359+(1-deg), int(360/deg), dtype='double');
    y = np.linspace(0, 179+(1-deg), int(180/deg), dtype='double');
    x1 = np.linspace(1*deg, 360, int(360/deg), dtype='double');
    y1 = np.linspace(1*deg, 180, int(180/deg), dtype='double');
    lambd,theta = np.meshgrid(x,y)  
    lambd1,theta1 = np.meshgrid(x1,y1)  
    a = np.sin(np.deg2rad(90-theta))-np.sin(np.deg2rad(90-theta1))
    b = (lambd1 - lambd)*np.pi/180
    
    
    # Area of each grid (360*720)
    area = (6378.137**2)*pow(10,6)*(np.multiply(a,b))        # units m^2
    tot_area = np.sum(np.sum(area))
    tws_m = np.zeros([m,lat_shape,lon_shape])
    for i in range(0,m,1):
        tws_m[i,:,:] = np.multiply(tws_val[i,:,:],area)
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
    fig,ax=plt.subplots(figsize=(15,5))
    alpha.plot(ax=ax,color='b');
    ax.set_box_aspect(0.33)
    ax.set_title('Time series for the basin',size=15)
    ax.set_ylabel('TWS anomaly in mm ',size=15)
    plt.tight_layout()
    
    return alpha