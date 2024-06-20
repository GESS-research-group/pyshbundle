# Applications of GRACE data to Hydrology
# Curator: Abhishek Mhamane

import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import geopandas as gpd
import rioxarray
from shapely.geometry import mapping
from pyshbundle.shutils import Gaussian
from pyshbundle.pysh_core import gshs


def TWSCalc(data, lmax: int, gs: float, r, m):
    """_summary_

    Args:
        data (np.ndarray): SC coefficients
        lmax (int): maximum degree
        gs (float): grid size
        r (_type_): _description_
        m (_type_): _description_
    
    Author:
        Vivek Yadav, Interdisciplinary Center for Water Research (ICWaR), Indian Institute of Science (IISc)
    """
    SC = data
    
    gfilter = Gaussian(lmax,r)
    grid_y = int(180/gs)
    grid_x = int(360/gs)
    tws_f = np.zeros([m,grid_y,grid_x], dtype ='longdouble')
    for i in tqdm(range(0,m,1)):
        field = SC[i,0:lmax+1,96-lmax:96+lmax+1]
        shfil = np.zeros([lmax+1,2*lmax+1])

        for j in range(0,2*lmax+1,1):
            shfil[:,j] = gfilter[:,0] * field[:,j]
        
        quant = 'water' 
        grd = 'cell'
        n = int(180/gs) 
        h = 0 
        jflag = 0
        
        ff = gshs(shfil, quant, grd, n, h, jflag)[0]
        
        ff = ff*1000    # convert units from m to mm
        tws_f[i,:,0:int(grid_x/2)] = ff[:,int(grid_x/2):]
        tws_f[i,:,int(grid_x/2):] = ff[:,0:int(grid_x/2)]   
    
    plt.imshow(tws_f[0,:,:])
    return(tws_f)

def apply_gaussian(sc_coeff, gaussian_coeff, lmax):
    
    # filtered SH Coeff
    shfil = np.zeros([lmax+1, 2 * lmax+1])

    # applying filter on substracted coeff
    for j in range(0, 2*lmax+1, 1):
        shfil[:,j] = gaussian_coeff[:,0] * sc_coeff[:,j]
    
    return shfil

def area_weighting(grid_resolution):
    """
     Calculate the area of each grid corresponding to the latitudes and longitudes.
 
     Parameters:
         grid_resolution (float): The resolution of the grid in grid_resolutionrees.
 
     Returns:
         numpy.ndarray: The area of each grid in square meters.
     """
    longitude_grid = np.linspace(0, 359+(1-grid_resolution), int(360/grid_resolution), dtype='double');
    latitude_grid = np.linspace(0, 179+(1-grid_resolution), int(180/grid_resolution), dtype='double');
    longitude_grid1 = np.linspace(1*grid_resolution, 360, int(360/grid_resolution), dtype='double');
    latitude_grid1 = np.linspace(1*grid_resolution, 180, int(180/grid_resolution), dtype='double');
    
    lambd,theta = np.meshgrid(longitude_grid,latitude_grid)  
    lambd1,theta1 = np.meshgrid(longitude_grid1,latitude_grid1)
    
    delta_latitude = np.sin(np.deg2rad(90-theta))-np.sin(np.deg2rad(90-theta1))
    delta_longitude = (lambd1 - lambd)*np.pi/180

    # Area of each grid
    area = (6378.137**2)*pow(10,6)*(np.multiply(delta_latitude,delta_longitude)) # units m^2
    return area

def Basinaverage(temp, gs, shp_basin, basin_area):

    from pyshbundle.hydro import area_weighting
    # area_weighting returns the area of each grid in m^2 for the grid resolution specified
    temp_weighted=temp.copy()
    temp_weighted['tws']=temp['tws']*area_weighting(gs)

    ''' add projection system to nc '''
    basin_tws = temp_weighted.rio.write_crs("EPSG:4326", inplace=True)
    basin_tws = basin_tws.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
            
    # mask data with shapefile
    basin_tws = basin_tws.rio.clip(shp_basin.geometry.apply(mapping), shp_basin.crs,drop=True)
    basin_avg_tws=basin_tws.sum(dim=('lon','lat'), skipna=True)/basin_area  #basin average tws

    basin_tws=basin_tws.drop_vars('spatial_ref')
    basin_avg_tws=basin_avg_tws.drop_vars('spatial_ref')

    return basin_tws, basin_avg_tws