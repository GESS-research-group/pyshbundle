# Applications of GRACE data to Hydrology
# Curator: Abhishek Mhamane
# Updated: Vivek, 2024-05-20, added the function Basinaverage, area_weighting.
# Updated: Vivek, 2024-07-28.

# - - - - - - - - - - - - - - 
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

import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import geopandas as gpd
import rioxarray
from shapely.geometry import mapping
from pyshbundle.shutils import Gaussian
from pyshbundle.pysh_core import gshs


def TWSCalc(data, lmax: int, gs: float, r:float, m: int):
    """
    Spherical Harmonics Synthesis for Total Water Storage (TWS) calculation.

    Calculate the total water storage (TWS) from spherical harmonics coefficients.
    Uses spherical harmonics synthesis to go from harmonics to gridded domain.

    Args:
        data (numpy.ndarray): Spherical harmonics coefficients in SC format.
        lmax (int): Maximum degree of the spherical harmonics coefficients.
        gs (float): Grid size.
        r (float): Half-width of the Gaussian filter.
        m (int): Number of time steps.

    Returns:
        (numpy.ndarray): Gridded TWS data.

    Uses:
        'Gaussian', 'gshs'

    Author:
        Vivek Yadav, Interdisciplinary Center for Water Research (ICWaR), Indian Institute of Science (IISc)
    """
    SC = data
    
    gfilter = Gaussian(lmax,r)
    grid_y = int(180/gs)
    grid_x = int(360/gs)
    tws_f = np.zeros([m,grid_y,grid_x], dtype ='float')
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
    """
    Apply Gaussian filter on the spherical harmonics coefficients.

    Args:
        sc_coeff (numpy.ndarray): Spherical harmonics coefficients in SC format.
        gaussian_coeff (numpy.ndarray): Gaussian filter weights.
        lmax (int): Maximum degree of the spherical harmonics coefficients.

    Returns:
        (numpy.ndarray): Filtered spherical harmonics coefficients in SC format.
    """
    
    # filtered SH Coeff
    shfil = np.zeros([lmax+1, 2 * lmax+1])

    # applying filter on substracted coeff
    for j in range(0, 2*lmax+1, 1):
        shfil[:,j] = gaussian_coeff[:,0] * sc_coeff[:,j]
    
    return shfil

def area_weighting(grid_resolution):
    """
    Calculate the area of each grid, globally, corresponding to the latitudes and longitudes.

    Args:
        grid_resolution (float): The resolution of the grid in degrees.

    Returns:
        (numpy.ndarray): A matrirx of area of each grid on the surface of the Earth in m^2.
    """
    
    longitude_grid = np.linspace(0, 359+(1-grid_resolution), int(360/grid_resolution), dtype='float');
    latitude_grid = np.linspace(0, 179+(1-grid_resolution), int(180/grid_resolution), dtype='float');
    longitude_grid1 = np.linspace(1*grid_resolution, 360, int(360/grid_resolution), dtype='float');
    latitude_grid1 = np.linspace(1*grid_resolution, 180, int(180/grid_resolution), dtype='float');
    
    lambd,theta = np.meshgrid(longitude_grid,latitude_grid)  
    lambd1,theta1 = np.meshgrid(longitude_grid1,latitude_grid1)
    
    delta_latitude = np.sin(np.deg2rad(90-theta))-np.sin(np.deg2rad(90-theta1))
    delta_longitude = (lambd1 - lambd)*np.pi/180

    # Area of each grid
    area = (6378.137**2)*pow(10,6)*(np.multiply(delta_latitude,delta_longitude)) # units m^2
    return area

def Basinaverage(temp, gs, shp_basin, basin_area):
    """
    Calculate the basin average of the total water storage (TWS) from the gridded TWS data.

    Applies area weighting to the gridded TWS data and then clips the data to the basin shapefile.
    Followed by summation of data over the latitude and longitude dimensions, divides it by basin
    area to get the basin average TWS.

    Args:
        temp (xarray.DataArray): Gridded total water storage data.
        gs (float): Grid size.
        shp_basin (geopandas.GeoDataFrame): Shapefile of the basin.
        basin_area (float): Area of the basin in square meters.

    Returns:
        basin_tws (xarray.DataArray): Total water storage data clipped to the basin.
        basin_avg_tws (xarray.DataArray): Basin average total water storage.
    """

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