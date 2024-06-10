''' Area weighting '''
import numpy as np

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