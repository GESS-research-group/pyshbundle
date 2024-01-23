# Applications of GRACE data to Hydrology
# Curator: Abhishek Mhamane

import salem
import numpy as np
import xarray as xr
from scipy import signal
import matplotlib.pyplot as plt
from tqdm import tqdm

from pyshbundle.shutils import Gaussian
from pyshbundle.core_process import GSHS

def PhaseCalc(fts, ffts):
    """calculates the phase difference between two time series based on the
    Hilbert transform method explained by Phillip et al.

    Args:
        fts (np.ndarray): _description_
        ffts (np.ndarray): _description_

    Returns:
        _type_: _description_
    
    References:
        1. Phillips, T., R. S. Nerem, B. Fox-Kemper, J. S. Famiglietti, and B. Rajagopalan (2012),
        The influence of ENSO on global terrestrial water storage using GRACE, Geophysical
        Research Letters, 39 (16), L16,705, doi:10.1029/2012GL052495.
    
    Author:
        Amin Shakya, Interdisciplinary Center for Water Research (ICWaR), Indian Institute of Science (IISc)
    """
    c = fts.shape[1]
    
    ps = np.zeros((1, c))
    
    filter_ = ~np.isnan(fts)
    filter__ = ~np.isnan(ffts)
    
    fts_ = fts[filter_] #Extract values and leave Nan
    ffts_ = ffts[filter__] #Extract values and leave Nan
    
    fts = fts_.reshape(int(fts_.shape[0]/c),c)
    ffts = ffts_.reshape(int(ffts_.shape[0]/c),c)
    
    rn = fts.shape[0]
    
    for i in range(c):
        # A = np.concatenate(np.ones((rn,1)), np.real(signal.hilbert(ffts[:, i])), np.imag(signal.hilbert(ffts[:, i]))) #design matrix
        
        A = np.array((np.ones((rn)), np.real(signal.hilbert(ffts[:, i])), np.imag(signal.hilbert(ffts[:, i])))).T
        
        A = A.astype('double')
        B = fts[:,i]
        B = B.astype('double')
        abc = np.linalg.lstsq(A.T @ A, A.T @ B)[0]
        
        ps[0,i] = np.arctan2(abc[3-1],abc[2-1])*(180/np.pi) #check indices and degree/radian
    return ps



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
        
        ff = GSHS(shfil, quant, grd, n, h, jflag)[0]
        
        ff = ff*1000
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



def BasinAvg(data, path: str, c_rs, m, gs):
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
    
    Author: 
        Vivek Yadav, Interdisciplinary Center for Water Research (ICWaR), Indian Institute of Science (IISc)
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