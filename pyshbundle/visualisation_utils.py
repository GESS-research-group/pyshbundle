# Visualisation Utilities for PySHBundle
# Author: Abhishek Mhamane, MS-Research Geoinformatics, IIT Kanpur (India)
# 

from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
import matplotlib.path as mpath
import matplotlib.pyplot as plt
import numpy as np
import cartopy
import cartopy.crs as ccrs

from pyshbundle import sc2cs, clm2sc
from pyshbundle import plm
import pyshbundle

def sctriplot(scmat: np.ndarray | list, lmax: int, units=''):
    """SCTRIPLOT plots the triangles of SH coefficients stored in SC-format. 
    If stored in CS format it is converted to SC-format before plotting.

    Args:
        fig (_type_): Matplotlib figure object
        scmat (np.ndarray | list): _descrMatrix of real SH coefficients in SC, CS or [l m Clm Slm] formats.iption_
        lmax (int): Maximum degree of spherical harmonic expansion
    
    Returns:
        _type_: Generates an image of the SC-formatted SH coefficients
    
    To Do:
        1. pass figure object from user
        2. or provide no plot changing options
    """

    fmt, l_max = checkformat(scmat)

    if fmt == 'clm':
        scmat = clm2sc(scmat)
    
    r, c = np.shape(scmat)

    if fmt == 'cs':
        scmat = sc2cs(scmat)
    elif (r > c) and (r == (2*lmax) + 1) and (c == lmax+1):
        scmat = scmat.T
    elif (min(r,c) != lmax+1) or (max(r,c) != (2*lmax + 1)):
        raise ValueError('Matrix neither confirms to SC-format, nor CS-format')
    
    fig, ax = plt.subplots(figsize=())

    # using masked array to avoid log(0) error
    im = ax.imshow(np.ma.log10(np.abs(scmat)), cmap='Spectral_r', extent=[-lmax, lmax+1, lmax+1, 0]) 

    ax.set_aspect('auto')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.xlabel("Order [m]")
    plt.ylabel("Degre [l]")

    fig.colorbar(im, ax=ax, label=f"[{units}]", orientation='vertical', pad=0.02)

    return fig, ax

def polar_plot(field, flag: str, dpi=100):

    if flag == 'greenland':
        extent = [-75, -5, 55, 85]

        # setting the 
        fig = plt.figure(1, figsize=(5, 5), dpi=dpi)
        ax = plt.axes(projection=ccrs.LambertConformal(central_latitude=72,central_longitude=-42.0))

        im = ax.imshow(field, origin='upper',  extent=extent, cmap='RdYlBu_r', transform=ccrs.PlateCarree(),)

        ax.set_extent(extent)
        gl = ax.gridlines(crs = ccrs.PlateCarree(), draw_labels=True, x_inline=False, y_inline=False, color='gray', alpha=0.9, linestyle='--')

        gl.top_labels = False
        ax.coastlines(resolution='50m')
        plt.colorbar(im)
    
    elif flag == 'antarctica':
        extent = [-180, 180, -85, -60]

        fig = plt.figure(1, figsize=(5, 5), dpi=dpi)
        # setting the projection for polar plot
        ax = plt.axes(projection=ccrs.SouthPolarStereo())

        # plotting the matrix field
        im = ax.imshow(field, origin='upper', extent=extent, cmap='RdYlBu_r', transform=ccrs.PlateCarree(),)

        # setting the gridlines and continental boundary
        ax.set_extent(extent, ccrs.PlateCarree())
        gl = ax.gridlines(crs = ccrs.PlateCarree(), draw_labels=True, x_inline=False, y_inline=False, color='gray', alpha=0.9, linestyle='--')
        #gl.top_labels = False
        
        ax.coastlines(resolution='50m')

        # to plot the circular plot boundary
        theta = np.linspace(0, 2*np.pi, 100)
        center, radius = [0.5, 0.5], 0.5
        verts = np.vstack([np.sin(theta), np.cos(theta)]).T
        circle = mpath.Path(verts * radius + center)

        ax.set_boundary(circle, transform=ax.transAxes)
    pass


def surface_spherical_hormonics(l: int, m: int):
    """_summary_

    Args:
        l (int): Degree
        m (int): Order

    Returns:
        _type_: _description_
    """

    # input handling
    assert l >= 0

    assert m >= 0

    assert m <= l

    # main code
    thetaRAD  = np.linspace(0,np.pi,37)
    lambdaRAD = np.linspace(0,2*np.pi,73)

    cosml = np.cos(m*lambdaRAD)
    sinml = np.sin(m*lambdaRAD)

    arr = np.zeros((1,1))
    arr[0] = l
    
    
    p = plm.plm(arr, m, thetaRAD, nargin=1, nargout=1)

    ylmc = p * cosml
    ylms = p * sinml

    return (ylmc, ylms)