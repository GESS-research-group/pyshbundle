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

def sc_triplot(scmat: np.ndarray, lmax: int, title: str, vmin, vmax):
    """_summary_

    Args:
        scmat (np.ndarray): _description_
        lmax (int): _description_
        title (str): _description_
        vmin (_type_): _description_
        vmax (_type_): _description_

    Returns:
        _type_: _description_
    """
    fig, ax = plt.subplots(1, 1, figsize=(25, 10))
    im = ax.imshow(np.ma.log10(abs(scmat)), extent=[-lmax, lmax, lmax, 0], cmap='Spectral_r',vmin=vmin, vmax=vmax)
    ax.grid()
    # plt.colorbar()
    x_vec = np.arange(-lmax, lmax+1, 6)
    y_vec = np.arange(lmax, -2, -6)

    x_st = 0*y_vec

    # vertical line
    ax.plot(x_st, y_vec, "black") 

    plt.xticks(x_vec,)
    plt.yticks(y_vec)
    plt.title(title)
    fig.colorbar(im,)    
    return ax

def cs_sqplot(csmat: np.ndarray, lmax: int, title: str, vmin, vmax):
    """

    Args:
        csmat (np.ndarray): _description_
        lmax (int): _description_
        title (str): _description_
        vmin (_type_): _description_
        vmax (_type_): _description_

    Returns:
        _type_: _description_
    """
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    im = ax.imshow(np.ma.log10(abs(csmat)), extent=[0, lmax, lmax, 0], cmap='Spectral_r',vmin=vmin, vmax=vmax)
    ax.grid()
    ax.set_aspect('equal')

    # plt.colorbar()
    x_vec = np.arange(0, lmax+1, 6)
    y_vec = np.arange(lmax, -2, -6)

    # diagonal line
    ax.plot(x_vec, np.flip(y_vec), "black") 

    # formating
    plt.xticks(x_vec,)
    plt.yticks(y_vec)
    plt.title(title)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=0.45)
    
    # Plot horizontal colorbar
    fig.colorbar(im, orientation="horizontal", cax=cax)
    
    return ax



def polar_plot(field, polar_loc: str, title, file_name=None, save_flag=False):
    """_summary_

    Args:
        field (_type_): _description_
        polar_loc (str): _description_
        title (_type_): _description_
        file_name (_type_, optional): _description_. Defaults to None.
        save_flag (bool, optional): _description_. Defaults to False.

    Returns:
        _type_: _description_
    """

    if polar_loc == 'greenland':
        extent = (-75, -5, 55, 85)

        # setting the 
        fig = plt.figure()
        ax = plt.axes(projection=ccrs.LambertConformal(central_latitude=72,central_longitude=-42.0))

        im = ax.imshow(field, origin='upper', cmap='Spectral', transform=ccrs.PlateCarree(), )

        ax.set_extent((-75, -5, 55, 85))
        gl = ax.gridlines(crs = ccrs.PlateCarree(), draw_labels=True, x_inline=False, y_inline=False, color='gray', alpha=0.9, linestyle='--')

        gl.top_labels = False
        ax.coastlines()
        plt.colorbar(im, orientation='vertical', shrink=1.0, pad=0.1,label=f"[...]")

        plt.title(f"{title}")
        if save_flag:
            plt.savefig(f"{file_name}.jpg")

    
    elif polar_loc == 'antarctica':
        extent = [-180, 180, -85, -60]

        fig = plt.figure(1, figsize=(7, 7))
        # setting the projection for polar plot
        ax = plt.axes(projection=ccrs.SouthPolarStereo())

        # plotting the matrix field
        im = ax.imshow(field, origin='upper', cmap='RdYlBu', transform=ccrs.PlateCarree(),)

        # setting the gridlines and continental boundary
        ax.set_extent(extent, ccrs.PlateCarree())
        gl = ax.gridlines(crs = ccrs.PlateCarree(), draw_labels=True, x_inline=False, y_inline=False, color='gray', alpha=0.9, linestyle='--')
        #gl.top_labels = False
        
        ax.coastlines()

        # to plot the circular plot boundary
        theta = np.linspace(0, 2*np.pi, 100)
        center, radius = [0.5, 0.5], 0.5
        verts = np.vstack([np.sin(theta), np.cos(theta)]).T
        circle = mpath.Path(verts * radius + center)

        ax.set_boundary(circle, transform=ax.transAxes)
        plt.colorbar(im, orientation='vertical', shrink=1.0, pad=0.1,label=f"[...]")
        plt.title(f"{title}")

        if save_flag:
            plt.savefig(f"{file_name}.jpg")
    return im

def mapfield(field, img_extent, title, name=None, colorbar_bounds=None, save_flag=False):
    # Plotting and Visualization
    
    fig = plt.figure(figsize=(16, 7.5))
    geo_ax = plt.axes(projection = ccrs.Robinson())

    

    # plot the data
    if colorbar_bounds is not None:
        vmin = colorbar_bounds[0]
        vmax = colorbar_bounds[1]
        im = geo_ax.imshow(field, origin='upper', extent=img_extent, cmap='RdYlBu', transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax)
    else:
        im = geo_ax.imshow(field, origin='upper', extent=img_extent, cmap='RdYlBu', transform=ccrs.PlateCarree(), )

    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    # setting gridlines
    gl = geo_ax.gridlines(crs = ccrs.PlateCarree(), draw_labels=True, x_inline=False, y_inline=False, color='gray', alpha=0.9, linestyle='--')
    # remove top x label
    gl.top_labels = False
    # change x label styles - font size ad colour
    gl.xlabel_style = {'size':12,}
    # left and right labels
    gl.left_labels = True
    gl.right_labels = False
    # coastlines
    geo_ax.coastlines()

    # Using new axes for colorbar
    

    plt.colorbar(im, shrink=0.845, orientation='vertical', pad=0.02,label=f"gravity [...]",)

    plt.title(f"{title}")
    if save_flag:
        plt.savefig(f"{name}.jpg")
    
    return fig, geo_ax

def ylm(l: int, m: int):
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


def ylm_plot(l: int, m: int):
    """_summary_

    Args:
        l (int): _description_
        m (int): _description_
    """
    ylmc, ylms = ylm(l, m)

    fig = plt.figure(figsize=(15, 7.5))
    ax = plt.axes(projection = ccrs.PlateCarree())

    lons = np.linspace(-180, 180, 73)
    lats = np.linspace(-90, 90, 37)

    x, y = np.meshgrid(lons, lats)

    if m >=0 :
        img_extent = (-180, 180, -90, 90)

        # plot the data
        im = ax.imshow(ylmc[:, 0, :], origin='upper', extent=img_extent, transform=ccrs.PlateCarree(), cmap="Spectral")
    else:
        #plt.contourf(x, y, ylms_00[:, 0, :], cmap='RdYlBu_r')
        im = ax.imshow(ylms[:, 0, :], origin='upper', extent=img_extent, transform=ccrs.PlateCarree(), cmap="Spectral")



    # setting gridlines
    gl = ax.gridlines(crs = ccrs.PlateCarree(), draw_labels=True, x_inline=False, y_inline=False, color='gray', alpha=0.9, linestyle='--')
    # remove top x label
    gl.top_labels = False
    # change x label styles - font size ad colour
    gl.xlabel_style = {'size':12,}
    # left and right labels
    gl.left_labels = True
    gl.right_labels = False
    # coastlines
    ax.coastlines()

    plt.colorbar(im, orientation='vertical', shrink=0.85, pad=0.02,label=f"[...]")

    plt.title(f"Visualization of Spherical Harmonics - degree: {l} order: {m}")

