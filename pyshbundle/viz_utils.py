# Visualisation Utilities for PySHBundle
# Author: Abhishek Mhamane, MS-Research Geoinformatics, IIT Kanpur (India)
# 2024-06-10, updated: Vivek Kumar Yadav, IISc Bengaluru
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

import calendar
from datetime import datetime
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.path as mpath
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import matplotlib.patches as mpatches

from pyshbundle.shutils import plm
from pyshbundle.pysh_core import gshs

def sc_triplot(scmat: np.ndarray, lmax: int, title: str, vmin, vmax):
    """
    Visualize the SH coefficients in SC triangular matrix format.

    Args:
        scmat (numpy.ndarray): SC matrix data (see clm2sc).
        lmax (int): Maximum degree of SH expansion.
        title (str): Title of the figure.
        vmin (float | int): Minimum value for the colorbar.
        vmax (float | int): Maximum value for the colorbar.

    Returns:
        (matplotlib.axes._axes.Axes): Plot axes.
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
    Visualize the SH coefficients in CS square matrix format.

    Args:
        csmat (numpy.ndarray): CS matrix data (see clm2cs or sc2cs).
        lmax (int): Maximum degree of SH expansion.
        title (str): Title of the figure.
        vmin (float): Minimum value for the colorbar.
        vmax (float): Maximum value for the colorbar.

    Returns:
        (matplotlib.axes._axes.Axes): Plot axes.
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
    """
    Visualize the polar regions of Greenland and Antarctica.

    Args:
        field (numpy.ndarray): The data field to visualize.
        polar_loc (str): The region to visualize, either 'greenland' or 'antarctica'.
        title (str): The title for the figure.
        file_name (str, optional): The file name along with the absolute path to save the figure. Defaults to None.
        save_flag (bool, optional): If True, the figure will be saved. Defaults to False.

    Returns:
        (matplotlib.axes._axes.Axes): Plot axes.
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
    """
    Visualize a field on a global map using the Robinson projection.

    Args:
        field (numpy.ndarray): The data field to visualize.
        img_extent (tuple): The extent of the image in the form (min_lon, max_lon, min_lat, max_lat).
        title (str): The title for the figure.
        name (str, optional): The file name along with the absolute path to save the figure. Defaults to None.
        colorbar_bounds (tuple, optional): The bounds for the colorbar in the form (vmin, vmax). Defaults to None.
        save_flag (bool, optional): If True, the figure will be saved. Defaults to False.

    Returns:
        fig (matplotlib.figure.Figure): Figure object.
        geo_ax (matplotlib.axes._axes.Axes): Plot axes.
    """
    # Plotting and Visualization
    
    fig = plt.figure(figsize=(16, 7.5))
    geo_ax = plt.axes(projection = ccrs.Robinson())

    

    # plot the data
    if colorbar_bounds is not None:
        vmin = colorbar_bounds[0]
        vmax = colorbar_bounds[1]
        im = geo_ax.imshow(field, origin='upper', extent=img_extent, cmap='Greens', transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax)
    else:
        im = geo_ax.imshow(field, origin='upper', extent=img_extent, cmap='Greens', transform=ccrs.PlateCarree(), )

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

    """
    Compute the spherical harmonics Ylm.

    Args:
        l (int): Degree, must be non-negative.
        m (int): Order, must be non-negative and less than or equal to l.
    Returns:
        tuple: A tuple containing two numpy arrays:
            - ylmc (numpy.ndarray): The real part of the spherical harmonics.
            - ylms (numpy.ndarray): The imaginary part of the spherical harmonics.
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
    
    
    p = plm(arr, m, thetaRAD, nargin=1, nargout=1)

    ylmc = p * cosml
    ylms = p * sinml

    return (ylmc, ylms)


def ylm_plot(l: int, m: int):
    """
    Plot the spherical harmonics Ylm.

    Args:
        l (int): Degree, must be non-negative.
        m (int): Order, must be non-negative and less than or equal to l.

    Returns:
        None
        Plot the spherical harmonics Ylm.
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

def gshs_prepare(lmax, gs, quant, grd, h, jflag, sc_coeff):
    """
    Prepare the grid for the given spherical harmonics coefficients.

    Args:
        lmax (int): Maximum degree of spherical harmonics.
        gs (float): Grid spacing in degrees.
        quant (int): Quantization level.
        grd (int): Grid type.
        h (float): Height parameter.
        jflag (int): Flag for the computation method.
        sc_coeff (numpy.ndarray): Spherical harmonics coefficients.

    Returns:
        (numpy.ndarray): The prepared grid field.
    """
    n = int(180/gs)
    
    grid_y = int(180/gs)
    grid_x = int(360/gs)

    ff = gshs(sc_coeff, quant, grd, n, h, jflag)[0]

    # rearranging
    field = np.zeros([grid_y,grid_x], dtype ='float')

    field[:,0:int(grid_x/2)] = ff[:,int(grid_x/2):]
    field[:,int(grid_x/2):] = ff[:,0:int(grid_x/2)]  
    
    return field

# Function to plot the calendar
def plot_calendar_months(datetime_object):
    """
    Plot a calendar for each year in the given list of datetime objects.

    Args:
        datetime_object (list): A list of datetime objects in the format '%Y-%m'.

    Returns:
        None

    This function takes a list of datetime objects and plots a calendar for each year in the list.
    The calendars are displayed in separate subplots, with each subplot representing a year.
    The function extracts the months and years from the datetime objects and determines the range of years to plot.
    For each year, the function creates a subplot and sets the title to the year.
    The function then highlights the months with replacement data by coloring the month names in the calendar.
    The color of the month names is 'lightblue' if the month is present in the given datetime objects, otherwise it is 'white'.
    """

    # Extract months and years from dictionary keys
    dates = [datetime.strptime(date, '%Y-%m').date() for date in datetime_object]
    months_years = {(date.year, date.month) for date in dates}
    
    # Determine the range of years to plot
    years = sorted({year for year, month in months_years})
    
    fig, axes = plt.subplots(nrows=len(years), ncols=1, figsize=(7, 1 * len(years)), dpi=300)

    if len(years) == 1:
        axes = [axes]
    
    for i, year in enumerate(years):
        ax = axes[i]
        ax.set_title(f'{year}')
        ax.set_xticks([])  # Hide y-axis
        ax.set_yticks([])  # Hide x-axis

        # Highlight months with replacement data
        for month in range(1, 13):
            if (year, month) in months_years:
                color = 'lightblue'
            else:
                color = 'white'
            ax.text(month - 1, 0, calendar.month_abbr[month], ha='center', va='center', 
                    bbox=dict(facecolor=color, edgecolor='black'))

        ax.set_xlim(-0.5, 11.5)
        ax.set_ylim(-0.5, 0.5)
        ax.grid(False)
        
    # # Add a legend
    white_patch = mpatches.Patch(color='white', label='Data unavailable')
    lightblue_patch = mpatches.Patch(color='lightblue', label='Data available')
    axes[0].legend(handles=[white_patch, lightblue_patch], loc='lower center', bbox_to_anchor=(0.85, 1.1), fontsize='9')
    plt.tight_layout()
    plt.show()