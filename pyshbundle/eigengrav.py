#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Created on Mon May 11 00:20:49 2020

# @author: Dr. Bramha Dutt Vishwakarma, Interdisciplinary Center for Water Research (ICWaR), Indian Institute of Science (IISc)

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

import numpy
from pyshbundle import GRACEpy as GB
from pyshbundle import GRACEconstants as GC


def eigengrav(lmax: int, fstr: str, h: float):
    """_summary_

    Args:
        lmax (int): Maximum degree of Spherical Coefficients
        fstr (str): gravity quantity, options: 'None', 'geoid', 'potential', 'gravity', 'tr', 'trr', 'slope'
                    'water', 'smd', 'height'
        h (float): _description_

    Returns:
        _type_: _description_

    Raises:
        TypeError: Enter a valid lmax value

    TO DO:
        Can we think about the raising a ValueError instead of instantly terminating the function
        Adding comments as variable names are not much descriptive
    """

    if type(lmax) == int:
        rows = 1
    else:
        rows = len(lmax)
    # rows = len(l)

    if rows > 1 or lmax < 0:
        raise TypeError("Enter a valid lmax value")

    r = GC.ae + h

    # lmax issue - using lmax as per used in shbundle
    # no reference for height

    if fstr == 'none':
        tf = numpy.ones((1, lmax+1))
    elif fstr == 'geoid':
        tf = numpy.ones((1, lmax+1)) * r
    elif fstr == 'potential':
        tf = numpy.ones((1, lmax+1)) * (GC.GM/r)
    elif fstr == 'gravity' or fstr == 'dg':
        tf = numpy.multiply(range(-1, lmax, 1), ((GC.GM/r/r) * 1e5))
        tf = tf.reshape((1, lmax+1))
    elif fstr == 'tr':
        tf = numpy.multiply(range(-1, -(lmax+2), -1), ((GC.GM/r/r) * 1e5))
        tf = tf.reshape((1, lmax+1))
    elif fstr == 'trr':
        tf = numpy.multiply(range(1, (lmax+2), 1),
                            range(2, (lmax + 3), 1))*((GC.GM/r/r) * 1e9)
        tf = tf.reshape((1, lmax+1))
    elif fstr == 'slope':
        tf = numpy.sqrt(numpy.multiply(
            range(0, lmax+1, 1), range(1, lmax+2, 1)))
        tf = tf.reshape((1, lmax+1))
    elif fstr == 'water':
        ln = GB.lovenr(lmax)
        tf = numpy.divide(numpy.multiply(
            5.517*r, numpy.add(range(0, 2*lmax + 1, 2), 1)), numpy.multiply(3, (1+ln)))
        tf = tf.reshape((1, lmax+1))
    elif fstr == 'smd':
        ln = GB.lovenr(lmax)
        tf = numpy.divide(numpy.multiply(
            5517*r, numpy.add(range(0, 2*lmax + 1, 2), 1)), numpy.multiply(3, (1+ln)))
        tf = tf.reshape((1, lmax+1))
    elif fstr == 'height':
        # not sure aobut heights - kept it unchanged ... abhishek
        kl, hl, ll = GB.lovenrPREM(90, 'CF')
        tf = numpy.divide(numpy.multiply(hl, (GC.ae*1000)), numpy.add(kl, 1))
    else:
        ValueError('Please choose a valid quantity for fstr')

    if h > 0:
        upConTerm = GB.upwcon(lmax, h)
        tf = numpy.multiply(tf, upConTerm)

    return(tf)
