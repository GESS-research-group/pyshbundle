#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Created on Sat May  9 18:49:45 2022
# This file contains some basic constants:
#  - physical
#  - geodetic (GRS80)
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

""" This script contains some of the major relavant Physical and Geodetic(GRS80) constants:

    + `clight` speed of light - $2.99792458e+8$ $m/s$
    + `G` Gravitational constant- $6.67259e-11$ $\frac{m^3} {kg \cdot s^2}$
    + `au` astronomical unit - $149.597870691e+9$ $m$

    + `ae` semi-major axis of ellipsoid `GRS 80`- $6378137$ m
    + `GM` geocentric grav. constant `GRS 80`- $3.986005e+14$ $\frac{m^3}{s^2}$
    + `J2` earth's dynamic form factor `GRS 80` - $1.08263e-3$ [unitless C20 unnormalized coefficient]
    + `Omega` mean ang. velocity `GRS 80` - $7.292115e-5 $\frac{rad}{s}$

    + `flat` flattening - $\frac{1}{298.257222101}$
"""

clight = 2.99792458e8	      # speed of light [m/s]
G = 6.67259e-11         # gravitational constant [m^3 /(kg s^2)]
au = 149.597870691e9     # astronomical unit [m]

# GRS80 defining constants:
ae = 6378137             # semi-major axis of ellipsoid [m]
GM = 3.986005e14         # geocentric grav. constant [m^3 / s^2]
J2 = 1.08263e-3          # earth's dyn. form factor (= -C20 unnormalized)
Omega = 7.292115e-5         # mean ang. velocity [rad/s]

# GRS80 derived constants:
flat = 1/298.257222101     # flattening
J4 = -0.237091222e-5     # -C40 unnormalized
J6 = 0.608347e-8        # -C60 unnormalized
J8 = -0.1427e-10         # -C80 unnormalized
