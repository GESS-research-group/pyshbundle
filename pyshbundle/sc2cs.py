#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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
    
# Authors: 
#    Dr. Bramha Dutt Vishwakarma, Interdisciplinary Center for Water Research (ICWaR), Indian Institute of Science (IISc)
#    Amin Shakya, Interdisciplinary Center for Water Research (ICWaR), Indian Institute of Science (IISc)

import numpy

def sc2cs(field):
    """converts the rectangular $(L+1) \times (2L+1)$ matrix FIELD, containing
    spherical harmonics coefficients in /S|C\ storage format into a 
    square (L+1)x(L+1) matrix in |C\S| format.

    Parameters:
        field (numpy.ndarray()):
            the rectangular (L+1)x(2L+1) matrix FIELD, containing
            spherical harmonics coefficients in /S|C\ storage format
        
    Returns: 
        cs (numpy.ndarray): 
            square (L+1)x(L+1) matrix in |C\S| format
    
    References:
        See the SHBundle docs or PySHBundle docs for more info about SH coeff. storage and retrival formats being implementd.

    Examples:
        >>> cs_fmt = sc2cs(field)
        TO DO: show suitable output
    """

    rows = len(field)
    cols = len(field[0])

    if (rows!=cols) and (cols!=2*rows - 1):
        sc2cs.exit("Input neither in cs nor in sc format")
    elif cols == rows:
        cs = field
    else:
        c    = field[:, rows-1:cols]
        st   = numpy.transpose(numpy.fliplr(field[:, 0:rows-1]))
        z    = numpy.zeros([1,rows])
        s    = numpy.concatenate((st, z), axis=0)
        cs   = numpy.add(c, s)
        
    return(cs)