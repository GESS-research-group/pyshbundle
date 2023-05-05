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

#    cs2sc(field) converts the square (L+1)x(L+1) matrix 'field', containing
#    spherical harmonics coefficients in |C\S| storage format into a 
#    rectangular (L+1)x(2L+1) matrix in  /S|C\format.
    
#    IN:
#    field .... the square (L+1)x(L+1) numpy matrix field , containing
#                   spherical harmonics coefficients in |C\S| storage format
#    OUT: 
#    sc ....... rectangular (L+1)x(2L+1) numpy matrix in  /S|C\format

#    @author: Amin Shakya, Interdisciplinary Center for Water Research (ICWaR), Indian Institute of Science (IISc)

import numpy

def cs2sc(field):
    """converts the square (L+1)x(L+1) matrix 'field', containing
    spherical harmonics coefficients in |C\S| storage format into a 
    rectangular (L+1)x(2L+1) matrix in  /S|C\ format.

    Args:
        field (np.ndarray): the square (L+1)x(L+1) numpy matrix field , containing
                   spherical harmonics coefficients in |C\S| storage format
    
    Returns:
        numpy.ndarray: Rectangular (L+1)x(2L+1) numpy matrix in  /S|C\ format

    Raises:
        TypeError: Input neither in cs nor in sc format
    
    Todo:
        + Rather use TypeError instead of base Exception
    
    Examples:
        >>> sc_shcoeff = cs2sc(cs_shcoeff)
        TO DO: write the output
    """
    rows = len(field)
    cols = len(field[0])

    if (rows != cols) and (cols != 2*rows - 1):
        raise TypeError("Input neither in cs nor in sc format")
    elif cols == 2*rows - 1:
        sc = field
    else:
        c    = numpy.tril(field)
        ut   = numpy.triu(field)
        i = numpy.identity(rows)
        i = 1-i
        s    = numpy.fliplr(numpy.transpose(numpy.multiply(ut, i, )))
        sc   = numpy.concatenate((s[:,1:rows], c), axis=1)
        
    return(sc)