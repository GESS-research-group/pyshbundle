# -*- coding: utf-8 -*-

import numpy

def cs2sc(field):

#     Created on Thu May  7 18:51:20 2020
#     CS2SC(FIELD) converts the square (L+1)x(L+1) matrix FIELD, containing
#     spherical harmonics coefficients in |C\S| storage format into a 
#     rectangular (L+1)x(2L+1) matrix in  /S|C\format.
    
#     IN:
#     field .... the square (L+1)x(L+1) matrix FIELD , containing
#                    spherical harmonics coefficients in |C\S| storage format
#     OUT: 
#     sc ....... rectangular (L+1)x(2L+1) matrix in  /S|C\format
    
# ----------------------------------------------------------------------------
#      project: GRACEpy
# ----------------------------------------------------------------------------
#     % author: Bramha Dutt Vishwakarma, University of Bristol
#     @author: bv18488

    rows = len(field)
    cols = len(field[0])

    if (rows!=cols) and (cols!=2*rows-1):
        raise Exception("Input neither in cs nor in sc format")
    elif cols==2*rows-1:
        sc = field
    else:
        c    = numpy.tril(field)
        ut   = numpy.triu(field)
        i = numpy.identity(rows)
        i = 1-i
        s    = numpy.fliplr(numpy.transpose(numpy.multiply(ut, i, )))
        sc   = numpy.concatenate((s[:,1:rows], c), axis=1)
        
    return(sc)