# -*- coding: utf-8 -*-

import numpy

def sc2cs(field):
    """
    Created on Thu May  7 18:51:20 2020
    % SC2CS(FIELD) converts the rectangular (L+1)x(2L+1) matrix FIELD, containing
    % spherical harmonics coefficients in /S|C\ storage format into a 
    % square (L+1)x(L+1) matrix in |C\S| format.
    %
    % IN:
    %    field .... the rectangular (L+1)x(2L+1) matrix FIELD, containing
    %               spherical harmonics coefficients in /S|C\ storage format
    %
    % OUT: 
    %    cs ....... square (L+1)x(L+1) matrix in |C\S| format
    
    % ----------------------------------------------------------------------------
    % project: GRACEpy 
    % ----------------------------------------------------------------------------
    % author: Bramha Dutt Vishwakarma, University of Bristol
    @author: bv18488
    """

    rows = len(field)
    cols = len(field[0])

    if (rows!=cols) and (cols!=2*rows-1):
        sc2cs.exit("Input neither in cs nor in sc format")
    elif cols==rows:
        cs = field
    else:
        c    = field[:,rows-1:cols]
        st   = numpy.transpose(numpy.fliplr(field[:,0:rows-1]))
        z    = numpy.zeros([1,rows])
        s    = numpy.concatenate((st, z), axis = 0)
        cs   = numpy.add(c, s)
        
    return(cs)