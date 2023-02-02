# -*- coding: utf-8 -*-
"""
Created on Mon May 11 00:20:49 2020

@author: bv18488
"""
# ----------------------------------------------------------------------------
# license:
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the  Free  Software  Foundation; either version 3 of the License, or
#    (at your option) any later version.
#  
#    This  program is distributed in the hope that it will be useful, but 
#    WITHOUT   ANY   WARRANTY;  without  even  the  implied  warranty  of 
#    MERCHANTABILITY  or  FITNESS  FOR  A  PARTICULAR  PURPOSE.  See  the
#    GNU General Public License for more details.
#  
#    You  should  have  received a copy of the GNU General Public License
#    along with Octave; see the file COPYING.  
#    If not, see <http://www.gnu.org/licenses/>.
# --------------------------------------------------------------------------
def eigengrav(lmax,fstr, h):
    import numpy
    import GRACEpy as GB
    import GRACEconstants as GC
    
    if type(lmax) == int:
        rows = 1
    else:
        rows = len(lmax)
#    rows = len(l)
    
    if rows>1 or lmax<0:
        raise Exception("please input a valid value for lmax")
        
        
    r = GC.ae + h
    
    
    if fstr == 'none':
        tf = numpy.ones((lmax+1, 1))
    elif fstr == 'geoid':
        tf = numpy.ones((lmax+1, 1)) * r
    elif fstr == 'potential':
        tf = numpy.ones((lmax+1, 1)) * (GC.GM/r)
    elif fstr == 'gravity' or fstr == 'dg':
        tf =  numpy.multiply( range(-1,lmax,1), ((GC.GM/r/r) * 1e5))
    elif fstr == 'tr':
        tf =  numpy.multiply( range(-1,-(lmax+2),-1), ((GC.GM/r/r) * 1e5))
    elif fstr == 'trr':
        tf =  numpy.multiply( range(1,(lmax+2),1), range(2,(lmax + 3),1))*((GC.GM/r/r) * 1e9)
    elif fstr == 'slope':
        tf = numpy.sqrt(numpy.multiply(range(0,lmax+1,1), range(1,lmax+2,1)) )
    elif fstr == 'water':
        ln = GB.lovenr(lmax)
        tf = numpy.divide(numpy.multiply(5.517*r, numpy.add(range(0,2*lmax + 1,2),1)) , numpy.multiply(3,(1+ln)))
    elif fstr == 'smd':
        ln = GB.lovenr(lmax)
        tf = numpy.divide(numpy.multiply(5517*r, numpy.add(range(0,2*lmax + 1,2),1)) , numpy.multiply(3,(1+ln)))
    elif fstr == 'height':
        kl,hl,ll = GB.lovenrPREM(90,'CF')
        tf  = numpy.divide(numpy.multiply(hl,(GC.ae*1000)),numpy.add(kl,1))
    else:
         eigengrav.exit('Please choose a valid quantity for fstr')  
         
    if h>0:
        upConTerm = GB.upwcon(lmax,h)
        tf = numpy.multiply(tf, upConTerm)
    
    return(tf)
   