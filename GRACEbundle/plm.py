#def funct(l, m = 0, thetaRAD = 0):
def plm(l,m,thetaRAD,nargin,nargout):

    
 # PLM Fully normalized associated Legendre functions for a selected order M
# Input as------
    '   l = np.array([np.arange(0,97,1)])            '
    '   m = 0                                        '
    '   thetaRAD = np.array([0,0.25,0.5,0.75,1])     '

 # HOW: 
 # '    p            = plm(l, m, thetaRAD,3,1)[:,:,0]		'	
 #    dp           = plm(l, m, thetaRAD,3,2)[1][:,:,0]
 #    ddp          = plm(l, m, thetaRAD,3,3)[2][:,:,0]
 #   

 # IN:
 #    l ........ degree (vector). Integer, but not necessarily monotonic.
 #               For l < m a vector of zeros will be returned.
 #    m ........ order (scalar). If absent, m = 0 is assumed.
 #    thetaRAD . co-latitude [rad] (vector)

 # OUT:
 #    p ........ Matrix with Legendre functions. The matrix has length(thetaRAD) 
 #               rows and length(l) columns, unless l or thetaRAD is scalar. 
 #               Then the output vector follows the shape of respectively l or 
 #               thetaRAD. 
 #    dp ....... Matrix with first derivative of Legendre functions. The matrix 
 #               has length(thetaRAD) rows and length(l) columns, unless l or 
 #               thetaRAD is scalar. Then the output vector follows the shape of
 #               respectively l or thetaRAD. 
 #    ddp ...... Matrix with second derivative of Legendre functions. The matrix 
 #               has length(thetaRAD) rows and length(l) columns, unless l or 
 #               thetaRAD is scalar. Then the output vector follows the shape of 
 #               respectively l or thetaRAD. 

 # SEE ALSO:
 #    LEGPOL, YLM, IPLM

 # REMARKS:
 #    Previous versions calculated the derivatives towards the latitude, 
 #    i. e. dP/d\phi are calculated. Please check your code in order to get 
 #    the derivatives correctly towards the co-latitude!
 #      ->  dP/d\thetaRAD      =   -dP/d\phi
 #      ->  d^2P/d\thetaRAD^2  =  d^2P/d\phi^2

 # -------------------------------------------------------------------------
 # project: SHBundle 
 # -------------------------------------------------------------------------
 # authors:
 #    Nico SNEEUW (NS), IAPG, TU-Munich
 #    Matthias WEIGELT (MW), DoGE, UofC
 #    Markus ANTONI (MA), GI, Uni Stuttgart 
 #    <bundle@gis.uni-stuttgart.de>
 # -------------------------------------------------------------------------
 # revision history:
 #    2013-01-29: MA, comments
 #    2013-01-23: MA, input argument thetaRAD [deg -> radian]
 #    2008-04-04: MW, extension for second derivative
 #    2004-11-24: MW, speed up calculation
 #    2004-08-13: MW, extension for first derivative
 #    1998-07-13: NS, Pmm non-recursive anymore
 #    1997-06-09: NS, help text brushed up
 #    1994-08-08: NS, initial version
 # -------------------------------------------------------------------------
 # license:
 #    This program is free software; you can redistribute it and/or modify
 #    it under the terms of the GNU General Public License as published by
 #    the  Free  Software  Foundation; either version 3 of the License, or
 #    (at your option) any later version.
  
 #    This  program is distributed in the hope that it will be useful, but 
 #    WITHOUT   ANY   WARRANTY;  without  even  the  implied  warranty  of 
 #    MERCHANTABILITY  or  FITNESS  FOR  A  PARTICULAR  PURPOSE.  See  the
 #    GNU General Public License for more details.
  
 #    You  should  have  received a copy of the GNU General Public License
 #    along with Octave; see the file COPYING.  
 #    If not, see <http://www.gnu.org/licenses/>.
    '@author: wslvivek'
 # -------------------------------------------------------------------------
 
    import sys
    import numpy as np
    if  min(l.shape) !=1:
        print('Degree l must be a vector (or scalar)') 
        sys.exit([])
    if  np.remainder(l,1).all() != 0:
        print('Vector l contains non-integers !!') 
        sys.exit([])
    # if 
    #     print('Order m must be a scalar') 
    #     sys.exit([])
    if  np.remainder(m,1) != 0:
        print('Order must be integer')
        sys.exit([])

     
# PRELIMINARIES
    lcol = len(l)
    trow = len(thetaRAD)
    lmax = int(max(l[0,:]))
    
    if lmax < m:
        p = np.zeros([len(thetaRAD), len(l)], dtype='longdouble')
        dp = np.zeros([len(thetaRAD), len(l)], dtype='longdouble')
        ddp = np.zeros([len(thetaRAD), len(l)], dtype='longdouble')
        sys.exit([])
    
    n = thetaRAD.size                                                               # number of latitudes
    t = thetaRAD[:]
    x = np.cos(t)
    y = np.sin(t)
    lvec = np.transpose(l)   
    lvec = np.intc(lvec)                                                       # l can now be used as running index
    
    if min(t).all() < 0 and max(t).all() > np.pi:
        print('Warning: Is co-latitude in radians ?')
    
    # Recursive computation of the temporary matrix ptmp, containing the Legendre
    # functions in its columns, with progressing degree l. The last column of
    # ptmp will contain zeros, which is useful for assignments when l < m.
    ptmp = np.zeros((n,lmax + 2 - m))
    if nargout >= 2:                                                                #  first derivative needs also P_{n,m+1} and P_{n,m-1}
        ptmp_m1 = np.zeros((n,lmax + 3 - m), dtype='longdouble')
        ptmp_p1 = np.zeros((n,lmax + 1 -m), dtype='longdouble')        
        dptmp = np.zeros((n,lmax + 2 - m), dtype='longdouble') 
    if nargout == 3:                                                                # second derivative needs also dP_{n,m+1} and dP_{n,m-1}
        dptmp_m1 = np.zeros((n,lmax + 3 -m), dtype='longdouble')
        dptmp_p1 = np.zeros((n,lmax + 1 -m), dtype='longdouble')
        ptmp_m2 = np.zeros((n,lmax + 4 -m), dtype='longdouble')                                         # but these first derivative need dP_{n,m+2} and dP_{n,m-2}
        ptmp_p2 = np.zeros((n,lmax - m), dtype='longdouble')
        ddptmp = np.zeros((n,lmax + 2 -m), dtype='longdouble')
     
    # sectorial recursion: PM (non-recursive, though)
    ptmp[:,0] = secrecur(m,y)
    if nargout >= 2:                                                                # frist derivative needs preceding and subsequent element
        if m > 0:    
            ptmp_m1[:,0] = secrecur(m-1,y)                                          # preceding elements
        if m < lmax: 
            ptmp_p1[:,0] = secrecur(m+1,y)                                          # subsequent elemtens
    if nargout == 3:                                                                # second derivative needs P_{n,m+2} and P_{n,m-2} as well
        if m > 1:           
            ptmp_m2[:,0] = secrecur(m-2,y)                                          # preceding elements
        if m < lmax-1: 
            ptmp_p2[:,0] = secrecur(m+2,y)                                          # subsequent elemtens
    
    # l-recursion: P
    ptmp = lrecur(ptmp,x,m,lmax);
    if nargout >= 2:                                                                # frist derivative needs preceding and subsequent element
        if m > 0:
            ptmp_m1 = lrecur(ptmp_m1,x,m-1,lmax)                                    # preceding elements
        if m < lmax:
            ptmp_p1 = lrecur(ptmp_p1,x,m+1,lmax)                                    # subsequent elemtens
    
    if nargout == 3:                                                                # second derivative needs P_{n,m+2} and P_{n,m-2} as well
        if m > 1:
            ptmp_m2 = lrecur(ptmp_m2,x,m-2,lmax)                                    # preceding elements
        if m < lmax-1:
            ptmp_p2 = lrecur(ptmp_p2,x,m+2,lmax)                                    # subsequent elemtens
            
    # now compute the derivatives 
    if nargout >= 2:                                                                # first derivative
        dptmp = derivALF(dptmp,ptmp_m1,ptmp_p1,m,lmax)
    if nargout == 3:                                                                # second derivative
        if m > 0:    
            dptmp_m1 = derivALF(dptmp_m1,ptmp,ptmp_m2,m-1,lmax)
        if m < lmax:
            dptmp_p1 = derivALF(dptmp_p1,ptmp,ptmp_p2,m+1,lmax)
        ddptmp = derivALF(ddptmp,dptmp_m1,dptmp_p1,m,lmax)
    
    
    # %--------------------------------------------------------------------
    # % The Legendre functions have been computed. What remains to be done, is to
    # % extract the proper columns from ptmp, corresponding to the vector lvec. 
    # % If l or thetaRAD is scalar the output matrix p reduces to a vector. It should
    # % have the shape of respectively thetaRAD or l in that case.
    # %--------------------------------------------------------------------
    lind       = (lvec < m)   #######################################################################################################################	 # index into l < m
    pcol       = lvec - m + 0			                                            # index into columns of ptmp
    pcol[lind] = np.ndarray((lmax-m+2-6)*np.ones((sum(sum(lind)),1)))	            # Now l < m points to last col.
    p      = ptmp[:,pcol]			                                                # proper column extraction 
    if nargout >= 2:
        dp =  dptmp[:,pcol]                                                         # proper column extraction 
    if nargout == 3: 
        ddp = ddptmp[:,pcol]                                                        # proper column extraction  
    if max(lvec.shape)==1  and min(thetaRAD.shape)==1 and (trow == 1):
        p = p.T
        if nargout >= 2:
            dp  = np.transpose(dp)
        if nargout == 3:
            ddp = np.transpose(ddp)
    if max(thetaRAD.shape)==1 and min(lvec.shape)==1  and (lcol == 1):
        p = p.T
        if nargout >= 2:
            dp  = dp.T  
        if nargout == 3:
            ddp = ddp.T

    if nargout == 1: return p
    if nargout == 2: return p,dp
    if nargout == 3: return p.dp.ddp
    
 # INLINE FUNCTIONS
 # function for the sectorial recursion, non-recursive though
def secrecur(m,y):
    import numpy as np
    if m == 0:
       fac = 1
    else:
       mm  = np.array([2*x for x in range(1,m+1)])
       fac = np.sqrt(2*np.prod((mm+1)/mm))
    out = fac*np.power(y,m)                                                         # The 1st column of ptmp
    return out


# % function for the l-recursion
def lrecur(inn,x,m,lmax):
    import numpy as np
    for ll in np.arange(int(m)+1,lmax+1,1):
       col   = ll - m +1			                                                # points to the next collumn of ptmp
       root1 = np.sqrt( (2*ll+1)*(2*ll-1)/((ll-m)*(ll+m)) ).real 
       root2 = np.sqrt( (2*ll+1)*(ll+m-1)*(ll-m-1) / ( (2*ll-3)*(ll-m)*(ll+m) ) ).real
    
       # % recursion 
       if ll == m+1:
           inn[:,col-1] = root1 *x*inn[:,col-2];
       else:
           inn[:,col-1] = root1 *x*inn[:,col-2] - root2 *inn[:,col-3] 
    return inn



# function to calculate the derivate
def derivALF(inn,miin,plin,m,lmax):
    import numpy as np
    l = np.arange(m,lmax+2,1)
    #l=l.reshape(l.shape[0],1)
    if m == 0:
        inn[:,0] = 0
        if lmax > m:             
            inn[:,1:] = plin*(-np.sqrt(    ((l[1:]+1)*l[1:]   /2).real))            # (-ones(n,1)*realsqrt((l(2:end)+1).*l(2:end)./2)).*plin            
    elif m == 1:
        inn[:,0] = miin[:,1]
        if lmax > m: 
            inn[:,1:] =  miin[:,2:]*(np.sqrt((l[1:]+1)*l[1:]/2).real) -0.5*plin*(np.sqrt((l[1:]-1)*(l[1:]+2)).real)
    elif m == lmax:
        inn[:,0] = np.sqrt(m/2*miin[:,1:]).real
    else:
        inn[:,0] = np.sqrt((m/2)*miin[:,1:]).real
        if lmax > m: 
            inn[:,1:] = 0.5*miin[:,2:]*np.sqrt((l[:,1:]+m)*(l[:,1:]-m+1)).real - 0.5*plin*(np.sqrt((l[:,1:]-m)*(l[:,1:]+m+1)).real)
    return inn
    