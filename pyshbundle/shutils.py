# Spherical Harmonics Utilities
# Curator: Abhishek Mhamane
# Updated: Vivek, 2024-10-04
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

import sys
import numpy as np
from scipy.fft import ifft
from scipy import sparse
from scipy.interpolate import PchipInterpolator

from pyshbundle import GRACEpy as GB
from pyshbundle import GRACEconstants as GC


def plm(l: np.array, m:int, thetaRAD, nargin, nargout): 
    """
    Fully normalized associated Legendre functions for a selected order M.

    Args:
        l (numpy.array): Degree, but not necessarily monotonic. For l < m a vector of zeros will be returned.
        m (int): Order. If absent, m = 0 is assumed.
        thetaRAD (numpy.array): Co-latitude in radians.
        nargin (int): Number of input arguments.
        nargout (int): Number of output arguments.

    Returns:
        (numpy.array): Fully normalized Legendre functions.
        (numpy.array): First derivative of the Legendre functions.
        (numpy.array): Second derivative of the Legendre functions.

    Author:
        Vivek Kumar Yadav, Interdisciplinary Center for Water Research (ICWaR), Indian Institute of Science (IISc).
    """

    if  min(l.shape) != 1:
        raise ValueError('Degree l must be a vector (or scalar)') 
    if  np.remainder(l,1).all() != 0:
        raise ValueError('Vector l contains non-integers !!') 
    if  np.remainder(m,1) != 0:
        raise ValueError('Order must be integer')

# PRELIMINARIES
    lcol = len(l)
    trow = len(thetaRAD)
    lmax = int(max(l[0,:]))
    
    if lmax < m:
        p = np.zeros([len(thetaRAD), len(l)], dtype='float')
        dp = np.zeros([len(thetaRAD), len(l)], dtype='float')
        ddp = np.zeros([len(thetaRAD), len(l)], dtype='float')
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
        ptmp_m1 = np.zeros((n,lmax + 3 - m), dtype='float')
        ptmp_p1 = np.zeros((n,lmax + 1 -m), dtype='float')        
        dptmp = np.zeros((n,lmax + 2 - m), dtype='float') 
    if nargout == 3:                                                                # second derivative needs also dP_{n,m+1} and dP_{n,m-1}
        dptmp_m1 = np.zeros((n,lmax + 3 -m), dtype='float')
        dptmp_p1 = np.zeros((n,lmax + 1 -m), dtype='float')
        ptmp_m2 = np.zeros((n,lmax + 4 -m), dtype='float')                                         # but these first derivative need dP_{n,m+2} and dP_{n,m-2}
        ptmp_p2 = np.zeros((n,lmax - m), dtype='float')
        ddptmp = np.zeros((n,lmax + 2 -m), dtype='float')
     
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
    
    
# --------------------------------------------------------------------
        # The Legendre functions have been computed. What remains to be done, is to
        # extract the proper columns from ptmp, corresponding to the vector lvec. 
        # If l or thetaRAD is scalar the output matrix p reduces to a vector. It should
        # have the shape of respectively thetaRAD or l in that case.
# --------------------------------------------------------------------
    lind       = (lvec < m)   	 # index into l < m
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

    if nargout == 1: 
        return p
    if nargout == 2: 
        return p,dp
    if nargout == 3: 
        return p, dp, ddp
    

  
def secrecur(m, y):
    """
    Helper Function for sectorial recursion.
    This function computes the sectorial recursion for given parameters.

    Args:
        m (int): The order of the recursion.
        y (numpy.ndarray): The input array for which the recursion is computed.

    Returns:
        numpy.ndarray: The result of the sectorial recursion.
    """
    if m == 0:
       fac = 1
    else:
       mm  = np.array([2*x for x in range(1, m+1)])
       fac = np.sqrt(2*np.prod((mm+1)/mm))
    out = fac*np.power(y,m)                                                         # The 1st column of ptmp
    return out


# % function for the l-recursion

def lrecur(inn, x, m, lmax):
    """
    Helper function for recursion.

    Args:
        inn (int): Input value representing the initial condition.
        x (int): The current value for the recursion.
        m (int): Order of the recursion.
        lmax (int): Maximum value for recursion.

    Returns:
        (int): Updated value after performing the recursion based on parameters.
    """
    for ll in np.arange(int(m)+1,lmax+1,1):
       col   = ll - m+1			                                                # points to the next collumn of ptmp
       root1 = np.sqrt( (2*ll+1)*(2*ll-1)/((ll-m)*(ll+m)) ).real 
       root2 = np.sqrt( (2*ll+1)*(ll+m-1)*(ll-m-1) / ( (2*ll-3)*(ll-m)*(ll+m) ) ).real
    
       # % recursion 
       if ll == m+1:
           inn[:, col-1] = root1 *x*inn[:, col-2]
       else:
           inn[:, col-1] = root1 *x*inn[:, col-2] - root2 *inn[:, col-3] 
    return inn



# function to calculate the derivate

def derivALF(inn, miin, plin, m, lmax):
    """
    Function to calculate the derivative of the associated Legendre functions.

    Args:
        inn (np.ndarray): Input array representing the initial condition.
        miin (np.ndarray): Array for the preceding elements in the recursion.
        plin (np.ndarray): Array for the subsequent elements in the recursion.
        m (int): Order of the associated Legendre functions.
        lmax (int): Maximum degree.

    Returns:
        (numpy.ndarray): Derivatives of the associated Legendre functions.
    """
    l = np.arange(m,lmax+2,1)
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


def iplm(l, m:int, theRAD, dt=-9999):
    """
    Integrals of the fully normalized associated Legendre functions over blocks for a selected order M.

    Args:
        l (numpy.array): Degree (vector). Integer, but not necessarily monotonic.
            For l < m a vector of zeros will be returned.
        m (int): Order of the Legendre function. If absent, m = 0 is assumed.
        theRAD (numpy.array): Co-latitude in radians.
        dt (int, optional): Integration block-size [rad] (scalar). Defaults to -9999.

    Returns:
        (numpy.ndarray): Matrix with integrated Legendre functions.
            Functions are integrated from theRAD(i)-dt/2 till theRAD(i)+dt/2.
            The matrix has length(TH) rows and length(L) columns, unless L 
            or TH is scalar. Then the output vector follows the shape of 
            respectively L or TH.

    Notes:
        The blocks at the pole might become too large under circumstances.
        This is not treated separately, i.e. unwanted output may appear.
        In case TH is scalar, dt will be 1 (arbitrarily).

    Uses:
        `plm`

    Author:
        Vivek Kumar Yadav, Interdisciplinary Center for Water Research (ICWaR), Indian Institute of Science (IISc)
    """

    if dt == -9999:
        dt = np.pi / 180 if len(theRAD) == 1 else theRAD[1] - theRAD[0]

    if min(l.shape) != 1:
        raise ValueError('Degree l must be a vector (or scalar)')

    if not np.all(np.remainder(l, 1) == 0):
        raise ValueError('Vector l contains non-integers !!')

    if not np.remainder(m, 1) == 0:
        raise ValueError('Order must be integer')

    if dt == 0:
        raise ValueError('DT cannot be zero')


    lcol = len(l)
    trow = len(theRAD)
    n = len(theRAD)
    theRAD.T
    if min(theRAD) < 0 or max(theRAD) > np.pi:
        raise ValueError('Is the co-latitude ''theta'' given in radian?')

    lmax = max(l[0])
    mfix = m
    lvec = np.transpose(l) 
    l = np.arange(mfix,lmax+1,1)
    
    # Initialization of cosine, sine and Plm functions
    stplus  = np.sin(theRAD+dt/2)
    stmin   = np.sin(theRAD-dt/2)
    ctplus  = np.cos(theRAD+dt/2)
    ctmin   = np.cos(theRAD-dt/2)
    plmplus = np.ones([n,lmax+1])
    plmmin = np.ones([n,lmax + 1])
    plmplus[:,l] = plm(np.array([l]),mfix,(theRAD + dt/2),3,1)[:,:,0]                  # Tesserals
    plmmin[:,l] = plm(np.array([l]),mfix,(theRAD - dt/2),3,1)[:,:,0] 
    if mfix > 0:
        m = np.arange(1,mfix + 1,1)
        mm = 2*m
        fac = np.sqrt(2*np.cumprod((mm+1)/mm))
        mgr, stp = np.meshgrid(m, stplus)
        fgr, stm = np.meshgrid(fac, stmin)
        plmplus[:, m] = fgr * np.power(stp, mgr)
        plmmin[:, m] = fgr * np.power(stm, mgr)
    ptmp = np.zeros([n, lmax +2 ])
    ptmp00 = np.cos(theRAD - dt/2) - ctplus
    ptmp11 = np.sqrt(3)/2 * (dt - ctplus* stplus + ctmin* stmin)
    ptmp10 = np.sqrt(3)/2 * (np.power(stplus,2) - np.power(stmin,2))
    ptmp[:,0] = ptmp00
    
    # Compute first the integrals of order m == 0
    if mfix == 0:
        ptmp[:,1] = ptmp10
        for l in range(2,lmax+1,1):              #loop over the degree l 
            rootnm = np.sqrt( (2*l+1)*(2*l-1)/np.power(l,2))
            root1nm = np.sqrt( (2*l-1)*(2*l-3)/np.power(l-1,2))
            ptmp[:,l] = rootnm/(l+1)*(((l-2)*ptmp[:,l-2]/root1nm).T + np.power(stplus,2)*plmplus[:,l-1].T - np.power(stmin,2)*plmmin[:,l-1].T )
    else:
        # Compute the integrals of order m > 0

        # First we compute the diagonal element IPmm (lmax == mfix)

        ptmp[:,1] = ptmp11
        for l in range(2,mfix+1,1):
            # print(l)
            rootmm = np.sqrt( (2*l+1)/(2*l) )
            root1mm = np.sqrt( (2*l-1)/(2*l-2))
            if l == 2:
                root1mm = np.sqrt(3)
            
            ptmp[:,l] = rootmm/(l+1)*( l*root1mm*ptmp[:,l-2].T - (ctplus*plmplus[:,l].T -ctmin*plmmin[:,l].T)/rootmm )
    #the arbitrary element IPlm ( computed only when lmax > mfix)        
        if lmax > mfix:
            l = mfix + 1
        #first we do the element IPlm, for which l - m = 1    
            rootnm = np.sqrt( (2*l+1)*(2*l-1)/(l+mfix)/(l-mfix))
            ptmp[:,l] = rootnm/(l+1)*(np.power(stplus,2)*plmplus[:,l-1].T - np.power(stmin,2)*plmmin[:,l-1].T)
        #now we do the rest
            for l in range(mfix+2,lmax+1,1):              #loop over the degree l
                rootnm = np.sqrt( (2*l+1) * (2*l-1) / (l+mfix) / (l-mfix) )
                root1nm = np.sqrt( (2*l-1) * (2*l-3) / (l-1+mfix) / (l-1-mfix) )
                ptmp[:,l] = rootnm/(l+1)*( (l-2)*ptmp[:,l-2].T/root1nm + np.power(stplus,2)*plmplus[:,l-1].T -np.power(stmin,2)*plmmin[:,l-1].T)

# --------------------------------------------------------------------
        # The integrated functions have been computed. What remains to be done, is to
        # extract the proper columns from ptmp, corresponding to the vector lvec. 
        # If l or theta is scalar the output matrix p reduces to a vector. It should
        # have the shape of respectively theta or l in that case.
# --------------------------------------------------------------------

# p     = zeros(n, length(lvec))
    lind = np.argwhere(lvec<mfix)[:,0]      #index into l < m
    pcol = lvec + 1                         #index into columns of ptmp
    pcol[lind] = (lmax + 2)*np.ones([len(lind),1])   #Now l < m points to last col
    p = ptmp[:,pcol[:,0]-1]                 #proper column extraction 
    
    if max(lvec.shape) == 1 and min(np.array([theRAD]).shape) == 1 and trow == 1:
        p = p.T
    if max(np.array([theRAD]).shape) == 1 and min(lvec.shape) == 1 and lcol == 1:
        p = p.T
    return p


def ispec(a,b = -9999):
    """
    Returns the function F from the spectra A and B.

    Args:
        a (numpy.ndarray): Cosine coefficients.
        b (numpy.ndarray, optional): Sine coefficients. Defaults to -9999.

    Returns:
        (numpy.ndarray): The function F computed from the spectra A and B.

    See Also:
        `spec`

    Author:
        Amin Shakya, Interdisciplinary Center for Water Research (ICWaR), Indian Institute of Science (IISc).
    """
    
    n2 = a.shape[0]
    a[0,:] = a[0, :]*2

    
    if (np.absolute(b[n2-1,:]) < 1e-10).all():
        n = 2 * n2 - 2     
        a[n2-1,:] = a[n2-1,:] * 2            
        fs = (a - 1j * b)/2
        fs  = (np.concatenate((fs,np.conj(fs[np.arange(n2-2,0,-1),:])), axis = 0))*max(n,1)

    else:
        n = 2 * n2 - 1                        
        fs = (a - 1j * b)/2
        fs = (np.concatenate((fs,np.conj(fs[np.arange(n2-1,0,-1),:])), axis = 0))*n

    f = np.real(ifft(fs.T).T)
    return f


def eigengrav(lmax: int, fstr: str, h: float):
    """
    Returns the isotropic spectral transfer (or: eigenvalues) of several gravity related quantities. 
    Upward continuation may be included.

    Args:
        lmax (int): Maximum degree of Spherical Coefficients.
        fstr (str): Denoting the functional under consideration:
            'none', 
            'geoid',
            'dg', 'gravity' ... gravity anomaly,
            'potential', 
            'tr' .............. gravity disturbance, 
            'trr' ............. (d^2/dr^2)
            'slope' ........... size of surface gradient, 
            'water' ........... equivalent water thickness, 
            'smd' ............. surface mass density.
            'height' .......... vertical displacements.
        h (float): Height above Earth mean radius [m].

    Returns:
        (np.ndarray): Transfer matrix. Size and shape equal to lmax. 
                    Units are respectively [none], [m], [mGal], [mGal], [E], [m^2/s^2], [rad], [m], [kg/m^2] [n x 1]

    Uses:
        upwcon, lovenr, uberall/constants, uberall/isint

    Raises:
        TypeError: Enter a valid lmax value.
    
    Author:
        Dr. Bramha Dutt Vishwakarma, Interdisciplinary Center for Water Research (ICWaR), Indian Institute of Science (IISc).
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
        tf = np.ones((1, lmax+1))
    elif fstr == 'geoid':
        tf = np.ones((1, lmax+1)) * r
    elif fstr == 'potential':
        tf = np.ones((1, lmax+1)) * (GC.GM/r)
    elif fstr == 'gravity' or fstr == 'dg':
        tf = np.multiply(range(-1, lmax, 1), ((GC.GM/r/r) * 1e5))
        tf = tf.reshape((1, lmax+1))
    elif fstr == 'tr':
        tf = np.multiply(range(-1, -(lmax+2), -1), ((GC.GM/r/r) * 1e5))
        tf = tf.reshape((1, lmax+1))
    elif fstr == 'trr':
        tf = np.multiply(range(1, (lmax+2), 1),
                            range(2, (lmax + 3), 1))*((GC.GM/r/r) * 1e9)
        tf = tf.reshape((1, lmax+1))
    elif fstr == 'slope':
        tf = np.sqrt(np.multiply(
            range(0, lmax+1, 1), range(1, lmax+2, 1)))
        tf = tf.reshape((1, lmax+1))
    elif fstr == 'water':
        ln = GB.lovenr(lmax)
        tf = np.divide(np.multiply(
            5.517*r, np.add(range(0, 2*lmax + 1, 2), 1)), np.multiply(3, (1+ln)))
        tf = tf.reshape((1, lmax+1))
    elif fstr == 'smd':
        ln = GB.lovenr(lmax)
        tf = np.divide(np.multiply(
            5517*r, np.add(range(0, 2*lmax + 1, 2), 1)), np.multiply(3, (1+ln)))
        tf = tf.reshape((1, lmax+1))
    elif fstr == 'height':
        kl, hl, ll = GB.lovenrPREM(90, 'CF')
        tf = np.divide(np.multiply(hl, (GC.ae*1000)), np.add(kl, 1))
    else:
        ValueError('Requested functional FSTR not available.')

    if h > 0:
        upConTerm = GB.upwcon(lmax, h)
        tf = np.multiply(tf, upConTerm)

    return(tf)


def grule(n: int):
    """
    Computes Gauss base points and weight factors using the algorithm.

    Args:
        n (int): Number of base points required.

    Returns:
        tuple: A tuple containing:
            - bp (numpy.ndarray): Cosine of the base points.
            - wf (numpy.ndarray): Weight factors for computing integrals and such.

    References:
        - 'Methods of Numerical Integration' by Davis and Rabinowitz, page 365, Academic Press, 1975.

    Author:
        Amin Shakya, Interdisciplinary Center for Water Research (ICWaR), Indian Institute of Science (IISc)
    """
    bp = np.zeros((n,1))
    wf = bp
    iter = 2
    m = np.floor((n+1)/2)
    e1 = n * (n+1)
    
    
    mm = 4*m - 1
    t = (np.pi / (4*n + 2)) * np.arange(3,mm+4,4)
    nn = (1 - (1 - 1/n)/(8*n*n))
    x0 = nn * np.cos(t)
    
    
    for i in range(iter):
        pkm1 = 1
        pk = x0
        
        for kk in range(n-1):
            k = kk + 2
            t1 = x0 * pk
            pkp1 = t1 - pkm1 - (t1-pkm1)/k  + t1
            pkm1=pk
            pk=pkp1
            
        den = 1 - x0*x0
        d1 = n * (pkm1 - x0*pk)
        dpn = d1/den
        
        
        d2pn = (2*x0*dpn - e1*pk) / den
        d3pn = (4*x0*d2pn + (2-e1)*dpn)/den
        d4pn = (6*x0*d3pn + (6-e1)*d2pn)/den
        u = pk/dpn
        v = d2pn/dpn
        h = -u * (1+(.5*u)*(v+u*(v*v - u*d3pn/(3*dpn))))
        p = pk + h*(dpn+(0.5*h)*(d2pn+(h/3)*(d3pn + 0.25*h*d4pn)))
        dp = dpn + h*(d2pn+(0.5*h)*(d3pn+h*d4pn/3))
        h = h-p/dp
        x0 = x0+h
    
    bp = -x0-h
    fx = d1 - h*e1*(pk+(h/2)*(dpn+(h/3)*(d2pn+(h/4)*(d3pn+(0.2*h)*d4pn))))
    wf = [2 * (1 - np.power(bp,2))]/(fx*fx)
                 
                
    for i in range(len(bp),n):
        bp = np.append(bp,[0])
        wf = np.append(wf,[0])
    
    if ((m)+(m)) != (n):
        m = m-1
    
    for i in range(1,int(m+1)):
        bp[-i] = -bp[i-1]
        wf[-i] = wf[i-1] 
    return bp, wf


def neumann(inn):
    """
    Returns the weights and nodes for Neumann's numerical integration.

    Args:
        inn (int or numpy.array): Base points (nodes) in the interval [-1;1].

    Returns:
        tuple: A tuple containing:
            - w (numpy.array): Quadrature weights.
            - x (numpy.array): Base points (nodes) in the interval [-1;1].

    Raises:
        TypeError: If the input argument is not an integer.
        ValueError: If there is an error in input dimensions.

    Remarks:
        * 1st N.-method: see Sneeuw (1994) GJI 118, pp 707-716, eq. 19.5.
        * 2nd N.-method: see uberall/GRULE.

    Todo:
        + TypeError is more relevant and shape error from numpy.

    Uses:
        `grule`, `plm`.

    Author:
        Amin Shakya, Interdisciplinary Center for Water Research (ICWaR), Indian Institute of Science (IISc).
    """

    try: #if input is an integer
        x, w = grule(inn)
    except: #if input is an array
        if(len(inn)==1): #2nd Neumann method
            x, w = grule(inn)
            if(np.not_equal(np.mod(x, 1), 0)): #Not integer
                raise TypeError("Integer input argument required")
            
            
        
        elif min(inn.shape) == 1: #1st Neumann method #Size gives 2 outputs for 2d array in matlab; for row and column
            x = inn
            theRAD = np.arccos(x) #x in radian
            l = np.array(list(range(len(x))))
            pp = plm(l, theRAD)
            
            rr = list([2])
            for i in len(x-1):
                rr.append(0)
            r = np.asarray(rr)
                
            w,resid,rank,s = np.linalg.lstsq(pp,r) #Solve system of equations; Double check this operation
            if(x.shape != w.shape):
                w = w.T
            
        else:
            raise ValueError("Error in input dimensions")
            # TO DO: Write more descriptive exception messages
    
    return w, x



def normalklm(lmax: int, typ: str = 'wgs84'):
    """
    Returns an ellipsoidal normal field consisting of normalized -Jn, n=0,2,4,6,8.

    Args:
        lmax (int): Maximum degree of the spherical harmonics.
        typ (str, optional): Ellipsoids can be either 'wgs84' (World Geodetic System 84), 
                             'grs80', or 'he' (hydrostatic equilibrium ellipsoid).

    Returns:
        (numpy.array): Normal field in CS-format (sparse array - [1, -J2, -J4, -J6, -J8]).

    Raises:
        TypeError: If `lmax` is not an integer.
        ValueError: If `lmax` is not positive.
        ValueError: If `typ` is an unknown ellipsoid type. Supports 'wgs84', 'grs80', and 'he'.

    References:
        1. J2, J4 values for hydrostatic equilibrium ellipsoid from Lambeck (1988)
           "Geophysical Geodesy", p.18.

    Author:
        Amin Shakya, Interdisciplinary Center for Water Research (ICWaR), Indian Institute of Science (IISc).
    """
    
    if type(lmax) != int:
        raise TypeError("lmax should be integer")
        
    if lmax < 0:
        raise ValueError("lmax should be positive")
        
    
    typ_ = typ.lower()
    if (typ_ == 'wgs84'):
        J2     =  1.08262982131e-3     #% earth's dyn. form factor (= -C20 unnormalized)
        J4     = -2.37091120053e-6    #% -C40 unnormalized
        J6     =  6.08346498882e-9     #% -C60 unnormalized
        J8     = -1.42681087920e-11    #% -C80 unnormalized
        jcoefs = np.array([1, -J2, -J4, -J6, -J8]).T.reshape(5,1)
        # as lmax + 2 is requires 
        l      = np.arange(0,min(lmax + 2,8 + 2), 2).T
        l.reshape(l.shape[0],1)
        
    elif (typ_ == 'grs80'):
        J2     =  1.08263e-3         # % earth's dyn. form factor (= -C20 unnormalized)
        J4     = -2.37091222e-6     #% -C40 unnormalized
        J6     =  6.08347e-9        #% -C60 unnormalized
        J8     = -1.427e-11         #% -C80 unnormalized
        jcoefs = np.array([1, -J2, -J4, -J6, -J8]).reshape(5,1)
        l      = np.arange(0,min(lmax + 2,8 + 2), 2).T
        l.reshape(l.shape[0],1)
        
    elif ((typ_ == 'he') or (typ_ == 'hydro')):
        J2     = 1.072618e-3		#% earth's dyn. form factor (= -C20 unnormalized)
        J4     = 0.2992e-5     	#% -C40 unnormalized
        jcoefs = np.array([1, -J2, -J4]).T.reshape(5,1)
        # adding (2) beacuse of arange function is only uptp last integer and not including last
        l      = np.arange(0,min(lmax + 2,4 + 2), 2).T
        l.reshape(l.shape[0],1)
        
    else:
        raise ValueError("Unknown type of ellipsoid:   ", typ)
    
    coefs = jcoefs[:len(l)].T / np.sqrt(2*l + 1)
#    coefs.reshape(coefs.shape[0],1)
    
    
    data = np.array(coefs)[0]
    row = np.array(l)
    col = np.zeros(len(l))
    # lmax = 96 then shape=(97, 97) -> consisitent with everything else
    nklm = sparse.coo_matrix((data,(row,col)),shape=(lmax+1,lmax+1)).toarray()
    return nklm


def Gaussian(L: int, cap: int):
    """Generates values for a Gaussian smoothing filter.

    The program delivers the spherical harmonic coefficients of a Gaussian
    smoothing filter. The coefficients are calculated according to Wahr et al. (1998)
    equation (34) and Swenson and Wahr equation (34).

    Args:
        L (int): Maximum degree of the spherical harmonics.
        cap (int): Half width of Gaussian smoothing function [km].

    Returns:
        (np.ndarray): Smoothing coefficients of the Gausiann filter.

    Raises:
        TypeError: If `L` is not an integer.
        ValueError: If `L` is less than or equal to 2.
        TypeError: If `cap` is not an integer.

    References:
        Wahr et al. (1998) equation (34) and Swenson and Wahr equation (34).

    Author:
        Amin Shakya, Interdisciplinary Center for Water Research (ICWaR), Indian Institute of Science (IISc).
    """
    
    #Check input
    if type(L) != int:
        raise TypeError('Degree must be integer')
        
    if L<2:
        raise ValueError('Maximum degree must be higher than 2')
        
    if type(cap) != int:
        raise TypeError('Cap size must be an integer')
        
    #Calculations
    W = np.zeros([L+1, 1])
    b = np.log(2)/(1 - np.cos(cap/6371))
    
    #Recursive calculation of the weighting coefficients
    W[0,0] = 1
    W[1,0] = np.power(10, np.log10( (1 + np.exp(-2*b))/(1-np.exp(-2*b)) - (1/b)))
    
    i = 1
    while i < L:        
        j = i + 1
        W[i+1][0] = W[i-1][0] - (2*(j-1) + 1)/b * W[i][0]
        if W[i+1, 0] > W[i] or W[i+1] < 0:
            W[i+1] = 0
        i = i + 1
    
    return W


def naninterp(X):
    """
    This function uses cubic interpolation to replace NaNs.

    Args:
        X (numpy.array): Array with NaN values.

    Returns:
        numpy.array: Cubic interpolated array.
    """
    
    ok = ~np.isnan(X)
    xp = ok.ravel().nonzero()[0] #Indices of xs with values
    fp = X[~np.isnan(X)]
    
    x  = np.isnan(X).ravel().nonzero()[0] #Indices of xs without values
    
    pchip = PchipInterpolator(xp,fp) #Initialize scipy PHCIP cubic interpolation
    X[np.isnan(X)] = pchip(x) #Interpolate Nan values in X
    
    return X
    