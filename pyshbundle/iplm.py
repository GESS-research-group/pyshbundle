# % IPLM Integrals of the fully normalized associated Legendre functions
# % over blocks for a selected order M. 
# %
# % HOW: p = iplm(l,m,theRAD)		- assumes dt = theRAD(2)-theRAD(1)
# %      p = iplm(l,m,theRAD,dt)
# %
# % IN:
# %    l ........ degree (vector). Integer, but not necessarily monotonic.
# %               For l < m a vector of zeros will be returned.
# %    m ........ order (scalar)
# %    theRAD ... co-latitude [rad] (vector)
# %    dt ....... integration block-size [rad] (scalar). Default: dt = theRAD(2)-theRAD(1)
# %
# % OUT: 
# %    p ........ Matrix with integrated Legendre functions.
# %               Functions are integrated from theRAD(i)-dt/2 till theRAD(i)+dt/2.
# %               The matrix has length(TH) rows and length(L) columns, unless L 
# %               or TH is scalar. Then the output vector follows the shape of 
# %               respectively L or TH. 
# % 
# % USES:
# %    PLM
# %
# % REMARKS:
# %    The blocks at the pole might become too large under circumstances.
# %    This is not treated separately, i.e. unwanted output may appear.
# %    In case TH is scalar, dt will be 1 (arbitrarily).

# % -------------------------------------------------------------------------
# % project: SHBundle 
# % -------------------------------------------------------------------------
# % authors:
# %    Dimitris TSOULIS (DT), IAPG, TU-Munich  
# %    Nico SNEEUW (NS), IAPG, TU-Munich
# %    Markus ANTONI (MA), GI, Uni Stuttgart
# %    <bundle@gis.uni-stuttgart.de>
# % -------------------------------------------------------------------------
# % revision history:
# %    2012-01-23: MR, input of plm in radian
# %    1999-02-??: NS, brush-up (layout, help text, inactive lines,...)
# %                    output variable (column extraction)
# %                    Plm's BEFORE for-loop -> large speed-up
# %                    variable redefinition (e.g. loop variable l)
# %    1998-12-??: DT, initial version 
# %----------------------------------------------------------------------------
# % license:
# %    This program is free software; you can redistribute it and/or modify
# %    it under the terms of the GNU General Public License as published by
# %    the  Free  Software  Foundation; either version 3 of the License, or
# %    (at your option) any later version.
# %  
# %    This  program is distributed in the hope that it will be useful, but 
# %    WITHOUT   ANY   WARRANTY;  without  even  the  implied  warranty  of 
# %    MERCHANTABILITY  or  FITNESS  FOR  A  PARTICULAR  PURPOSE.  See  the
# %    GNU General Public License for more details.
# %  
# %    You  should  have  received a copy of the GNU General Public License
# %    along with Octave; see the file COPYING.  
# %    If not, see <http://www.gnu.org/licenses/>.
# % ----------------------------------------------------------------------------

# % diagnostics and preliminaries
# %narginchk(3, 4) % error(nargchk(3,4,nargin))

import numpy as np
def iplm(l,m,theRAD,dt=-9999):
    import numpy as np
    import sys
    from . import plm
    if dt==-9999:
        if len(theRAD) == 1:
            dt = np.pi/180;
        else:
            dt = theRAD[1] - theRAD[0]
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
    # if min(dt.shape) !=1:
    #     print('Block size DT must be scalar.') 
        sys.exit([])
    if dt == 0: 
        print('DT cannot be zero') 
        sys.exit([])
        
    # init
    lcol = len(l)
    trow = len(theRAD)
    n = len(theRAD)
    theRAD.T
    if min(theRAD) < 0 or max(theRAD) > np.pi:
        print('Is the co-latitude ''theta'' given in radian?')
        sys.exit([])
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
    ptmp00 = np.cos(theRAD - dt/2) - ctplus;
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


# The integrated functions have been computed. What remains to be done, is to
# extract the proper columns from ptmp, corresponding to the vector lvec. 
# If l or theta is scalar the output matrix p reduces to a vector. It should
# have the shape of respectively theta or l in that case.

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