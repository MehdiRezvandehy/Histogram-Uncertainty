import math  
import numpy as np 
from numba import jit 
import scipy
import scipy.linalg
import random
import pandas as pd
from scipy.stats import norm
def Spatial_Bootstrap(x,y,var,variog,seed,nsim,outfl):
    """ Spatial Bootstrap for 2D uncertainty quantificaion of spatial data. 
    x     : x coordinate (numpy array) 
    y     : y coordinate (numpy array)
    var   : variable (numpy array)
    variog: variogram parameter (GSLIB format)
    nsim  : number of simulaton
    
    output:
    1- yvec_r : uncertainty in distribution of data
    2- mean   : uncertainty in distribution of mean
    """

    random.seed(seed)
    # load the variogram fo data
    nst = variog['nst']
    cc = np.zeros(nst); aa = np.zeros(nst); it = np.zeros(nst)
    ang = np.zeros(nst); anis = np.zeros(nst)
    
    c0 = variog['nug'];
    cc[0] = variog['cc1']; it[0] = variog['it1']; ang[0] = variog['azi1'];
    aa[0] = variog['hmaj1']; anis[0] = variog['hmin1']/variog['hmaj1'];
    if nst == 2:
        cc[1] = variog['cc2']; it[1] = variog['it2']; ang[1] = variog['azi2'];
        aa[1] = variog['hmaj2']; anis[1] = variog['hmin2']/variog['hmaj2'];
    
    # Initialize the covariance subroutine and cbb at the same time:
    rotmat, maxcov = setup_rotmat2(c0,nst,it,cc,ang)
    PMX = 9999.0
    cbb = cova2(0.0,0.0,0.0,0.0,nst,c0,PMX,cc,aa,it,ang,anis,rotmat,maxcov)
    
    # Allocate the needed memory:
    len_=len(x)
    corr=np.zeros((len_,len_))
    low=np.zeros((len_,len_))
    wvec=np.zeros((len_,nsim))
    yvec_n=np.zeros((len_,nsim))
    
    for i in range(len_):
        for j in range(len_):
            co=cova2(x[i], y[i], x[j], y[j], nst, c0, PMX, cc, aa, it, ang, anis, rotmat, maxcov)
            if(i==j): co=co+0.01 # Add small value to diagonal elements to make matrix positive definite
            corr[i,j] = co
            corr[j,i] = co
    
    # LU decomposition
    try: 
        L=scipy.linalg.cholesky(corr, lower=True, overwrite_a=True)
    except ValueError:
        pass
    
    # Generate uncorrelated Gaussian simulation at data locations
    mu=0; sigma=1
    for i in range(len_):
        Dist = np.random.normal(mu, sigma, nsim)
        wvec[i,:]=Dist
    
    # Generate correlated Gaussian simulation at data locations
    yvec_r=[]
    for k in range(nsim):
        tmp=(np.matmul(L,wvec[:,k]))
        tmp_=[]
        if(outfl):
            outfl_=outfl+'_'+str(k+1)
            f1=open(outfl_,'w')
            txt='data \n'+'1 \n'+'Value \n'
        f1.write(str(txt))
        for j in range(len(tmp)):
            yvec_n[j,k]=tmp[j]
            prob=norm.cdf(tmp[j])
            quantle=np.quantile(var, prob, axis=0, keepdims=True)[0]
            if(outfl): f1.write("{0:.4f}".format(quantle)+"\n")
            tmp_.append(quantle)
        f1.close()    
        yvec_r.append(tmp_)    
    
    # Calculate distribution of the mean
    mean=[]
    for k in range(nsim):
        mean.append(np.mean(yvec_r[k]))
                    
    print(' Effective number of data= ', np.var(var)/np.var(mean))    
    print(' Mean of mean=',np.round(np.mean(mean),4),', Variance of the mean = ', np.round(np.var(mean),4))    
    
    return yvec_r, mean


def setup_rotmat2(c0,nst,it,cc,ang):
    DTOR=3.14159265/180.0; EPSLON=0.000000; PI=3.141593
# The first time around, re-initialize the cosine matrix for the
# variogram structures:
    rotmat = np.zeros((4,nst))
    maxcov = c0
    for js in range(0,nst):
        azmuth = (90.0-ang[js])*DTOR
        rotmat[0,js] =  math.cos(azmuth)
        rotmat[1,js] =  math.sin(azmuth)
        rotmat[2,js] = -1*math.sin(azmuth)
        rotmat[3,js] =  math.cos(azmuth)
        if it[js] == 4:
            maxcov = maxcov + 9999.9
        else:
            maxcov = maxcov + cc[js]
    return rotmat, maxcov

@jit(nopython=True)

def cova2(x1, y1, x2, y2, nst, c0, pmx, cc, aa, it, ang, anis, rotmat, maxcov):
    """The reimplementation is by Michael Pyrcz, Associate Professor,
    the University of Texas at Austin."""
    """Calculate the covariance associated with a variogram model specified by
    a nugget effect and nested variogram structures.
    :param x1: x coordinate of first point
    :type x1: float
    :param y1: y coordinate of first point
    :type y1: float
    :param x2: x coordinate of second point
    :type x2: float
    :param y2: y coordinate of second point
    :type y2: float
    :param nst: number of nested structures (maximum of 4)
    :type nst: int
    :param c0: isotropic nugget constant (TODO: not used)
    :type c0: float
    :param pmx: Maximum variogram value needed for kriging when using power
                model. pmx is a unique value used for all nested structures
                that use the power model, so pmx should be chosen to account
                for the largest structure that uses the power model.
    :type pmx: float
    :param cc: multiplicative factor of each nested structure
    :type cc: array
    :param aa: parameter `a` of each nested structure
    :type aa: array
    :param it: Integer value indicating type of variogram model
             for values 0,1,2,..., nst
             it[value] == 1: Spherical model 
                 (aa[value] == `a` is the range, cc[value] is the contribution)
             it[value] == 2: Exponential model 
                 (aa[value] == `a`, 3a is the practical range, 
                 cc[value] is the contribution)
             it[value] == 3: Gaussian model 
                 (aa[value] == `a`, a*sqrt(3)  is the practical range), 
                 cc[value] is the contribution)
             it[value] == 4: Power model 
                 (aa[value] == `a` is the power such that 0 < a < 2,
                 if linear, then a == 1, and cc[value] is the slope)
    :type it: array
    :param ang: azimuth angle measured in degrees clockwise from positive 
                y-diretion for each variogram structure: not used
                (accounted for in anis)
    :type ang: array
    :param anis: Anistropy factors that apply after rotations
    :type anis: array
    :param rotmat: rotation matrices
    :type rotmat: array
    :param maxcov: maximum covariance value
    :type maxcov: float
    :return: covariance of a nested variogram model described by the inputs
    :type return: float
    """
    EPSLON = 0.000001

    # Check for very small distance
    dx = x2 - x1
    dy = y2 - y1

    if (dx * dx + dy * dy) < EPSLON:
        cova2_ = maxcov
        return cova2_

    # Non-zero distance, loop over all the structures
    cova2_ = 0.0
    for js in range(0, nst):
        # Compute the appropriate structural distance
        dx1 = dx * rotmat[0, js] + dy * rotmat[1, js]
        dy1 = (dx * rotmat[2, js] + dy * rotmat[3, js]) / anis[js]
        h = math.sqrt(max((dx1 * dx1 + dy1 * dy1), 0.0))
        if it[js] == 1:
            # Spherical model
            hr = h / aa[js]
            if hr < 1.0:
                cova2_ = cova2_ + cc[js] * (1.0 - hr * (1.5 - 0.5 * hr * hr))
        elif it[js] == 2:
            # Exponential model
            cova2_ = cova2_ + cc[js] * np.exp(-3.0 * h / aa[js])
        elif it[js] == 3:
            # Gaussian model
            hh = -3.0 * (h * h) / (aa[js] * aa[js])
            cova2_ = cova2_ + cc[js] * np.exp(hh)
        elif it[js] == 4:
            # Power model
            cov1 = pmx - cc[js] * (h ** aa[js])
            cova2_ = cova2_ + cov1
    return cova2_

def GSLIB_View(File):
    """a fucntion to convert ASCII flat file compatible with Geo-EAS to Pandas Data Frame"""
    with open(File, 'r') as f:
        next(f)
        txt=[] 
        tmp_val=[]
        value=[]
        ii=0
        for line in f:
            if(ii==0):
                p = line.split()
                no=int(p[0])
            elif ii<=no:
                p = line.replace('\n', '').split()
                txt.append(p[0])  
            else: 
                tmp=[]
                p = line.split()
                for j in range(len(txt)):
                    tmp.append(float(p[j]))  
                tmp_val.append(tmp)
            ii+=1 
        for j in range(len(txt)):
            tmp=[]
            for k in range(len(tmp_val)):
                tmp.append(tmp_val[k][j])
            value.append(tmp)            
    f.close()
    pd_=["" for x in range(len(txt))]
    pd_=[]
    for j in range(len(txt)):
            pd_.append(pd.DataFrame(({txt[j]:value[j]})))  
    result = pd.concat(pd_, axis=1)
    return result