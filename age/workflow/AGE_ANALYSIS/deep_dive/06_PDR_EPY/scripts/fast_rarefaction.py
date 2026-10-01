"""Numerically stable sparse exact hypergeometric entropy contributions."""
import os
from pathlib import Path
import numpy as np
from scipy.special import gammaln
def contributions(n,counts,depth):
    k=np.asarray(counts,dtype=np.float64)[:,None]
    j=np.arange(1,depth+1,dtype=np.float64)[None,:]
    def choose(a,b):
        return gammaln(a+1)-gammaln(b+1)-gammaln(a-b+1)
    valid=(k>=j)&(n-k>=depth-j)
    with np.errstate(invalid='ignore',over='ignore'):
        lp=choose(k,j)+choose(n-k,depth-j)-choose(float(n),float(depth))
        probability=np.where(valid,np.exp(lp),0)
    value=(probability*(-(j/depth)*np.log2(j/depth))).sum(1)
    value[(counts==0)|(counts==n)]=0
    return value
