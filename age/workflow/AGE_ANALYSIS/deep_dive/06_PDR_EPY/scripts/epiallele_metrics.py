"""Four consecutive reference-CpG windows; U=0, M=1, H=2."""
import os
from pathlib import Path
import numpy as np
from scipy.special import xlogy
CODES = np.arange(81)
DIGITS = (CODES[:,None] // np.array([27,9,3,1])) % 3
BINARY_CODE = ((DIGITS>0)*np.array([8,4,2,1])).sum(1)
LOG81 = np.log2(81)

def counts_from_calls(indices, states):
    indices=np.asarray(indices); states=np.asarray(states)
    if len(indices)!=len(states): raise ValueError("Length mismatch")
    if len(indices)>1 and not np.all(np.diff(indices)>0): raise ValueError("Indices must be strictly increasing")
    result={}
    for j in range(len(indices)-3):
        z=states[j:j+4]
        if indices[j+3]-indices[j]!=3 or np.any(z<0): continue
        if np.any(z>2): raise ValueError("Invalid state")
        code=int(z@np.array([27,9,3,1]))
        result.setdefault(int(indices[j]),np.zeros(81,np.uint32))[code]+=1
    return result

def metrics(counts):
    c=np.asarray(counts)
    if c.shape[-1]!=81 or np.any(c<0): raise ValueError("Need 81 nonnegative counts")
    if np.any(c!=np.floor(c)): raise ValueError("Counts must be integral")
    c=c.astype(np.float64)
    b=np.stack([c[...,BINARY_CODE==k].sum(-1) for k in range(16)],axis=-1)
    n=c.sum(-1); safe=np.maximum(n,1)
    def entropy(a):
        p=a/safe[...,None]
        return -xlogy(p,p).sum(-1)/np.log(2)
    h2=entropy(b);h3=entropy(c)
    pdr=(n-b[...,0]-b[...,15])/safe
    out=dict(informative_molecules=n,PDR=pdr,binary_entropy=h2/4,
             three_state_entropy=h3/LOG81,DeltaEntropy=h3/LOG81-h2/4,
             hidden_information_bits=h3-h2)
    for k in out:
        if k!='informative_molecules':out[k]=np.where(n>0,out[k],np.nan)
    return out

def rarefaction_lookup(n, depth):
    """Exact expected plug-in entropy contributions for sampling without replacement."""
    from scipy.stats import hypergeom
    k=np.arange(n+1)[:,None];j=np.arange(1,depth+1)[None,:]
    terms=-(j/depth)*np.log2(j/depth)
    return (hypergeom.pmf(j,n,k,depth)*terms).sum(1)

def rarefied_entropy(counts,depth=20):
    c=np.asarray(counts);n=c.sum(1).astype(int)
    b=np.stack([c[:,BINARY_CODE==k].sum(1) for k in range(16)],axis=1)
    h2=np.full(len(c),np.nan);h3=h2.copy()
    for nn in np.unique(n[n>=depth]):
        mask=n==nn; lut=rarefaction_lookup(int(nn),depth)
        h2[mask]=lut[b[mask].astype(int)].sum(1)/4
        h3[mask]=lut[c[mask].astype(int)].sum(1)/LOG81
    return h2,h3
