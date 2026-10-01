import os
from pathlib import Path
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
from scipy.stats import rankdata,spearmanr
R=Path(__file__).resolve().parents[3];A=R/'AGE_ANALYSIS';N=A/'next_generation';SEED=20260917
COLORS={'Glioblastoma':'#B65046','Meningioma':'#147D92'}
MODCOLORS={'U':'#cfcfcf','M':'#497ba6','H':'#df9d29','missing':'#eeeeee'}
def table(df,path):df.to_csv(N/path,sep='\t',index=False,na_rep='NA',float_format='%.10g')
def save(x,path):(N/path).write_text(json.dumps(x,indent=2,allow_nan=True)+'\n')
def inv(eta):return np.where(eta<0,21*np.exp(np.minimum(eta,0))-1,21*eta+20)
def inputs():
 p=pd.read_csv(N/'00_audit/frozen_previous/01_tables/age_master_table.tsv',sep='\t');c=pd.read_csv(R/'data/reference/horvath2013/clock353_hg38.tsv',sep='\t');ref=pd.read_csv(R/'data/reference/horvath2013/clock353_sesame_reference.csv').set_index('Probe_ID').loc[c.CpGmarker,'median'].to_numpy();counts=np.stack([np.load(R/f'results/age_hox_joint/tables/{x}_clock353_counts.npz')['counts'] for x in p.sample_id]);return p,c,ref,counts

def corr(x,y):
 a=rankdata(x,axis=-1);b=rankdata(y,axis=-1);a=a-a.mean(axis=-1,keepdims=True);b=b-b.mean(axis=-1,keepdims=True);den=np.sqrt((a*a).sum(axis=-1)*(b*b).sum(axis=-1));return np.divide((a*b).sum(axis=-1),den,out=np.full(np.broadcast_shapes(a.shape[:-1],b.shape[:-1]),np.nan),where=den>0)
def draws(groups,nrep,rng):
 groups=np.asarray(groups);return np.column_stack([rng.choice(np.flatnonzero(groups==g),(nrep,sum(groups==g))) for g in np.unique(groups)])
def association(x,y,groups,adjust=None,nboot=5000,nperm=9999,seed=SEED):
 x=np.array(x,float);y=np.array(y,float);g=np.array(groups);ok=np.isfinite(x)&np.isfinite(y)
 if adjust is not None:ok &=np.isfinite(np.array(adjust,float))
 x=x[ok];y=y[ok];g=g[ok];z=None if adjust is None else np.array(adjust,float)[ok];n=len(x);rng=np.random.default_rng(seed)
 def stat(a,b,z=None):
  if z is None:return corr(a,b)
  a=rankdata(a,axis=-1);b=rankdata(b,axis=-1);z=rankdata(z,axis=-1)
  a-=a.mean(axis=-1,keepdims=True);b-=b.mean(axis=-1,keepdims=True);z-=z.mean(axis=-1,keepdims=True)
  zz=(z*z).sum(axis=-1,keepdims=True);a-=np.divide((a*z).sum(axis=-1,keepdims=True),zz,out=np.zeros_like(zz),where=zz>0)*z;b-=np.divide((b*z).sum(axis=-1,keepdims=True),zz,out=np.zeros_like(zz),where=zz>0)*z
  den=np.sqrt((a*a).sum(axis=-1)*(b*b).sum(axis=-1));return np.divide((a*b).sum(axis=-1),den,out=np.full(np.shape(den),np.nan),where=den>0)
 rho=float(stat(x,y,z));ix=draws(g,nboot,rng);boot=stat(x[ix],y[ix],None if z is None else z[ix]);ci=np.nanquantile(boot,[.025,.975]);loo=np.array([stat(np.delete(x,i),np.delete(y,i),None if z is None else np.delete(z,i)) for i in range(n)])
 if z is None:
  perms=np.array([rng.permutation(n) for _ in range(nperm)]);null=corr(x,y[perms])
 else:
  # Freedman-Lane permutation of rank residuals, nuisance model includes intercept.
  Z=np.column_stack([np.ones(n),rankdata(z)]);H=Z@np.linalg.pinv(Z);rx=rankdata(x);ry=rankdata(y);ex=rx-H@rx;ey=ry-H@ry;perms=np.array([rng.permutation(n) for _ in range(nperm)])
  if len(np.unique(z))==2:
   perms=np.tile(np.arange(n),(nperm,1))
   for value in np.unique(z):
    ix=np.flatnonzero(z==value);perms[:,ix]=np.array([rng.permutation(ix) for _ in range(nperm)])
  yp=(H@ry)[None,:]+ey[perms];res=yp-yp@H.T;null=(res@ex)/np.sqrt((res*res).sum(1)*(ex@ex))
 p=(1+sum(abs(null)>=abs(rho)-1e-12))/(nperm+1)
 return dict(n=n,rho=rho,ci_low=float(ci[0]),ci_high=float(ci[1]),p=p,loo_min=float(np.nanmin(loo)),loo_max=float(np.nanmax(loo)),max_abs_loo_change=float(np.nanmax(abs(loo-rho))),bootstrap_valid=int(np.isfinite(boot).sum())),loo

def bh(p):
 p=np.asarray(p,float);out=np.full(len(p),np.nan);ok=np.flatnonzero(np.isfinite(p));order=ok[np.argsort(p[ok])];out[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(order)/np.arange(1,len(order)+1))[::-1])[::-1]);return out
