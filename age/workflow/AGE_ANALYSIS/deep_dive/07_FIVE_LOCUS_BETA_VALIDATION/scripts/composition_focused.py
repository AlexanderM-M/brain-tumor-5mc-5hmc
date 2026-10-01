"""Focused pattern decomposition; counts and matching are unchanged."""
import os
from pathlib import Path
from focused_beta import *
EXTRA=['PDR_M_only','PDR_H_only','PDR_MH_mixed','fully_modified_MH_mixed','discordant_top2_share','discordant_entropy_normalized']
def extra(c):
 p=c/c.sum(1)[:,None];mod=DIGITS>0;dis=(mod.sum(1)>0)&(mod.sum(1)<4);m=(DIGITS==1).any(1);h=(DIGITS==2).any(1)
 out={name:p[:,sel].sum(1) for name,sel in [
 ('PDR_M_only',dis&m&~h),('PDR_H_only',dis&h&~m),('PDR_MH_mixed',dis&m&h),
 ('fully_modified_MH_mixed',mod.all(1)&m&h)]}
 b=np.stack([c[:,BINARY_CODE==k].sum(1) for k in range(1,15)],axis=1)
 n=b.sum(1);q=np.divide(b,n[:,None],out=np.zeros_like(b,dtype=float),where=n[:,None]>0)
 out['discordant_top2_share']=np.sort(q,axis=1)[:,-2:].sum(1)
 out['discordant_entropy_normalized']=-xlogy(q,q).sum(1)/np.log(14)
 out['has_discordance']=n>0
 return out
def main():
 for s in cohort().sample_id:
  dest=W/f'{s}_composition.json'
  if dest.exists():continue
  z=np.load(W/f'{s}_focused.npz');idx=z['indices'];c=z['counts'];order=np.argsort(idx);rows=[];nulls={}
  records=json.loads((W/f'{s}_analysis.json').read_text())
  for rec in records:
   if rec['scenario']!='primary' or not rec['patient_eligible']:continue
   loc=rec['locus'];pairs=np.load(W/f'{s}_{loc}_primary_matches.npz');tar=pairs['targets'];ct=pairs['controls']
   ids=np.unique(np.r_[tar,ct.ravel()]);ids=ids[ids!=2**32-1];jj=order[np.searchsorted(idx[order],ids)];assert np.all(idx[jj]==ids)
   f=extra(c[jj]);pos={int(k):j for j,k in enumerate(ids)}
   tid=np.array([pos[int(k)] for k in tar]);cid=[np.array([pos[int(k)] for k in v if k!=2**32-1]) for v in ct]
   row=dict(sample_id=s,locus=loc,scenario='primary')
   for name in EXTRA:
    if name.startswith('discordant'):
     # Conditional entropy/top2 is only defined where discordant molecules exist.
     keep=f['has_discordance'][tid]&np.array([f['has_discordance'][v].any() for v in cid])
     tt=tid[keep];cv=[v[f['has_discordance'][v]] for j,v in enumerate(cid) if keep[j]]
    else:tt=tid;cv=cid
    obs=f[name][tt];controls=[f[name][v] for v in cv]
    row[name+'_target']=float(np.mean(obs));row[name+'_control']=float(np.mean([v.mean() for v in controls]))
    row[name]=row[name+'_target']-row[name+'_control']
    qq=np.stack([np.sort(v)[np.minimum(((np.arange(1000)+.5)*len(v)//1000).astype(int),len(v)-1)] for v in controls]).mean(0);qq-=qq.mean()
    nulls[f'{loc}|primary|{name}']=qq.astype(np.float32)
   rows.append(row)
  save(dest,rows);np.savez_compressed(W/f'{s}_composition_nulls.npz',**nulls)
  print(s,'composition complete',flush=True)
if __name__=='__main__':main()
