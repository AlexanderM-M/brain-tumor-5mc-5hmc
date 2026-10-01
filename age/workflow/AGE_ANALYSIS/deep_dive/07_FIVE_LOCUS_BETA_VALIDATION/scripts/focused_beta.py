"""Focused five-locus beta validation. Reads existing counts; never extracts reads."""
import os
from pathlib import Path
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']: os.environ[k]='1'
import sys,json,hashlib,gzip,argparse
from pathlib import Path
import numpy as np,pandas as pd
from scipy.spatial import cKDTree
from scipy.special import xlogy
R=Path(os.environ['AGE_WORKSPACE'])
D=R/'AGE_ANALYSIS/deep_dive'; O=D/'06_PDR_EPY'
I=D/'.intermediate/local_epialleles'; W=D/'.intermediate/five_locus_beta'
F=D/'07_FIVE_LOCUS_BETA_VALIDATION'
sys.path.insert(0,str(O/'scripts'))
from analyze_local_epialleles import cohort,clock,ref,bh,rare,FIVE
from epiallele_metrics import DIGITS,BINARY_CODE,LOG81
SEED=20260923
METRICS=['PDR','binary_bits','three_state_bits','hidden_bits','M_bits','H_bits','binary_rarefied_bits','three_state_rarefied_bits','PDR_beta_residual','binary_bits_beta_residual','PDR_independent_site_residual']
def save(p,x):p.write_text(json.dumps(x,indent=2,default=lambda v:v.item() if isinstance(v,np.generic) else str(v))+'\n')
def entropy(c):
 n=c.sum(1);p=c/np.maximum(n[:,None],1)
 return -xlogy(p,p).sum(1)/np.log(2)
def features(c):
 n=c.sum(1,dtype=np.uint32);b=np.stack([c[:,BINARY_CODE==k].sum(1) for k in range(16)],axis=1)
 mcode=((DIGITS==1)*[8,4,2,1]).sum(1);hcode=((DIGITS==2)*[8,4,2,1]).sum(1)
 mm=np.stack([c[:,mcode==k].sum(1) for k in range(16)],axis=1)
 hh=np.stack([c[:,hcode==k].sum(1) for k in range(16)],axis=1)
 beta=c@((DIGITS>0).mean(1))/n;mf=c@((DIGITS==1).mean(1))/n;hf=c@((DIGITS==2).mean(1))/n
 marg=(c@((DIGITS>0).astype(float)))/n[:,None]
 pdr=1-(b[:,0]+b[:,15])/n;h2=entropy(b);h3=entropy(c)
 return dict(beta=beta,M_fraction=mf,H_fraction=hf,PDR=pdr,binary_bits=h2,three_state_bits=h3,hidden_bits=h3-h2,M_bits=entropy(mm),H_bits=entropy(hh),PDR_independent_site_residual=pdr-(1-np.prod(marg,axis=1)-np.prod(1-marg,axis=1)))
def setup():
 if (F/'validation/prespecified_protocol.json').exists():raise RuntimeError('Frozen protocol exists; use completed focused caches.')
 for p in [W,F,F/'scripts',F/'tables',F/'figures',F/'validation']:p.mkdir(parents=True,exist_ok=True)
 plan=dict(seed=SEED,scope='Five fixed +/-1kb regions only; existing four-CpG count caches',
 beta_tolerances={'primary':.025,'tight':.01,'wide':.05},primary_depth=20,sensitivity_depth=10,
 matching='Same patient; exact 12-class midpoint gene/CGI context; coverage, quartet span, +/-1kb CpG count ratio 0.8-1.25; beta absolute caliper. Nearest up to20 controls, >=5 distinct 2kb genomic blocks. Exclude all353 clock +/-1kb regions and ENCODE blacklist. No outcome used in selection.',
 candidates='Outcome-blind reference sampling, up to50000 per chromosome stratified across target contexts and coarse target density/span range; fixed seed. Feasibility is conditional on this sampled universe.',
 region='All4CpGs within marker +/-1kb. >=5 matched windows and >=70% eligible windows per patient. >=6 patients for inference. Equal windows within each patient, equal patients.',
 null='9999 matched-control draws. Within each patient, use the same uniform quantile across all overlapping target windows (comonotonic conservative dependence sensitivity); independent uniforms across patients. Empirical centered distributions, two-sided p with +1 correction. Also report fully synchronous patient quantiles as a stronger dependence stress test.',
 inference='5000 patient bootstraps, 95% percentile CI, patient leave-one-out effects and empirical p; BH separately across five loci for each endpoint/scenario; also joint10-test PDR/binary BH.',
 residualization='Control-trained beta spline via cubic truncated powers at0.1,0.25,0.5,0.75,0.9 plus log coverage/span/density and context indicators. Compare target and matched-control residuals. Additionally subtract independent-CpG PDR from each window using its four marginal betas.',
 entropy='Bits and normalized values, exact20/10-molecule rarefaction. Chemical binary M and H indicator entropies are descriptive and not additive; hidden H3-H2 is conditional M/H identity information, not causality.',
 weight='Arbitrary genomic windows have no clock coefficient. Primary test does not fabricate one. Report separate feasibility of other clock regions within absolute coefficient ratio0.8-1.25 and exact central context; do not relax strict calipers.',
 interpretation='Cross-sectional tissue mixture; no clone, selection, evolution or regulatory-function inference.',
 output='Three final TSV tables; <=2 figures; README. Reproducibility files in scripts/validation and private intermediates.')
 save(F/'validation/prespecified_protocol.json',plan)
 # Source-integrity manifests are maintained by the workspace operator.
 sites,off=ref();cl=clock().set_index('CpGmarker');rng=np.random.default_rng(SEED)
 targets={};contexts=set();dens=[];spans=[]
 for loc in FIVE:
  row=cl.loc[loc];ch=int(row.hg38_chromosome[3:]);a,b=map(int,off[ch-1:ch+1]);p=sites[a:b]
  ids=np.flatnonzero((p[:-3]>=row.hg38_start-1000)&(p[3:]<=row.hg38_start+1000))+a
  targets[loc]=ids.tolist();ann=np.load(I/f'annotation_chr{ch}.npz')
  contexts.update(ann['context'][ids-a].tolist());dens.extend(ann['density'][ids-a].tolist());spans.extend(ann['span'][ids-a].tolist())
 save(W/'targets.json',targets)
 mask=np.load(I/'genomic_blacklist_window_mask.npy',mmap_mode='r')
 for ch in range(1,23):
  a,b=map(int,off[ch-1:ch+1]);ann=np.load(I/f'annotation_chr{ch}.npz')
  ctx=ann['context'];den=ann['density'];span=ann['span']
  good=(~ann['excluded_clock'])&(~mask[a:b-3])&(den>=min(dens)*.8)&(den<=max(dens)*1.25)&(span>=min(spans)*.8)&(span<=max(spans)*1.25)
  picks=[]
  for c in sorted(contexts):
   ix=np.flatnonzero(good&(ctx==c));cap=50000//len(contexts)
   if len(ix)>cap:ix=rng.choice(ix,cap,replace=False)
   picks.extend(ix+a)
  tar=np.concatenate([np.array(v,dtype=int) for k,v in targets.items() if int(cl.loc[k].hg38_chromosome[3:])==ch]) if any(int(cl.loc[k].hg38_chromosome[3:])==ch for k in targets) else np.array([],int)
  np.savez_compressed(W/f'candidates_chr{ch}.npz',indices=np.unique(np.r_[picks,tar]).astype(np.uint32))
 print('Protocol and outcome-blind candidate universe frozen',flush=True)
def extract(sample):
 out=W/f'{sample}_focused.npz'
 if out.exists():return
 sites,off=ref();parts=[];counts=[]
 for ch in range(1,23):
  wanted=np.load(W/f'candidates_chr{ch}.npz')['indices'];z=np.load(I/sample/f'chr{ch}.npz')
  idx=z['indices'];j=np.searchsorted(idx,wanted);ok=j<len(idx);ok[ok]&=idx[j[ok]]==wanted[ok]
  cc=z['counts'][j[ok]];ix=wanted[ok];depth=cc.sum(1,dtype=np.uint32);keep=depth>=10;cc=cc[keep];ix=ix[keep];depth=depth[keep]
  ann=np.load(I/f'annotation_chr{ch}.npz');local=ix.astype(int)-int(off[ch-1])
  part=pd.DataFrame(dict(indices=ix,chromosome=ch,position0=sites[ix],n=depth))
  for k in ['span','density','context','excluded_clock']:part[k]=ann[k][local]
  mz=np.load(I/sample/f'metrics_chr{ch}.npz');mi=mz['indices'];jj=np.searchsorted(mi,ix);assert np.all(mi[jj]==ix)
  for dep in [10,20]:
   for name,factor in [('binary',4),('three_state',LOG81)]:part[f'{name}_rarefied_bits{dep}']=mz[f'{name}_rarefied{dep}'][jj]*factor
  parts.append(part);counts.append(cc)
  print(sample,'focused cache',ch,len(cc),flush=True)
 x=pd.concat(parts,ignore_index=True);c=np.concatenate(counts)
 np.savez_compressed(out,counts=c,**{k:x[k].to_numpy() for k in x})
def design(x):
 b=x.beta.to_numpy();v=[np.ones(len(x)),b,b*b,b*b*b]
 v += [np.maximum(b-t,0)**3 for t in [.1,.25,.5,.75,.9]]
 v += [np.log(x[k].to_numpy()) for k in ['n','span','density']]
 v += [(x.context.to_numpy()==k).astype(float) for k in range(1,12)]
 return np.column_stack(v)
def analyze(sample):
 dest=W/f'{sample}_analysis.json'
 if dest.exists():return
 z=np.load(W/f'{sample}_focused.npz');c=z['counts'];x=pd.DataFrame({k:z[k] for k in z.files if k!='counts'})
 for k,v in features(c).items():x[k]=v

 targets=json.loads((W/'targets.json').read_text());mask=np.load(I/'genomic_blacklist_window_mask.npy',mmap_mode='r')
 x['masked']=mask[x.indices.to_numpy()]
 # Only fitting to controls; the five regions never train the beta-response curve.
 base=(~x.excluded_clock)&(~x.masked)&(x.density>0)
 fit=base&(x.n>=20);dd=design(x)
 for name in ['PDR','binary_bits']:
  coef=np.linalg.lstsq(dd[fit],x.loc[fit,name].to_numpy(),rcond=1e-10)[0]
  x[name+'_beta_residual']=x[name]-dd@coef
 records=[];freq=[];nulls={}
 rng=np.random.default_rng(SEED)
 configs=[('primary',20,.025),('tight_beta',20,.01),('wide_beta',20,.05),('depth10',10,.025)]
 for scenario,depth,tol in configs:
  x['binary_rarefied_bits']=x[f'binary_rarefied_bits{depth}'];x['three_state_rarefied_bits']=x[f'three_state_rarefied_bits{depth}']
  pool=x.loc[base&(x.n>=depth)].copy()
  trees={};poolids={}
  for context,g in pool.groupby('context'):
   ids=g.index.to_numpy();vec=np.column_stack([g.beta.to_numpy()/tol]+[np.log(g[k].to_numpy())/np.log(1.25) for k in ['n','span','density']])
   trees[context]=cKDTree(vec);poolids[context]=ids
  for loc in FIVE:
   target=x.loc[x.indices.isin(targets[loc])&(x.n>=depth)&(~x.masked)]
   pairs=[];distances=[]
   for idx,t in target.iterrows():
    ctx=t.context
    if ctx not in trees:continue
    v=np.r_[t.beta/tol,np.log(t[['n','span','density']].to_numpy(float))/np.log(1.25)]
    sel=trees[ctx].query_ball_point(v,1+1e-12,p=np.inf)
    if len(sel)<5:continue
    ids=poolids[ctx][sel];g=x.loc[ids]
    dv=np.column_stack([(g.beta.to_numpy()-t.beta)/tol]+[np.log(g[k].to_numpy()/t[k])/np.log(1.25) for k in ['n','span','density']])
    order=np.argsort(np.sum(dv*dv,axis=1),kind='stable');chosen=[];blocks=set()
    for kk in order:
     ii=ids[kk];row=x.loc[ii];block=(int(row.chromosome),int(row.position0)//2000)
     if block in blocks:continue
     blocks.add(block);chosen.append(ii)
     if len(chosen)==20:break
    if len(chosen)>=5:
     pairs.append((int(idx),np.array(chosen,int)))
     distances.extend(np.abs(x.loc[chosen,'beta'].to_numpy()-t.beta).tolist())
   eligible=len(target);passed=len(pairs)>=5 and len(pairs)>=.7*eligible
   rec=dict(sample_id=sample,locus=loc,scenario=scenario,eligible_windows=eligible,matched_windows=len(pairs),patient_eligible=passed,
    mean_abs_beta_mismatch=float(np.mean(distances)) if distances else None,
    median_controls=float(np.median([len(v) for _,v in pairs])) if pairs else 0)
   records.append(rec)
   if not passed:continue
   # Save selected reference indices for complete audit; no fresh molecule extraction.
   np.savez_compressed(W/f'{sample}_{loc}_{scenario}_matches.npz',
    targets=np.array([int(x.loc[j,'indices']) for j,_ in pairs],np.uint32),
    controls=np.array([np.pad(x.loc[v,'indices'].to_numpy(np.uint32),(0,20-len(v)),constant_values=2**32-1) for _,v in pairs]))
   tid=np.array([j for j,_ in pairs]);flat=np.unique(np.concatenate([v for _,v in pairs]))
   for field in METRICS:
    obs=x.loc[tid,field].to_numpy();means=np.array([x.loc[v,field].mean() for _,v in pairs])
    # Quantile-coupled control draws preserve maximal positive dependence among overlapping windows.
    qq=np.stack([np.sort(x.loc[v,field].to_numpy())[np.minimum((np.arange(1000)+.5)*len(v)//1000,len(v)-1).astype(int)] for _,v in pairs]).mean(0)
    qq-=qq.mean()
    rec[field]=float(np.mean(obs-means));rec[field+'_target']=float(np.mean(obs));rec[field+'_control']=float(np.mean(means))
    nulls[f'{loc}|{scenario}|{field}']=qq.astype(np.float32)
   rec['beta_target']=float(x.loc[tid,'beta'].mean())
   rec['beta_control']=float(np.mean([x.loc[v,'beta'].mean() for _,v in pairs]))
   for name in ['M_fraction','H_fraction']:
    rec[name]=float(x.loc[tid,name].mean());rec[name+'_control']=float(np.mean([x.loc[v,name].mean() for _,v in pairs]))
   if scenario=='primary':
    for alphabet,code,ncat in [('binary',BINARY_CODE,16),('chemical',np.arange(81),81)]:
     for arm in ['target','control']:
      vectors=[]
      for j,v in pairs:
       use=np.array([j]) if arm=='target' else v
       cc=c[use];pp=cc/cc.sum(1)[:,None]
       vv=np.array([pp[:,code==k].sum(1).mean() for k in range(ncat)])
       vectors.append(vv)
      avg=np.mean(vectors,axis=0)
      assert np.isclose(avg.sum(),1)
      for k,value in enumerate(avg):
       pattern=format(k,'04b') if alphabet=='binary' else ''.join('0MH'[q] for q in DIGITS[k])
       freq.append(dict(sample_id=sample,locus=loc,arm=arm,alphabet=alphabet,pattern=pattern,frequency=float(value),matched_windows=len(pairs)))
 np.savez_compressed(W/f'{sample}_nulls.npz',**nulls)
 pd.DataFrame(freq).to_csv(W/f'{sample}_frequencies.tsv',sep='\t',index=False)
 save(dest,records);print(sample,'analysis complete',flush=True)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('step',choices=['setup','sample']);ap.add_argument('--sample')
 a=ap.parse_args()
 if a.step=='setup':setup()
 else:extract(a.sample);analyze(a.sample)
if __name__=='__main__':main()
