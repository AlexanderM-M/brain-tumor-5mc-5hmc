"""Genome-wide local-window summaries and patient-level inference only."""
import os
from pathlib import Path
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import sys,json,gzip,hashlib,functools,fcntl,inspect
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import spearmanr,rankdata
from scipy.spatial import cKDTree
from epiallele_metrics import metrics,rarefaction_lookup,LOG81,BINARY_CODE
from fast_rarefaction import contributions
from local_epialleles import R,D,O,I,save,SEED
FIELDS=['PDR','binary_entropy','three_state_entropy','DeltaEntropy','hidden_information_bits']
TEST_FIELDS=['PDR','binary_entropy','three_state_entropy','hidden_information_bits',
             'binary_rarefied','three_state_rarefied','hidden_rarefied_bits']
FIVE=['cg07158339','cg22947000','cg01511567','cg19761273','cg13460409']
def write(x,p):x.to_csv(p,sep='\t',index=False,na_rep='NA',float_format='%.10g')
def cohort():return pd.read_csv(D/'FINAL_INTEGRATION/tables/paired_patient_clock_WGS.tsv',sep='\t')
def ref():return [np.load(R/'data/reference'/p,mmap_mode='r') for p in ['cpg_positions.npy','cpg_offsets.npy']]
def clock():return pd.read_csv(R/'data/reference/horvath2013/clock353_hg38.tsv',sep='\t')
def interval_mask(pos, intervals):
    change=np.zeros(len(pos)+1,np.int32)
    for a,b in intervals:
        j,k=np.searchsorted(pos,[a,b]);change[j]+=1;change[k]-=1
    return np.cumsum(change[:-1])>0
def annotate():
    sites,off=ref();cl=clock();records=[]
    for ch in range(1,23):
        p=I/f'annotation_chr{ch}.npz'
        if p.exists():continue
        a,b=map(int,off[ch-1:ch+1]);pos=sites[a:b].astype(np.int64)
        midpoint=(pos[:-3]+pos[3:])//2
        span=pos[3:]-pos[:-3]+2
        density=np.searchsorted(pos,midpoint+1001)-np.searchsorted(pos,midpoint-1000)
        rr=json.loads(gzip.decompress((R/f'data/reference/ucsc_refseq_snapshot/chr{ch}.json.gz').read_bytes()))['ncbiRefSeq']
        body=interval_mask(midpoint,[(t['txStart'],t['txEnd']) for t in rr])
        prom=[]
        for t in rr:
            tss=t['txStart'] if t['strand']=='+' else t['txEnd']-1
            prom.append((max(0,tss-2000),tss+501) if t['strand']=='+' else (max(0,tss-500),tss+2001))
        gene=body.astype(np.uint8);gene[interval_mask(midpoint,prom)]=2
        xx=json.loads(gzip.decompress((R/f'data/reference/story_context_snapshot/cpgIslandExt_chr{ch}.json.gz').read_bytes()))
        assert xx['start']==0 and xx['end']>=int(pos[-1]),'Incomplete island annotation'
        islands=[(z['chromStart'],z['chromEnd']) for z in xx['cpgIslandExt']]
        island=np.zeros(len(midpoint),np.uint8)
        for width,code in [(4000,1),(2000,2),(0,3)]:
            island[interval_mask(midpoint,[(max(0,x-width),y+width) for x,y in islands])]=code
        c=cl.loc[cl.hg38_chromosome==f'chr{ch}']
        # Exact interval overlap, rather than a midpoint-only exclusion.
        excluded=np.zeros(len(midpoint),bool)
        for t in c.itertuples():
            excluded|=(pos[:-3]<=int(t.hg38_start)+1000)&(pos[3:]>=int(t.hg38_start)-1000)
        np.savez_compressed(p,span=span.astype(np.uint32),density=density.astype(np.uint32),
                            context=(gene*4+island),excluded_clock=excluded)
        records.append(dict(chromosome=f'chr{ch}',reference_windows=len(midpoint)))
    print('Autosomal genomic annotations ready',flush=True)
@functools.lru_cache(maxsize=1500)
def lut(n,depth):return rarefaction_lookup(n,depth)
def rare(c,depth):
    n=c.sum(1).astype(int)
    b=np.stack([c[:,BINARY_CODE==k].sum(1) for k in range(16)],axis=1)
    h2=np.full(len(c),np.nan);h3=h2.copy()
    for nn in np.unique(n[n>=depth]):
        mask=n==nn
        if nn<=512:
            v=lut(int(nn),depth)
            h2[mask]=v[b[mask].astype(int)].sum(1)/4
            h3[mask]=v[c[mask].astype(int)].sum(1)/LOG81
        else:
            cc=c[mask];bb=b[mask]
            keys,inv=np.unique(np.r_[cc.ravel(),bb.ravel()],return_inverse=True)
            v=contributions(int(nn),keys,depth)
            h3[mask]=v[inv[:cc.size]].reshape(cc.shape).sum(1)/LOG81
            h2[mask]=v[inv[cc.size:]].reshape(bb.shape).sum(1)/4
    return h2,h3
def _prepare_sample(sample):
    assert (I/sample/'complete.json').exists(),'Extraction incomplete'
    for ch in range(1,23):
        p=I/sample/f'metrics_chr{ch}.npz'
        if p.exists():continue
        z=np.load(I/sample/f'chr{ch}.npz');n=z['counts'].sum(1,dtype=np.uint32)
        keep=n>=10;ix=z['indices'][keep];c=z['counts'][keep]
        arrays={k:np.empty(len(c),np.float32) for k in FIELDS}
        for depth in [10,20]:
            arrays[f'binary_rarefied{depth}']=np.empty(len(c),np.float32)
            arrays[f'three_state_rarefied{depth}']=np.empty(len(c),np.float32)
            arrays[f'hidden_rarefied_bits{depth}']=np.empty(len(c),np.float32)
        for a in range(0,len(c),40000):
            b=min(a+40000,len(c));v=metrics(c[a:b])
            for k in FIELDS:
                value=v[k]
                if k=='hidden_information_bits':value=np.where(np.abs(value)<1e-12,0,np.maximum(0,value))
                arrays[k][a:b]=value
            for depth in [10,20]:
                h2,h3=rare(c[a:b],depth)
                arrays[f'binary_rarefied{depth}'][a:b]=h2
                arrays[f'three_state_rarefied{depth}'][a:b]=h3
                hidden=h3*LOG81-h2*4
                arrays[f'hidden_rarefied_bits{depth}'][a:b]=np.where(np.abs(hidden)<1e-12,0,np.maximum(0,hidden))
        tmp=p.with_suffix('.tmp.npz')
        np.savez_compressed(tmp,indices=ix,informative_molecules=n[keep],**arrays);tmp.replace(p)
        print(sample,'metrics',ch,flush=True)
    return sample
def prepare_sample(sample):
    folder=I/sample
    lock=(folder/'metrics.lock').open('w')
    fcntl.flock(lock,fcntl.LOCK_EX)
    source=inspect.getsource(_prepare_sample)+inspect.getsource(rare)+inspect.getsource(lut.__wrapped__)+inspect.getsource(metrics)+inspect.getsource(contributions)
    sig=hashlib.sha256(source.encode()).hexdigest()
    meta=folder/'metrics_signature.json'
    if meta.exists():
        previous=json.loads(meta.read_text())['sha256']
        if previous!=sig:
            transition=json.loads((O/'validation/optimization_transition.json').read_text())
            assert previous in transition['accepted_legacy_metric_signatures'], 'Stale local metric cache'
            assert json.loads((O/'validation/sparse_rarefaction_equivalence.json').read_text())['passed']
            save(meta,dict(sha256=sig,validated_equivalent_previous_sha256=previous))
    else:save(meta,dict(sha256=sig))
    result=_prepare_sample(sample)
    save(folder/'metrics_complete.json',dict(sample_id=sample,complete=True,sha256=sig))
    return result

def common_panel():
    sites,off=ref();samples=cohort().sample_id
    totals={10:{20:0,16:0},20:{20:0,16:0}}
    for ch in range(1,23):
        a,b=map(int,off[ch-1:ch+1]);support=np.zeros((2,b-a),np.uint8)
        for sample in samples:
            z=np.load(I/sample/f'metrics_chr{ch}.npz');ix=z['indices'].astype(int)-a;n=z['informative_molecules']
            for j,depth in enumerate([10,20]):support[j,ix[n>=depth]]+=1
        np.savez_compressed(I/f'support_chr{ch}.npz',support=support)
        for j,depth in enumerate([10,20]):
            for gate in [16,20]:totals[depth][gate]+=int((support[j]>=gate).sum())
    panels={depth:dict(required_patients=20 if totals[depth][20]>=500 else 16,
                       available_windows=totals[depth][20] if totals[depth][20]>=500 else totals[depth][16],
                       complete_common_windows=totals[depth][20],
                       minimum_patient_windows=500 if totals[depth][20]>=500 else 1000)
            for depth in [10,20]}
    save(O/'validation/coverage_panel.json',panels);return panels
def metric_frame(sample):
    sites,off=ref();parts=[]
    for ch in range(1,23):
        z=np.load(I/sample/f'metrics_chr{ch}.npz')
        x=pd.DataFrame({k:z[k] for k in z.files})
        idx=x['indices'].to_numpy().astype(int)-int(off[ch-1])
        ann=np.load(I/f'annotation_chr{ch}.npz')
        for k in ['span','density','context','excluded_clock']:x[k]=ann[k][idx]
        x['chromosome']=ch;x['position0']=sites[x['indices'].to_numpy()]
        sup=np.load(I/f'support_chr{ch}.npz')['support']
        x['support10']=sup[0,idx];x['support20']=sup[1,idx]
        parts.append(x)
    return pd.concat(parts,ignore_index=True)
def bh(p):
    p=np.asarray(p,float);out=np.full(len(p),np.nan);ix=np.flatnonzero(np.isfinite(p))
    order=ix[np.argsort(p[ix])];n=len(order)
    if n:out[order]=np.minimum(1,np.minimum.accumulate((p[order]*n/np.arange(1,n+1))[::-1])[::-1])
    return out
def correlation(x,y,cov=None):
    rx=rankdata(x);ry=rankdata(y)
    if cov is not None:
        design=np.column_stack([np.ones(len(x))]+[rankdata(cov[:,j]) for j in range(cov.shape[1])])
        rx=rx-design@np.linalg.lstsq(design,rx,rcond=None)[0]
        ry=ry-design@np.linalg.lstsq(design,ry,rcond=None)[0]
    if np.std(rx)<1e-12 or np.std(ry)<1e-12:return np.nan
    return float(np.corrcoef(rx,ry)[0,1])
def assoc(frame,x,y,adjust=False,boots=5000):
    cols=[x,y]+(['diagnosis_GBM','median_informative_molecules'] if adjust else [])
    q=frame[cols].dropna();n=len(q)
    base=dict(n=n,effect=np.nan,CI_low=np.nan,CI_high=np.nan,p=np.nan,
              loo_low=np.nan,loo_high=np.nan,loo_same_sign_fraction=np.nan,adjusted=adjust)
    if n<6 or q[x].nunique()<2 or q[y].nunique()<2:return base
    xx=q[x].to_numpy(float);yy=q[y].to_numpy(float)
    cov=q[['diagnosis_GBM','median_informative_molecules']].to_numpy(float) if adjust else None
    r=correlation(xx,yy,cov);rng=np.random.default_rng(SEED)
    draws=[];loo=[]
    for _ in range(boots):
        ix=rng.integers(0,n,n);draws.append(correlation(xx[ix],yy[ix],cov[ix] if adjust else None))
    for i in range(n):
        keep=np.arange(n)!=i;loo.append(correlation(xx[keep],yy[keep],cov[keep] if adjust else None))
    draws=np.asarray(draws);draws=draws[np.isfinite(draws)]
    if not adjust:p=float(spearmanr(xx,yy).pvalue)
    else:
        # Residual rank permutation; exploratory small-cohort adjustment.
        d=np.column_stack([np.ones(n)]+[rankdata(cov[:,j]) for j in range(cov.shape[1])])
        rx=rankdata(xx);ry=rankdata(yy);rx-=d@np.linalg.lstsq(d,rx,rcond=None)[0];ry-=d@np.linalg.lstsq(d,ry,rcond=None)[0]
        perm=np.array([rng.permutation(n) for _ in range(4999)])
        h=d@np.linalg.pinv(d);res=ry[perm];res-=res@h.T
        den=np.linalg.norm(res,axis=1)*np.linalg.norm(rx)
        null=np.divide(res@rx,den,out=np.full(len(res),np.nan),where=den>1e-12)
        p=(1+np.sum(np.abs(null)>=abs(r)))/5000
    base.update(effect=r,CI_low=float(np.quantile(draws,.025)),CI_high=float(np.quantile(draws,.975)),p=p,
                loo_low=float(np.nanmin(loo)),loo_high=float(np.nanmax(loo)),
                loo_same_sign_fraction=float(np.mean(np.sign(loo)==np.sign(r))))
    return base
def paired_difference(values):
    x=np.asarray(values,float);x=x[np.isfinite(x)];n=len(x)
    out=dict(n=n,effect=np.nan,CI_low=np.nan,CI_high=np.nan,p=np.nan,loo_low=np.nan,loo_high=np.nan,loo_same_sign_fraction=np.nan)
    if n<6:return out
    rng=np.random.default_rng(SEED);effect=x.mean();boot=x[rng.integers(0,n,(5000,n))].mean(1)
    signs=rng.choice([-1,1],(49999,n));null=(signs*x).mean(1)
    loo=np.array([np.delete(x,i).mean() for i in range(n)])
    out.update(effect=effect,CI_low=np.quantile(boot,.025),CI_high=np.quantile(boot,.975),
               p=(1+(np.abs(null)>=abs(effect)-1e-14).sum())/50000,
               loo_low=loo.min(),loo_high=loo.max(),loo_same_sign_fraction=np.mean(np.sign(loo)==np.sign(effect)))
    return out
def diagnosis_difference(frame,field):
    import itertools
    x=frame.loc[frame.diagnosis_GBM==1,field].dropna().to_numpy()
    y=frame.loc[frame.diagnosis_GBM==0,field].dropna().to_numpy();a=len(x);b=len(y)
    out=dict(n=a+b,n_GBM=a,n_MEN=b,effect=np.nan,CI_low=np.nan,CI_high=np.nan,p=np.nan,loo_low=np.nan,loo_high=np.nan,loo_same_sign_fraction=np.nan)
    if min(a,b)<3:return out
    effect=np.median(x)-np.median(y);allx=np.r_[x,y];n=len(allx);null=[]
    for comb in itertools.combinations(range(n),b):
        mask=np.zeros(n,bool);mask[list(comb)]=True
        null.append(np.median(allx[~mask])-np.median(allx[mask]))
    rng=np.random.default_rng(SEED)
    boot=np.median(x[rng.integers(0,a,(5000,a))],axis=1)-np.median(y[rng.integers(0,b,(5000,b))],axis=1)
    loo=np.r_[[np.median(np.delete(x,i))-np.median(y) for i in range(a)],
              [np.median(x)-np.median(np.delete(y,i)) for i in range(b)]]
    out.update(effect=effect,CI_low=np.quantile(boot,.025),CI_high=np.quantile(boot,.975),
               p=np.mean(np.abs(null)>=abs(effect)-1e-14),loo_low=loo.min(),loo_high=loo.max(),
               loo_same_sign_fraction=np.mean(np.sign(loo)==np.sign(effect)))
    return out
def sample_summary(sample,panels):
    x=metric_frame(sample);cl=clock();summaries=[];loci=[];cdf=[]
    for depth in [20,10]:
        z=x.loc[x.informative_molecules>=depth].copy()
        z['binary_rarefied']=z[f'binary_rarefied{depth}'];z['three_state_rarefied']=z[f'three_state_rarefied{depth}']
        z['hidden_rarefied_bits']=z[f'hidden_rarefied_bits{depth}']
        panel=panels[depth];common=z.loc[z[f'support{depth}']>=panel['required_patients']]
        row=dict(sample_id=sample,coverage_threshold=depth,valid_windows=len(z),common_panel_windows=len(common),
                 common_panel_required_patients=panel['required_patients'],
                 association_eligible=len(common)>=panel['minimum_patient_windows'])
        for scope,frame in [('genome',z),('common',common)]:
            row[f'{scope}_median_informative_molecules']=float(frame.informative_molecules.median())
            for metric in list(dict.fromkeys(FIELDS+TEST_FIELDS)):
                quantiles=frame[metric].quantile([.1,.25,.5,.75,.9])
                for q,value in quantiles.items():
                    row[f'{scope}_{metric}_q{int(100*q)}']=float(value)
        summaries.append(row)
        for metric in ['PDR','binary_entropy','three_state_entropy']:
            values=np.sort(z[metric].to_numpy())
            for value in np.linspace(0,1,101):
                cdf.append(dict(sample_id=sample,coverage_threshold=depth,metric=metric,value=value,cumulative_fraction=float(np.searchsorted(values,value,side='right')/len(values)) if len(values) else np.nan))
        # A fixed deterministic random control pool; selection does not depend on heterogeneity.
        eligible=z.loc[(~z.excluded_clock)&(z.density>0)]
        pool=eligible.sample(n=min(500000,len(eligible)),random_state=SEED).reset_index(drop=True)
        features=np.log(pool[['informative_molecules','span','density']].to_numpy(float))/np.log(1.25)
        trees={};groups={}
        for context,g in pool.groupby('context'):
            ids=g.index.to_numpy();groups[context]=ids;trees[context]=cKDTree(features[ids])
        for marker in FIVE:
            t=cl.loc[cl.CpGmarker==marker].iloc[0];ch=int(t.hg38_chromosome[3:]);pos=int(t.hg38_start)
            # Require all four sites to lie inside the +/-1kb target region.
            target=z.loc[(z.chromosome==ch)&(z.position0>=pos-1000)&(z.position0+z.span<=pos+1002)]
            pairs=[];used=set()
            for rr in target.itertuples():
                context=rr.context
                if context not in trees:continue
                query=np.log([rr.informative_molecules,rr.span,rr.density])/np.log(1.25)
                ids=trees[context].query_ball_point(query,r=1,p=np.inf)
                if len(ids)<5:continue
                ids=groups[context][ids]
                dist=((features[ids]-query)**2).sum(1)
                picked=ids[np.argsort(dist,kind='stable')[:10]]
                controls=pool.iloc[picked];used.update(controls['indices'].tolist())
                pair={k:getattr(rr,k)-float(controls[k].median()) for k in TEST_FIELDS}
                pairs.append(pair)
                assert (controls.context==context).all()
                for j,key in enumerate(['informative_molecules','span','density']):
                    assert np.max(np.abs(np.log(controls[key].to_numpy(float)/getattr(rr,key))))<=np.log(1.25)+1e-10
            valid=len(pairs)>=5 and len(pairs)>=.7*len(target)
            row=dict(sample_id=sample,coverage_threshold=depth,CpGmarker=marker,eligible_target_windows=len(target),
                     matched_target_windows=len(pairs),matched_fraction=len(pairs)/len(target) if len(target) else np.nan,
                     distinct_control_windows=len(used),match_eligible=valid)
            for k in TEST_FIELDS:
                row[k]=float(target[k].median()) if len(target)>=5 else np.nan
                row['matched_difference_'+k]=float(np.median([p[k] for p in pairs])) if valid else np.nan
            loci.append(row)
    del x
    return summaries,loci,cdf
def statistics(summary,loci):
    a=cohort();rows=[];locrows=[]
    for depth in [20,10]:
        z=summary.loc[(summary.coverage_threshold==depth)&summary.association_eligible].copy()
        for k in TEST_FIELDS:z[k]=z[f'common_{k}_q50']
        z['median_informative_molecules']=z.common_median_informative_molecules
        z=z.merge(a,on='sample_id',validate='one_to_one')
        for field in TEST_FIELDS:
            for y in ['primary_purity','DNAm_age_acceleration','relative_abs_log2_gt_0p2_fraction']:
                for adjust in [False,True]:
                    row=assoc(z,field,y,adjust);row.update(coverage_threshold=depth,metric=field,outcome=y,test='partial_Spearman' if adjust else 'Spearman')
                    rows.append(row)
            row=diagnosis_difference(z,field);row.update(coverage_threshold=depth,metric=field,outcome='GBM_minus_MEN',test='exact_label_permutation_median');rows.append(row)
        ll=loci.loc[loci.coverage_threshold==depth].merge(a,on='sample_id',validate='many_to_one')
        for marker,g in ll.groupby('CpGmarker'):
            for field in TEST_FIELDS:
                row=paired_difference(g['matched_difference_'+field]);row.update(coverage_threshold=depth,CpGmarker=marker,metric=field,outcome='target_minus_matched_controls',test='patient_sign_flip_mean');locrows.append(row)
            for field in ['PDR','binary_entropy','three_state_entropy']:
                for y in ['primary_purity','DNAm_age_acceleration']:
                    row=assoc(g,field,y);row.update(coverage_threshold=depth,CpGmarker=marker,metric=field,outcome=y,test='Spearman');locrows.append(row)
    stats=pd.DataFrame(rows);ls=pd.DataFrame(locrows)
    for x in [stats,ls]:
        x['BH_FDR']=np.nan
        for depth in [20,10]:
            mask=x.coverage_threshold==depth;x.loc[mask,'BH_FDR']=bh(x.loc[mask,'p'])
        x['robust_numerical_support']=(x.BH_FDR<.05)&(x.CI_low*x.CI_high>0)&(x.loo_same_sign_fraction==1)
    write(stats,O/'tables/patient_associations.tsv');write(ls,O/'tables/five_region_tests.tsv')
def run():
    annotate()
    samples=cohort().sample_id.tolist()
    import concurrent.futures
    with concurrent.futures.ProcessPoolExecutor(max_workers=8) as pool:
        for sample in pool.map(prepare_sample,samples):print(sample,'metric cache ready',flush=True)
    panels=common_panel();summaries=[];loci=[];cdf=[]
    for sample in samples:
        ss,ll,cc=sample_summary(sample,panels);summaries+=ss;loci+=ll;cdf+=cc
        print(sample,'summary and matched controls ready',flush=True)
    summary=pd.DataFrame(summaries);locus=pd.DataFrame(loci)
    write(summary,O/'tables/cohort_local_epialleles.tsv');write(locus,O/'tables/five_region_patient_comparisons.tsv')
    write(pd.DataFrame(cdf),O/'tables/local_window_CDF.tsv')
    statistics(summary,locus)
    save(O/'validation/analysis_complete.json',dict(complete=True,samples=len(samples),random_seed=SEED))
if __name__=='__main__':
    if sys.argv[1]=='annotate':annotate()
    elif sys.argv[1]=='prepare':prepare_sample(sys.argv[2])
    elif sys.argv[1]=='run':run()
    else:raise SystemExit('Use annotate, prepare SAMPLE, or run')
