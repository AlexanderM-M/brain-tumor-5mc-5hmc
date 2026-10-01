"""Retained WGS scientific qualification, purity association and five-CpG analysis.
Starts from de-identified retained segmentation/QC evidence; no sequence decoding.
"""
import os
from pathlib import Path
import sys,json,hashlib,importlib.util,collections
import numpy as np,pandas as pd
R=Path(os.environ['AGE_WORKSPACE']);D=R/'AGE_ANALYSIS/deep_dive';N=R/'AGE_ANALYSIS/next_generation';O=D/'FINAL_INTEGRATION'
W=Path(os.environ['AGE_WGS_PROCESSED']);SEED=20260921


def write(df,path):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    df.to_csv(path,sep='\t',index=False,na_rep='NA',float_format='%.10g')

def save(obj,path):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,indent=2,default=lambda x:x.item() if hasattr(x,'item') else str(x))+'\n')

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m

def fit_stability():
    """Chromosome-block bootstrap of retained segments, conditional on allele-derived search range.
    Reproduces original SAVANA objective/local-minimum/viability/ranking rules; no WGS recounting.
    """
    from scipy.ndimage import minimum_filter
    q=pd.read_csv(D/'04_WGS_genomic_backbone/tables/cohort_WGS_QC.tsv',sep='\t')
    segments=pd.read_csv(O/'provenance/WGS_segment_evidence.tsv',sep='\t')
    for i,row in q.loc[q.savana_fitted_purity.notna()].iterrows():
        s=row.sample_id;cache=O/'provenance'/f'{s}_final_fit_stability.json'
        if cache.exists():
            for k,v in json.loads(cache.read_text()).items():q.loc[i,k]=v
            continue
        g=segments.loc[segments.sample_id==s].copy()
        # Original author table follows chr1,chr2,...; preserve exact segment order from original table.
        ab=pd.read_csv(W/s/'savana_cna'/f'{s}_segmented_absolute_copy_number.tsv',sep='\t')
        g=g.set_index('segment_id').loc[ab.segment_id].reset_index()
        rel=2**g.log2ratio.to_numpy();weights=g.covered_bp.to_numpy()/10000
        ch=g.chromosome.str[3:].astype(int).to_numpy()-1
        p0=float(str(row.preliminary_cellularity).split(',')[-1])
        pseq=np.round(np.arange(.01,1.005,.01),2)
        plseq=np.round(np.arange(1.5,5.005,.01),2)
        purs,plos=np.meshgrid(pseq,plseq,indexing='ij');pf=purs.ravel();lf=plos.ravel()
        acn=lf[:,None]+(rel[None,:]-1)*(lf[:,None]+2/pf[:,None]-2)
        ints=np.rint(acn).astype(np.int16);diff=np.abs(acn-ints)
        E=np.stack([np.sum((diff[:,ch==c]**2)*weights[ch==c],axis=1) for c in range(22)],axis=1)
        CW=np.array([weights[ch==c].sum() for c in range(22)])
        def pick(mult,prior):
            allowed=(pseq>=round(max(0,prior-.1),2)-1e-8)&(pseq<=round(min(1,prior+.1),2)+1e-8)
            d=np.sqrt((E@mult)/(CW@mult)).reshape(len(pseq),len(plseq));d[~allowed]=np.inf
            loc=(d<=minimum_filter(d,size=3,mode='constant',cval=np.inf))&np.isfinite(d)
            if not loc.any():return np.nan,np.nan
            best=d[loc].min();ix=np.flatnonzero((loc&(d<best*1.25)).ravel())
            ix=sorted(ix,key=lambda k:(round(float(d.ravel()[k]),3),lf[k]))
            ws=weights*mult[ch];den=ws.sum()
            for k in ix:
                if ws[ints[k]<=0].sum()/den>.1 or ws[diff[k]<.25].sum()/den<.5:continue
                counts=collections.Counter()
                for val,mc in zip(ints[k],mult[ch]):
                    if mc:counts[int(val)]+=int(mc)
                common=counts.most_common(2)
                if len(common)<2 or abs(common[0][0]-common[1][0])>1:continue
                return float(pf[k]),float(lf[k])
            return np.nan,np.nan
        base=pick(np.ones(22,dtype=int),p0)
        if not np.allclose(base,[row.savana_fitted_purity,row.savana_fitted_ploidy],atol=.011):
            save(dict(sample_id=s,baseline=base,retained=[row.savana_fitted_purity,row.savana_fitted_ploidy],status='baseline_reproduction_failed'),cache)
            print(s,'BASELINE MISMATCH',base,flush=True);continue
        rng=np.random.default_rng(SEED+100+i);boot=np.array([pick(rng.multinomial(22,np.ones(22)/22),p0) for _ in range(2000)])
        loo=np.array([pick(1-np.eye(22,dtype=int)[c],p0) for c in range(22)])
        valid=np.isfinite(boot).all(axis=1);ci=np.nanquantile(boot[:,0],[.025,.975])
        prior_lo=float(row.preliminary_block_bootstrap_CI_low);prior_hi=float(row.preliminary_block_bootstrap_CI_high)
        prior_sens=np.array([pick(np.ones(22,dtype=int),float(p)) for p in np.linspace(prior_lo,prior_hi,11)])
        max_loo=float(np.nanmax(abs(loo[:,0]-row.savana_fitted_purity)))
        prior_change=float(np.nanmax(abs(prior_sens[:,0]-row.savana_fitted_purity)))
        passed=bool(valid.mean()>=.95 and np.isfinite(loo).all() and np.isfinite(prior_sens).all() and ci[1]-ci[0]<=.20+1e-10 and max_loo<=.10+1e-10 and prior_change<=.10+1e-10 and row.near_best_purity_spread<=.1 and row.savana_fitted_ploidy not in [1.5,5] and row.preliminary_stability_screen_pass)
        d=dict(sample_id=s,final_fit_baseline_reproduced=True,final_fit_bootstrap_status='2000_chromosome_blocks_conditional_allele_prior_with_prior_CI_sensitivity',final_purity_bootstrap_low=ci[0],final_purity_bootstrap_high=ci[1],final_purity_bootstrap_width=ci[1]-ci[0],final_fit_bootstrap_valid_fraction=valid.mean(),final_purity_max_chromosome_omission_change=max_loo,final_purity_max_prior_CI_change=prior_change,final_fit_stability_screen_pass=passed,primary_purity=row.savana_fitted_purity if passed else None,primary_ploidy=row.savana_fitted_ploidy if passed else None,primary_purity_status='qualified_conditional_WGS_fit' if passed else 'unresolved_fit_instability',scientific_fit_class='accepted_QC_qualified_fit' if passed else 'valid_CNV_unresolved_purity_ploidy')
        save(d,cache)
        write(pd.DataFrame(boot,columns=['purity','ploidy']),O/'provenance'/f'{s}_fit_bootstrap.tsv')
        write(pd.DataFrame(loo,columns=['purity','ploidy']).assign(omitted_chromosome=[f'chr{c}' for c in range(1,23)]),O/'provenance'/f'{s}_fit_omission.tsv')
        for k,v in d.items():q.loc[i,k]=v
        write(q,D/'04_WGS_genomic_backbone/tables/cohort_WGS_QC.tsv')
        print(s,'final CI',np.round(ci,3),'valid',round(valid.mean(),3),'omission',round(max_loo,3),'qualified',passed,flush=True)
    write(q,D/'04_WGS_genomic_backbone/tables/cohort_WGS_QC.tsv')

def association(df,x,y,adjust=(),seed=SEED):
    from scipy.stats import rankdata
    a=df[[x,y]+list(adjust)].dropna().to_numpy(float);n=len(a)
    if n<6:return dict(n=n,rho=None,CI_low=None,CI_high=None,loo_min=None,loo_max=None,status='insufficient_patients',adjustment=','.join(adjust) or 'none')
    def stat(a):
        v=rankdata(a,axis=0);vx=v[:,0];vy=v[:,1]
        Z=np.column_stack([np.ones(len(a)),v[:,2:]])
        vx=vx-Z@np.linalg.lstsq(Z,vx,rcond=None)[0];vy=vy-Z@np.linalg.lstsq(Z,vy,rcond=None)[0]
        den=np.linalg.norm(vx)*np.linalg.norm(vy)
        return np.dot(vx,vy)/den if den>1e-10 else np.nan
    r=stat(a);rng=np.random.default_rng(seed)
    b=np.array([stat(a[rng.integers(0,n,n)]) for _ in range(5000)])
    ci=np.nanquantile(b,[.025,.975]);loo=np.array([stat(np.delete(a,i,axis=0)) for i in range(n)])
    return dict(n=n,rho=r,CI_low=ci[0],CI_high=ci[1],loo_min=np.nanmin(loo),loo_max=np.nanmax(loo),bootstrap_valid=int(np.isfinite(b).sum()),status='evaluated',adjustment=','.join(adjust) or 'none')

def integrate():
    q=pd.read_csv(D/'04_WGS_genomic_backbone/tables/cohort_WGS_QC.tsv',sep='\t')
    a=pd.read_csv(N/'final/tables/age_patient_summary.tsv',sep='\t')
    a=a.merge(q[['sample_id','savana_fitted_purity','savana_fitted_ploidy','primary_purity','primary_ploidy','primary_purity_status','relative_abs_log2_gt_0p2_fraction']],on='sample_id',validate='one_to_one')
    a['diagnosis_GBM']=(a.tumour_type=='Glioblastoma').astype(int)
    bounds=q.set_index('sample_id').preliminary_cellularity.astype(str).str.split(',').str[-1].astype(float)
    a['purity_at_search_boundary']=[bool(np.isfinite(p) and min(abs(p-round(max(0,b-.1),2)),abs(p-round(min(1,b+.1),2)))<.005) for p,b in zip(a.savana_fitted_purity,bounds.loc[a.sample_id])]
    q['purity_at_search_boundary']=q.sample_id.map(a.set_index('sample_id').purity_at_search_boundary)
    write(q,D/'04_WGS_genomic_backbone/tables/cohort_WGS_QC.tsv')
    stats=[]
    for name,subset,col in [('primary_qualified_WGS',a,'primary_purity'),('sensitivity_all_17_algorithm_fits',a,'savana_fitted_purity'),('sensitivity_GBM_only',a.loc[a.diagnosis_GBM==1],'primary_purity'),('sensitivity_exclude_search_boundary',a.loc[~a.purity_at_search_boundary],'primary_purity')]:
        for adj in [(),('chronological_age',),('clock_CpG_coverage_fraction',),('chronological_age','clock_CpG_coverage_fraction'),('diagnosis_GBM',)]:
            if name.endswith('GBM_only') and 'diagnosis_GBM' in adj:continue
            d=association(subset,col,'DNAm_age_acceleration',adj);d.update(analysis=name,x=col,y='DNAm_age_acceleration');stats.append(d)
    for x,y in [('primary_purity','tumour_content'),('savana_fitted_purity','tumour_content'),('tumour_content','DNAm_age_acceleration')]:
        d=association(a,x,y);d.update(analysis='purity_comparison',x=x,y=y);stats.append(d)
    write(pd.DataFrame(stats),O/'tables/purity_age_statistics.tsv')
    a['ONT_brain_age']=np.nan;a['ONT_brain_DeltaAge']=np.nan;a['ONT_feature_coverage']=np.nan
    a['ONT_status']='blocked_exact_pretrained_model_features_scaling_unavailable'
    write(a,O/'tables/paired_patient_clock_WGS.tsv')
    print(pd.DataFrame(stats).to_string(index=False),flush=True)

def five():
    # CPU affinity is managed by the calling environment.
    sys.path.insert(0,str(D/'06_independent_preparation/scripts'))
    f=module('five_null_post',D/'06_independent_preparation/scripts/five_null.py')
    c=pd.read_csv(R/'data/reference/horvath2013/clock353_hg38.tsv',sep='\t')
    a=pd.read_csv(O/'tables/paired_patient_clock_WGS.tsv',sep='\t')
    ref=pd.read_csv(R/'data/reference/horvath2013/clock353_sesame_reference.csv').set_index('Probe_ID').loc[c.CpGmarker,'median'].to_numpy()
    counts=np.stack([np.load(R/f'results/age_hox_joint/tables/{s}_clock353_counts.npz')['counts'] for s in a.sample_id])
    dp=counts.sum(axis=2);beta=np.divide(counts[:,:,1:].sum(axis=2),dp,out=np.zeros(dp.shape,float),where=dp>0);beta=np.where(dp>=5,beta,ref)
    w=c.CoefficientTraining.to_numpy();ids=list(c.CpGmarker);selected=['cg07158339','cg22947000','cg01511567','cg19761273','cg13460409'];ix=[ids.index(x) for x in selected]
    effect=f.replacement_effect(beta,ref,w,ix)
    assert np.allclose(effect['full_age'],a.Age_BS,atol=1e-7)
    a['five_raw_signed_score']=(beta[:,ix]*w[ix]).sum(axis=1)
    a['remaining_raw_signed_score']=(beta*w).sum(axis=1)-a.five_raw_signed_score
    a['five_reference_contrast_score']=effect['score_contribution']
    a['five_reference_contrast_years']=effect['contribution_years']
    a['remaining_DeltaAge']=effect['counterfactual_age']-a.chronological_age
    a['five_signed_fraction_of_offset']=np.divide(a.five_reference_contrast_years,a.DNAm_age_acceleration,out=np.full(len(a),np.nan),where=abs(a.DNAm_age_acceleration)>1e-8)
    previous=pd.read_csv(D/'05_five_locus_signal/tables/five_locus_patient_contributions.tsv',sep='\t').set_index('sample_id').loc[a.sample_id]
    assert np.allclose(a.five_reference_contrast_years,previous.five_locus_distortion_contribution,atol=1e-7)
    ann=pd.read_csv(N/'tables/clock_locus_chemical_modes.tsv',sep='\t').set_index('CpGmarker').loc[ids]
    feats=[dict(id=ids[j],weight=w[j],median_depth=float(np.median(dp[:,j])),eligible=bool(np.all(dp[:,j]>=5) and w[j]!=0),island_context=ann.iloc[j].CpG_island_context,genic_context=ann.iloc[j].genic_context) for j in range(len(ids))]
    pools={k:f.matching_pools(feats,selected,k) for k in ['strict','island_only','calipers_only']}
    sets=f.draw_matched_sets(pools['calipers_only'])
    null=f.conditional_null(beta,ref,w,ids,selected,sets)
    write(pd.DataFrame(dict(replicate=np.arange(1,len(sets)+1),CpGs=[','.join(s) for s in sets],mean_replacement_effect_years=null['null'])),O/'tables/five_CpG_caliper_null.tsv')
    rows=[]
    for context,p in pools.items():
        for target,candidates in p.items():rows.append(dict(matching=context,target=target,n_candidates=len(candidates),candidate_CpGs=','.join(candidates)))
    write(pd.DataFrame(rows),O/'provenance/five_matching_pools.tsv')
    stats=[]
    for y in ['DNAm_age_acceleration','remaining_DeltaAge','primary_purity','savana_fitted_purity','ONT_brain_DeltaAge']:
        d=association(a,'five_reference_contrast_years',y);d.update(analysis='five_locus_contribution',x='five_reference_contrast_years',y=y,part_whole_coupling=y=='DNAm_age_acceleration');stats.append(d)
    write(pd.DataFrame(stats),O/'tables/five_CpG_associations.tsv')
    write(a,O/'tables/paired_patient_clock_WGS.tsv')
    save(dict(mean_full_offset=a.DNAm_age_acceleration.mean(),mean_five_effect=a.five_reference_contrast_years.mean(),remaining_mean_offset=a.remaining_DeltaAge.mean(),ratio_cohort_means=a.five_reference_contrast_years.mean()/a.DNAm_age_acceleration.mean(),strict_matching_feasible=all(pools['strict'].values()),pool_sizes={k:{x:len(v) for x,v in p.items()} for k,p in pools.items()},null_seed=f.SEED,null_draws=len(sets),null_unique_sets=len(set(tuple(sorted(s)) for s in sets)),null_quantiles=np.quantile(null['null'],[.025,.5,.975]).tolist(),conditional_tail=null['conditional_tail'],selection_adjusted=False,context_matched=False),O/'provenance/five_results.json')
    print((O/'provenance/five_results.json').read_text(),flush=True)

def ploidy_review():
    q=pd.read_csv(D/'04_WGS_genomic_backbone/tables/cohort_WGS_QC.tsv',sep='\t')
    for i,row in q.loc[q.savana_fitted_purity.notna()].iterrows():
        b=pd.read_csv(O/'provenance'/f'{row.sample_id}_fit_bootstrap.tsv',sep='\t')
        lo,hi=b.ploidy.quantile([.025,.975])
        q.loc[i,'final_ploidy_bootstrap_low']=lo;q.loc[i,'final_ploidy_bootstrap_high']=hi
        q.loc[i,'alternative_ploidy_bootstrap_fraction']=float((abs(b.ploidy-row.savana_fitted_ploidy)>.5).mean())
        # Distinct copy-number baselines in the bootstrap are evidence against a stable joint fit.
        if hi-lo>.5:
            q.loc[i,'primary_purity']=np.nan;q.loc[i,'primary_ploidy']=np.nan
            q.loc[i,'primary_purity_status']='unresolved_joint_fit_ploidy_ambiguity'
            q.loc[i,'scientific_fit_class']='valid_CNV_unresolved_purity_ploidy'
            q.loc[i,'final_fit_stability_screen_pass']=False
    q['outcome_class']=np.where(q.technical_status!='PASS','genuine_technical_failure',np.where(q.primary_purity.notna(),'accepted_QC_qualified_purity_ploidy_fit','valid_CNV_unresolved_purity_ploidy'))
    base=pd.read_csv(N/'tables/age_patient_summary.tsv',sep='\t').set_index('sample_id')
    q['original_aligned_depth']=q.sample_id.map(base.overall_aligned_depth)
    q['alternative_ASCAT_status']='below_documented_tumour_only_WGS_hg38_50X_preset;no_matched_normal;unvalidated_custom_genotype_model'
    write(q,D/'04_WGS_genomic_backbone/tables/cohort_WGS_QC.tsv')
    save(dict(reason='Joint-fit QC must consider distinct bootstrap ploidy modes, not purity width alone.',ploidy_ambiguity_rule='95% ploidy interval spans >0.5 copies; operational audit rule added after fitting QC, not optimized against age association',consequence='GBM-13 excluded from primary joint-fit-qualified analysis; previously reported 14-case association remains sensitivity evidence'),O/'provenance/ploidy_review_amendment.json')
    integrate();five()

if __name__=='__main__':
    allowed={'fit_stability':fit_stability,'integrate':integrate,'five':five,'ploidy_review':ploidy_review}
    allowed[sys.argv[1]]()
