#!/usr/bin/env python3
"""Run only the fixed five-locus external validation from the verified cache."""
import os
from pathlib import Path
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
from curate_inputs import load
from read_archive import ROOT, PROJECT, LOCI
import targeted_statistics as st


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--cache',type=Path,default=Path(os.environ['AGE_EXTERNAL_CACHE']))
    args=ap.parse_args()
    nm,nv,cohorts,patients,ont,audit=load(args.cache)
    seed=20260930
    rng=np.random.default_rng(seed)
    normal_rows,gbm_rows,normal_fits=[],[],{}
    coeff=pd.read_csv(PROJECT/'data/reference/horvath2013/clock353_hg38.tsv',sep='\t').set_index('CpGmarker').CoefficientTraining
    refs=pd.read_csv(PROJECT/'data/reference/horvath2013/clock353_sesame_reference.csv').set_index('Probe_ID')['median']
    assert len(coeff)==353 and coeff.notna().all()
    clocks={}
    for name,c in cohorts.items():
        beta=c['beta'].reindex(columns=coeff.index)
        complete=beta.notna().mean(axis=1).ge(.95) & beta[LOCI].notna().all(axis=1) & c['meta'].age.notna()
        complete &= c['meta'].age.ge(18)
        beta=beta.loc[complete]
        # The supplied GEO cohorts are complete; no target/reference imputation is used.
        if len(beta):
            beta=beta.fillna(refs)
            assert beta.notna().all().all()
            eta=.696+beta@coeff
            eta5=.696+beta.assign(**{p:refs[p] for p in LOCI})@coeff
            ages=c['meta'].loc[complete,'age']
            full=st.inv_age(eta)-ages
            rest=st.inv_age(eta5)-ages
            contrib={p:st.inv_age(eta)-st.inv_age(eta-coeff[p]*(beta[p]-refs[p])) for p in LOCI}
            clocks[name]=(complete[complete].index,full,rest,contrib)
        audit.setdefault(name,{})['clock_qualified_n']=int(complete.sum())
        audit[name]['observed_clock_fraction_min']=float(c['beta'].notna().mean(axis=1).min())
        audit[name]['observed_clock_fraction_max']=float(c['beta'].notna().mean(axis=1).max())
    for p in LOCI:
        good=nv[p].notna()
        x=nm.loc[good,'age'].to_numpy(float);y=nv.loc[good,p].to_numpy(float)
        sp=st.spearman(x,y,rng)
        slope,intercept,boot=st.theil_bootstrap(x,y,rng)
        normal_fits[p]=(slope,intercept,boot,x,y)
        row=dict(cohort='GSE74193',CpGmarker=p,status='ANALYSED',n_analysed=len(x),
                 age_min_years=x.min(),age_max_years=x.max(),spearman_rho=sp['rho'],
                 rho_CI_low=sp['rho_CI_low'],rho_CI_high=sp['rho_CI_high'],age_association_p=sp['p'],
                 robust_beta_slope_per_decade=slope*10,slope_CI_low=st.ci(boot[:,0]*10)[0],
                 slope_CI_high=st.ci(boot[:,0]*10)[1],normal_beta_median=float(np.median(y)))
        nuisance=pd.concat([nm.loc[good,['male','neuron']],pd.get_dummies(nm.loc[good,'plate'],drop_first=True,dtype=float)],axis=1)
        adjusted=st.adjusted_rank_correlation(x,y,nuisance.to_numpy(float),rng)
        for k,v in adjusted.items():row['sex_neuron_plate_adjusted_'+k]=v
        for tumour,label in [('Glioblastoma','ONT_GBM'),('Meningioma','ONT_MEN')]:
            ids=patients.index[patients.tumour_type.eq(tumour)]
            d=st.trajectory_deviation(patients.loc[ids,'chronological_age'],ont.loc[ids,p],x,slope,intercept,boot,rng)
            for k,v in d.items():row[label+'_'+k]=v
        row['notes']='Detection P<=0.01; bestQC; dropsample FALSE; adult controls; one array per BrNum. Adjusted rank-correlation P uses approximate t inference.'
        normal_rows.append(row)
        print('NORMAL',p,sp['rho'],sp['p'],'slope/decade',slope*10,flush=True)
        for name,c in cohorts.items():
            meta,beta=c['meta'],c['beta']
            if c['control_beta'] is not None:
                tids=meta.index.difference(c.get('exclude_from_unpaired',[]),sort=False)
                t=beta.loc[tids,p]; control=c['control_beta'][p]
                contrast='within_study_non_tumour_control'
            else:
                tids=meta.index[(meta.age >= min(x)) & (meta.age <= max(x))]
                t=beta.loc[tids,p]
                age_min,age_max=meta.loc[tids,'age'].min(),meta.loc[tids,'age'].max()
                control=nv.loc[good & nm.age.between(age_min,age_max),p]
                contrast='cross_study_GSE74193_controls_in_tumour_age_range'
            d=st.group_difference(t,control,rng)
            r=dict(cohort=name,CpGmarker=p,comparison=contrast,status='ANALYSED' if d['n_GBM'] else 'UNAVAILABLE_PROBE',
                   eligible_primary_GBM=len(meta),**d,notes=c['notes'])
            age_d=st.trajectory_deviation(meta.age,beta[p],x,slope,intercept,boot,rng)
            for k,v in age_d.items():r['normal_trajectory_'+k]=v
            r['normal_trajectory_reference_n']=len(x)
            if name=='GSE60274':
                cm=c['control_meta']; cb=c['control_beta'][p]
                lo=max(meta.age.min(),cm.age.min());hi=min(meta.age.max(),cm.age.max())
                tok=meta.age.between(lo,hi);cok=cm.age.between(lo,hi)
                rob=st.robust_age_group(meta.loc[tok,'age'],beta.loc[tok,p],cm.loc[cok,'age'],cb.loc[cok],rng)
                for k,v in rob.items():r['within_study_age_adjusted_'+k]=v
            if name in clocks:
                ids,full,rest,contr=clocks[name]
                for outcome,values in [('full_acceleration',full),('leave_five_out_acceleration',rest)]:
                    cor=st.spearman(contr[p],values,rng)
                    for k,v in cor.items():r['contribution_vs_'+outcome+'_'+k]=v
                r['clock_status']='FIXED_353_COMPLETE_NO_IMPUTATION'
            else:
                r['clock_status']='NOT_ESTIMATED: <95% clock coverage and/or missing target; no imputed target validation'
            gbm_rows.append(r)
            print('GBM',name,p,'delta',d['delta_beta'],'age deviation',age_d['delta'],flush=True)
    normal=pd.DataFrame(normal_rows)
    normal['age_association_BH_q']=st.bh(normal.age_association_p)
    normal['sex_neuron_plate_adjusted_BH_q']=st.bh(normal.sex_neuron_plate_adjusted_p)
    gbm=pd.DataFrame(gbm_rows)
    for name,idx in gbm.groupby('cohort',sort=False).groups.items():
        assert len(idx)==5
        gbm.loc[idx,'BH_q']=st.bh(gbm.loc[idx,'p'])
        for out in ['full_acceleration','leave_five_out_acceleration']:
            col='contribution_vs_'+out+'_p'
            gbm.loc[idx,col.replace('_p','_BH_q')]=st.bh(gbm.loc[idx,col])
    summary=[]
    for p in LOCI:
        n=normal.set_index('CpGmarker').loc[p]
        g=gbm[gbm.CpGmarker.eq(p)].set_index('cohort')
        significant=[];directions=[]
        for cohort in ['GSE60274','GSE195684','TCGA_GBM']:
            r=g.loc[cohort]
            if np.isfinite(r.delta_beta) and r.BH_q<.05 and r.delta_CI_low*r.delta_CI_high>0:
                significant.append(cohort);directions.append(np.sign(r.delta_beta))
        age_assoc=n.age_association_BH_q<.05 and n.slope_CI_low*n.slope_CI_high>0
        replicated=len(significant)>=2 and len(set(directions))==1
        cls='inconclusive'
        if age_assoc and replicated:
            cls='concordant with normal ageing' if np.sign(n.robust_beta_slope_per_decade)==directions[0] else 'opposite to normal ageing'
        row=dict(CpGmarker=p,classification=cls,normal_age_associated=bool(age_assoc),
                 normal_age_direction='increasing' if age_assoc and n.robust_beta_slope_per_decade>0 else 'decreasing' if age_assoc else 'not established',
                 normal_age_BH_q=n.age_association_BH_q,normal_slope_beta_per_decade=n.robust_beta_slope_per_decade,
                 sex_neuron_plate_adjusted_BH_q=n.sex_neuron_plate_adjusted_BH_q,
                 external_GBM_replicated=bool(replicated),supporting_cohorts=';'.join(significant),
                 tumour_shift_direction='hypermethylated' if replicated and directions[0]>0 else 'hypomethylated' if replicated else 'inconclusive',
                 ONT_GBM_age_adjusted_delta_beta=n.ONT_GBM_delta,
                 ONT_GBM_CI_low=n.ONT_GBM_low,ONT_GBM_CI_high=n.ONT_GBM_high,
                 notes='Classification describes direction only; concordance does not imply normal ageing explains the magnitude. Cross-study comparisons and small control groups limit specificity; non-significance does not establish tumour specificity.')
        for cohort in ['GSE60274','GSE195684','TCGA_GBM']:
            r=g.loc[cohort]
            row[cohort+'_age_adjusted_delta_beta']=r.normal_trajectory_delta
            row[cohort+'_age_adjusted_CI_low']=r.normal_trajectory_low
            row[cohort+'_age_adjusted_CI_high']=r.normal_trajectory_high
        summary.append(row)
    summary=pd.DataFrame(summary)
    # Stage output in temporary cache; publishing occurs after numerical/figure QC.
    normal.to_csv(args.cache/'normal_brain_validation.tsv',sep='\t',index=False)
    gbm.to_csv(args.cache/'GBM_validation.tsv',sep='\t',index=False)
    summary.to_csv(args.cache/'five_locus_summary.tsv',sep='\t',index=False)
    np.savez_compressed(args.cache/'normal_bootstrap.npz',**{p:np.column_stack([v[2][:,0],v[2][:,1]]) for p,v in normal_fits.items()})
    audit['analysis_settings']={'seed':seed,'bootstrap_replicates':st.B,'permutations':st.P,
        'detection_p_max':.01,'age_min':18,'clock_coverage_min':.95,'numpy':np.__version__,
        'pandas':pd.__version__,'scipy':scipy.__version__}
    (args.cache/'analysis_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print('COMPLETE',flush=True)


if __name__=='__main__':main()
