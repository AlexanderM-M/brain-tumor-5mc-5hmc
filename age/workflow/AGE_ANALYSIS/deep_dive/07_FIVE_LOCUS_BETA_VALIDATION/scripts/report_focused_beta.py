"""Aggregate existing beta-matched patient caches with the retained statistics."""
import os
from pathlib import Path
from focused_beta import *
from scipy.special import xlogy
from composition_focused import EXTRA


def infer():
 samples=cohort().sample_id.tolist()
 allrec=pd.DataFrame([r for s in samples for r in json.loads((W/f'{s}_analysis.json').read_text())])
 comps=pd.DataFrame([r for s in samples for r in json.loads((W/f'{s}_composition.json').read_text())])
 allrec=allrec.merge(comps,on=['sample_id','locus','scenario'],how='left')
 nullcache={s:{**dict(np.load(W/f'{s}_nulls.npz')),**dict(np.load(W/f'{s}_composition_nulls.npz'))} for s in samples}
 rng=np.random.default_rng(SEED);draws=rng.integers(0,1000,(20,9999))
 rows=[];loo_records=[]
 for (scenario,loc),g0 in allrec.groupby(['scenario','locus'],sort=False):
  g=g0.loc[g0.patient_eligible].copy()
  for field in METRICS+EXTRA:
   if field in EXTRA and scenario!='primary':continue
   row=dict(locus=loc,scenario=scenario,metric=field,n_patients=len(g),patients=';'.join(g.sample_id),
    eligible_patient_candidates=int((g0.eligible_windows>=5).sum()),strict_matching_failed_patients=';'.join(g0.loc[(g0.eligible_windows>=5)&(~g0.patient_eligible),'sample_id']),
    total_eligible_windows=int(g0.eligible_windows.sum()),total_matched_windows=int(g0.matched_windows.sum()),
    mean_abs_beta_mismatch=g.mean_abs_beta_mismatch.mean(),beta_target=g.get('beta_target',pd.Series(dtype=float)).mean(),beta_control=g.get('beta_control',pd.Series(dtype=float)).mean(),
    effect=np.nan,CI_low=np.nan,CI_high=np.nan,empirical_p=np.nan,
    effect_normalized=np.nan,CI_low_normalized=np.nan,CI_high_normalized=np.nan,
    target=np.nan,control=np.nan,loo_effect_low=np.nan,loo_effect_high=np.nan,loo_same_sign_fraction=np.nan,loo_max_p=np.nan,synchronous_patient_p=np.nan)
   if len(g)>=6:
    eff=g[field].to_numpy();n=len(eff);effect=eff.mean()
    boot=eff[rng.integers(0,n,(5000,n))].mean(1)
    qs=np.stack([nullcache[s][f'{loc}|{scenario}|{field}'] for s in g.sample_id])
    null=np.stack([qs[j,draws[samples.index(s)]] for j,s in enumerate(g.sample_id)])
    mean_null=null.mean(0);p=(1+(np.abs(mean_null)>=abs(effect)-1e-12).sum())/10000
    le=[];lp=[]
    for j,s in enumerate(g.sample_id):
     e=np.delete(eff,j).mean();nul=np.delete(null,j,axis=0).mean(0);pv=(1+(np.abs(nul)>=abs(e)-1e-12).sum())/10000
     le.append(e);lp.append(pv)
    for s in samples:
     j=g.sample_id.tolist().index(s) if s in g.sample_id.tolist() else None
     loo_records.append(dict(scenario=scenario,locus=loc,metric=field,omitted_patient=s,p=p if j is None else lp[j]))
    factor=4 if field.startswith('binary') or field in ['M_bits','H_bits'] else LOG81 if field.startswith('three_state') else 1
    row.update(effect=effect,CI_low=np.quantile(boot,.025),CI_high=np.quantile(boot,.975),empirical_p=p,
      effect_normalized=effect/factor,CI_low_normalized=np.quantile(boot,.025)/factor,CI_high_normalized=np.quantile(boot,.975)/factor,
      target=g[field+'_target'].mean(),control=g[field+'_control'].mean(),target_normalized=g[field+'_target'].mean()/factor,control_normalized=g[field+'_control'].mean()/factor,normalization_denominator=factor,
      loo_effect_low=min(le),loo_effect_high=max(le),loo_same_sign_fraction=np.mean(np.sign(le)==np.sign(effect)),loo_max_p=max(lp),
      synchronous_patient_p=(1+(np.abs(qs.mean(0))>=abs(effect)).sum())/1001)
   rows.append(row)
 stats=pd.DataFrame(rows);loo=pd.DataFrame(loo_records)
 # Always a family of five hypotheses, with unavailable hypotheses assigned p=1.
 for _,ix in stats.groupby(['scenario','metric']).groups.items():
  stats.loc[ix,'BH_FDR_five_loci']=bh(stats.loc[ix,'empirical_p'].fillna(1))
 for scenario in stats.scenario.unique():
  ix=stats.index[(stats.scenario==scenario)&stats.metric.isin(['PDR','binary_bits'])]
  stats.loc[ix,'BH_FDR_joint10_PDR_binary']=bh(stats.loc[ix,'empirical_p'].fillna(1))
 for _,ix in loo.groupby(['scenario','metric','omitted_patient']).groups.items():
  p=loo.loc[ix,'p'].to_numpy();order=np.argsort(p);v=np.minimum.accumulate((p[order]*5/np.arange(1,len(p)+1))[::-1])[::-1]
  q=np.empty(len(p));q[order]=np.minimum(1,v);loo.loc[ix,'q']=q
 for idx,row in stats.iterrows():
  l=loo.loc[(loo.scenario==row.scenario)&(loo.locus==row.locus)&(loo.metric==row.metric)]
  stats.loc[idx,'LOO_max_FDR_five_loci']=l.q.max() if len(l) else np.nan
  stats.loc[idx,'LOO_FDR_pass_fraction']=float((l.q<.05).mean()) if len(l) else np.nan
 stats.to_csv(F/'tables/five_locus_statistics.tsv',sep='\t',index=False,na_rep='NA',float_format='%.9g')
 frequencies=pd.concat([pd.read_csv(W/f'{s}_frequencies.tsv',sep='\t',dtype={'pattern':str}) for s in samples if (W/f'{s}_frequencies.tsv').stat().st_size>2],ignore_index=True)
 agg=frequencies.groupby(['locus','arm','alphabet','pattern'],sort=False).frequency.agg(['mean','count']).reset_index().rename(columns={'mean':'frequency','count':'n_patients'})
 agg['sample_id']='cohort_equal_patient_mean'
 frequencies['n_patients']=1
 frequencies=pd.concat([frequencies,agg],ignore_index=True)
 frequencies.to_csv(F/'tables/epiallele_frequencies.tsv',sep='\t',index=False,na_rep='NA',float_format='%.9g')
 save(F/'validation/patient_matching_and_effects.json',allrec.to_dict('records'))
 save(F/'validation/leave_one_patient_out.json',loo.to_dict('records'))
 print(stats.loc[(stats.scenario=='primary')&stats.metric.isin(['PDR','binary_bits','hidden_bits','M_bits','H_bits','PDR_beta_residual','binary_bits_beta_residual']),['locus','metric','n_patients','effect','CI_low','CI_high','BH_FDR_five_loci','LOO_max_FDR_five_loci']].to_string(index=False),flush=True)

if __name__=='__main__': infer()
