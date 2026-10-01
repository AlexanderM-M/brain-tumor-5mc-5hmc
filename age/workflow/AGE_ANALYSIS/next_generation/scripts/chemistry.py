import os
from pathlib import Path
from common import *
from scipy.stats import rankdata
from scipy.cluster.hierarchy import linkage,fcluster
from scipy.spatial.distance import pdist

def mode(m,h,total,w,common,scale=1.):
 out=[]
 for j in range(len(w)):
  a,b,t=m[:,j],h[:,j],total[:,j];sm,sh,st=np.std(a),np.std(b),np.std(t);hm=np.mean(abs(21*w[j]*b));frac=np.sum(abs(b))/max(np.sum(abs(t)),1e-12)
  if not common[j] or max(sm,sh,st)<.02*scale:k='LOW-INFORMATION'
  elif min(sm,sh)>=.03*scale and float(corr(a,b))<=-.7 and st<=.5*scale*max(sm,sh):k='RECIPROCAL-COMPOSITION'
  elif hm>=.5*scale and frac>=.2*scale:k='5hmC-SENSITIVE'
  elif frac<=.1*scale and sm>=.02*scale:k='5mC-DOMINANT'
  elif st>=.05*scale and abs(float(corr(a,t)))>=.8:k='TOTAL-MODIFICATION DRIVEN'
  else:k='UNRESOLVED'
  out.append(k)
 return np.array(out)

def run():
 p,c,ref,counts=inputs();w=c.CoefficientTraining.to_numpy();depth=counts.sum(2);observed=depth>=5;v=np.divide(counts,depth[:,:,None],out=np.full(counts.shape,np.nan,float),where=depth[:,:,None]>0);u,m,h=v.transpose(2,0,1);total=m+h;assert np.nanmax(abs(u+m+h-1))<1e-12
 bs=np.where(observed,total,ref);mm=np.where(observed,m,ref);eta=.696+bs@w;age=inv(eta);age_m=inv(.696+mm@w);assert max(abs(age-p.predicted_DNAm_age))<1e-6;assert max(abs(age_m-p.predicted_age_5mC_only))<1e-6
 delta=age-p.chronological_age.to_numpy();contrast=21*w*(bs-ref);mc=np.where(observed,21*w*m,0);hc=np.where(observed,21*w*h,0);imp=np.where(observed,0,21*w*ref);assert (eta>=0).all() and ((.696+mm@w)>=0).all();assert np.max(abs(age-(34.616+mc.sum(1)+hc.sum(1)+imp.sum(1))))<1e-10
 pp=p.copy();pp['Age_BS']=age;pp['Age_M']=age_m;pp['hydroxymethylation_clock_effect']=age-age_m;pp['absolute_measured_5mC_signal_fraction']=np.sum(abs(mc),axis=1)/np.sum(abs(mc)+abs(hc),axis=1);table(pp,'tables/patient_base.tsv')
 rows=[]
 for i,label in enumerate(p.sample_id):
  for j,probe in enumerate(c.CpGmarker):rows.append(dict(sample_id=label,CpGmarker=probe,clock_weight=w[j],valid_depth=int(depth[i,j]),observed=observed[i,j],U=u[i,j],M=m[i,j],H=h[i,j],total_modification=total[i,j],modification_balance=m[i,j]/total[i,j] if total[i,j]>=.1 else np.nan,clock_input=bs[i,j],M_contribution_eta=mc[i,j]/21,H_contribution_eta=hc[i,j]/21,total_contribution_eta=(mc[i,j]+hc[i,j])/21,reference_imputed_eta=imp[i,j]/21,M_contribution_years=mc[i,j],H_contribution_years=hc[i,j],reference_contrast_years=contrast[i,j],chemical_status='measured' if observed[i,j] else 'imputed chemistry unresolved'))
 table(pd.DataFrame(rows),'tables/clock_chemical_decomposition.tsv');common=observed.all(0);md=mode(m,h,total,w,common);lo=mode(m,h,total,w,common,.8);hi=mode(m,h,total,w,common,1.2)
 ann=pd.read_csv(N/'00_audit/frozen_previous/01_tables/clock_CpG_contributions.tsv',sep='\t');ann=ann[ann.population=='Pooled'].set_index('CpGmarker').loc[c.CpGmarker]
 loci=pd.DataFrame(dict(CpGmarker=c.CpGmarker,chromosome=c.hg38_chromosome,start_0based=c.hg38_start,clock_weight=w,observed_patients=observed.sum(0),median_depth=np.median(depth,axis=0),minimum_depth=depth.min(0),mean_M=np.nanmean(np.where(observed,m,np.nan),0),mean_H=np.nanmean(np.where(observed,h,np.nan),0),mean_total=np.nanmean(np.where(observed,total,np.nan),0),mean_M_effect_years=mc.mean(0),mean_H_effect_years=hc.mean(0),mean_reference_contrast_years=contrast.mean(0),SD_contribution_years=np.std(contrast,axis=0),chemical_mode=md,mode_threshold_stable=(md==lo)&(md==hi),mode_lo=lo,mode_hi=hi,promoter_genes=ann.promoter_genes.to_numpy(),gene_body_genes=ann.gene_body_genes.to_numpy()))
 # All eligible locus-feature tests, reduce part-whole arithmetic by excluding tested reference contrast from outcome.
 tests=[];rng=np.random.default_rng(SEED+1)
 for j in range(len(w)):
  for feature,xx in [('total',total[:,j]),('M',m[:,j]),('H',h[:,j]),('balance',np.where(total[:,j]>=.1,m[:,j]/np.maximum(total[:,j],1e-12),np.nan))]:
   good=observed[:,j]&np.isfinite(xx);n=good.sum()
   if n<18 or np.std(xx[good])<.02:continue
   x=xx[good];y=(delta-contrast[:,j])[good];rho=float(corr(x,y));perms=np.array([rng.permutation(n) for _ in range(4999)]);null=corr(x,y[perms]);pv=(1+sum(abs(null)>=abs(rho)-1e-12))/5000;loo=np.array([corr(np.delete(x,i),np.delete(y,i)) for i in range(n)])
   tests.append(dict(CpGmarker=c.CpGmarker.iloc[j],feature=feature,n=int(n),rho=rho,raw_offset_rho=float(corr(x,delta[good])),p=pv,loo_min=float(np.nanmin(loo)),loo_max=float(np.nanmax(loo)),min_abs_loo=float(np.nanmin(abs(loo))),loo_same_sign=bool(np.all(np.sign(loo)==np.sign(rho)))))
 test=pd.DataFrame(tests);test['FDR']=bh(test.p);test['highlight']=(abs(test.rho)>=.5)&(test.FDR<=.05)&(test.min_abs_loo>=.4)&test.loo_same_sign
 # Intervals only for shortlisted effects, avoid hundreds of unsupported final claims.
 test['ci_low']=np.nan;test['ci_high']=np.nan
 for k,r in test[test.highlight].iterrows():
  j=list(c.CpGmarker).index(r.CpGmarker);xx={'total':total[:,j],'M':m[:,j],'H':h[:,j],'balance':np.where(total[:,j]>=.1,m[:,j]/np.maximum(total[:,j],1e-12),np.nan)}[r.feature];xx=np.where(observed[:,j],xx,np.nan);st,_=association(xx,delta-contrast[:,j],p.tumour_type,nperm=999);test.loc[k,['ci_low','ci_high']]=[st['ci_low'],st['ci_high']]
 table(test,'tables/locus_associations.tsv');best=test.sort_values(['FDR','p']).drop_duplicates('CpGmarker');loci=loci.merge(best[['CpGmarker','feature','rho','p','FDR','loo_min','loo_max','highlight','ci_low','ci_high']],on='CpGmarker',how='left');table(loci,'tables/clock_locus_chemical_modes.tsv')
 # Null control: one fixed independent patient permutation; same family and missingness rules.
 permutation=rng.permutation(len(p));nullps=[]
 for r in test.itertuples():
  j=list(c.CpGmarker).index(r.CpGmarker);xx={'total':total[:,j],'M':m[:,j],'H':h[:,j],'balance':np.where(total[:,j]>=.1,m[:,j]/np.maximum(total[:,j],1e-12),np.nan)}[r.feature];good=observed[:,j]&np.isfinite(xx);x=xx[good];y=(delta-contrast[:,j])[permutation][good];rr=float(corr(x,y));perms=np.array([rng.permutation(len(x)) for _ in range(999)]);nu=corr(x,y[perms]);nullps.append((1+sum(abs(nu)>=abs(rr)-1e-12))/1000)
 save(dict(tests=len(test),highlighted_tests=int(test.highlight.sum()),highlighted_loci=int(test[test.highlight].CpGmarker.nunique()),fixed_permuted_family_BH05=int(sum(bh(nullps)<=.05)),sum_to_one_max_error=float(np.nanmax(abs(u+m+h-1))),common_clock_loci=int(common.sum()),chemical_modes=pd.Series(md[common]).value_counts().to_dict(),threshold_stable_common=int(sum(loci.mode_threshold_stable&common)),mean_H_effect=float((age-age_m).mean()),H_effect_range=[float(min(age-age_m)),float(max(age-age_m))],absolute_measured_M_fraction_median=float(pp.absolute_measured_5mC_signal_fraction.median())),'qc/chemical_validation.json')
 # Patient clustering validation: common 30 highest contribution SD, independent of p.
 chosen=np.flatnonzero(common);chosen=chosen[np.argsort(-np.std(contrast[:,chosen],axis=0))[:30]];X=contrast[:,chosen].T;X=(X-X.mean(1,keepdims=True))/np.maximum(X.std(1,keepdims=True),1e-9);base=fcluster(linkage(pdist(X.T,'correlation'),method='average'),2,criterion='maxclust');support=[];co=np.zeros((20,20))
 for _ in range(500):
  ix=rng.integers(0,len(chosen),len(chosen));cl=fcluster(linkage(pdist(X[ix].T,'correlation'),method='average'),2,criterion='maxclust');cl=cl if np.mean(cl==base)>=.5 else 3-cl;support.append(cl==base);co+=cl[:,None]==cl[None,:]
 support=np.mean(support,axis=0);stable=bool(np.mean(support>=.8)>=.8)
 save(dict(loci=c.CpGmarker.iloc[chosen].tolist(),k=2,bootstrap_replicates=500,assignment_support=dict(zip(p.sample_id,map(float,support))),stable_gate=stable),'qc/landscape_stability.json');np.savez_compressed(N/'intermediate/chemical_arrays.npz',M=m,H=h,U=u,total=total,observed=observed,contrast=contrast,depth=depth,chosen=chosen,cluster=base,co_clustering=co/500)
 print('CHEMISTRY COMPLETE',len(test),'tests;',test.highlight.sum(),'highlighted',flush=True)
if __name__=='__main__':run()
