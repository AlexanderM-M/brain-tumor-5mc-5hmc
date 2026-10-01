import os
from pathlib import Path
from common import *
from scipy.sparse import csr_matrix
import concurrent.futures

def metrics(z,n,repeats,rng):
 # Identical molecule subsets and fixed CpGs for ternary and binary encodings.
 ix=np.argsort(rng.random((repeats,len(z))),axis=1)[:,:n];a=z[ix];code3=a@np.array([1,3,9]);code2=(a>0)@np.array([1,2,4]);freq3=np.bincount((code3+27*np.arange(repeats)[:,None]).ravel(),minlength=repeats*27).reshape(repeats,27)/n;freq2=np.bincount((code2+8*np.arange(repeats)[:,None]).ravel(),minlength=repeats*8).reshape(repeats,8)/n
 h3=-(freq3*np.log2(np.maximum(freq3,1e-300))).sum(1);h2=-(freq2*np.log2(np.maximum(freq2,1e-300))).sum(1);assert min(h3-h2)>-1e-12
 return np.array([h3.mean(),h2.mean(),(h3-h2).mean(),freq3.max(1).mean(),freq2.max(1).mean()])

def sample(label):
 if not globals().get('FORCE',False) and (N/f'qc/{label}_sampling.json').exists() and (N/f'intermediate/{label}_triplet_metrics.tsv').exists() and (N/f'intermediate/{label}_age_bootstrap.npz').exists() and 'triplet_H_fraction' in (N/f'intermediate/{label}_triplet_metrics.tsv').open().readline():return
 p,c,ref,counts=inputs();i=list(p.sample_id).index(label);w=c.CoefficientTraining.to_numpy();d=counts[i].sum(1);obs=d>=5;v=np.divide(counts[i],d[:,None],out=np.zeros(counts[i].shape,float),where=d[:,None]>0);baseM=np.where(obs,v[:,1],ref);baseB=np.where(obs,v[:,1]+v[:,2],ref)
 q=np.load(N/f'intermediate/{label}_molecules.npz');idx=q['indices'];row=q['row'];col=q['col'];state=q['state'];nread=len(q['read_hash']);clockcol=np.searchsorted(idx,c.reference_index.to_numpy());mapping=np.full(len(idx),-1);mapping[clockcol]=np.arange(353);good=mapping[col]>=0;rr=row[good];cc=mapping[col[good]];ss=state[good];uniq,rr=np.unique(rr,return_inverse=True);K=len(uniq);den=csr_matrix((np.ones(len(rr)),(rr,cc)),shape=(K,353));M=csr_matrix(((ss==1).astype(float),(rr,cc)),shape=(K,353));H=csr_matrix(((ss==2).astype(float),(rr,cc)),shape=(K,353));assert np.array_equal(np.asarray(den.sum(0)).ravel(),d)
 rng=np.random.default_rng(SEED+100+i);out=[];thin=[];zero=0;zero_thin=0
 for procedure,reps in [('bootstrap',2000),('downsample',500)]:
  results=[]
  for start in range(0,reps,100):
   weights=rng.poisson(1,(100,K)) if procedure=='bootstrap' else (rng.random((100,K))<min(1,10/np.median(d[obs])))
   dd=den.T.dot(weights.T).T;mm=M.T.dot(weights.T).T;hh=H.T.dot(weights.T).T;bm=np.divide(mm,dd,out=np.tile(v[:,1],(100,1)),where=dd>0);bh=np.divide(hh,dd,out=np.tile(v[:,2],(100,1)),where=dd>0);aM=inv(.696+np.where(obs,bm,ref)@w);aB=inv(.696+np.where(obs,bm+bh,ref)@w);results.append(np.column_stack([aM,aB,aB-aM]));zeros=int(((dd==0)&obs).sum())
   if procedure=='bootstrap':zero+=zeros
   else:zero_thin+=zeros
  if procedure=='bootstrap':out=np.vstack(results)
  else:thin=np.vstack(results)
 np.savez_compressed(N/f'intermediate/{label}_age_bootstrap.npz',bootstrap=out,downsample=thin)
 trip=pd.read_csv(N/'tables/molecule_triplet_definition.tsv',sep='\t');matrix=csr_matrix((state+1,(row,col)),shape=(nread,len(idx)));rows=[]
 for t in trip.itertuples():
  cols=np.searchsorted(idx,list(map(int,t.indices.split(','))));x=matrix[:,cols].toarray()-1;x=x[(x>=0).all(1)];n=len(x);record=dict(sample_id=label,CpGmarker=t.CpGmarker,complete_molecules=n,triplet_H_fraction=float(np.mean(x==2)) if n else np.nan)
  if n>=12:
   vals=metrics(x,12,256,rng);v8=metrics(x,8,256,rng)
   for k,a in zip(['entropy3_bits','entropyBS_bits','hidden_entropy_bits','major3_fraction','majorBS_fraction'],vals):record[k]=float(a)
   record.update(entropy3_normalised=vals[0]/np.log2(27),entropyBS_common_normalisation=vals[1]/np.log2(27),hidden_entropy8_bits=v8[2],entropy3_8_bits=v8[0])
   if n>=16:
    v16=metrics(x,16,256,rng);record.update(hidden_entropy16_bits=v16[2],entropy3_16_bits=v16[0])
  rows.append(record)
 table(pd.DataFrame(rows),f'intermediate/{label}_triplet_metrics.tsv');save(dict(sample_id=label,unique_clock_molecules=K,molecules_linking_multiple_clock_CpGs=int(sum(np.asarray((den>0).sum(1)).ravel()>1)),bootstrap_zero_denominators=zero,bootstrap_observed_locus_replicates=2000*int(obs.sum()),downsample_zero_denominators=zero_thin,downsample_observed_locus_replicates=500*int(obs.sum()),bootstrap_replicates=2000,downsample_probability=float(min(1,10/np.median(d[obs])))) ,f'qc/{label}_sampling.json')
 print(label,'sampling/entropy complete',flush=True)

def summarize():
 p,c,ref,counts=inputs();a=np.stack([np.load(N/f'intermediate/{l}_age_bootstrap.npz')['bootstrap'] for l in p.sample_id]);thin=np.stack([np.load(N/f'intermediate/{l}_age_bootstrap.npz')['downsample'] for l in p.sample_id]);base=pd.read_csv(N/'tables/patient_base.tsv',sep='\t')
 for j,k in enumerate(['Age_M','Age_BS','hydroxymethylation_clock_effect']):
  base[k+'_CI_low']=np.quantile(a[:,:,j],.025,axis=1);base[k+'_CI_high']=np.quantile(a[:,:,j],.975,axis=1);base[k+'_downsample_mean']=thin[:,:,j].mean(1)
 base['DeltaAge_CI_low']=base.Age_BS_CI_low-base.chronological_age;base['DeltaAge_CI_high']=base.Age_BS_CI_high-base.chronological_age;base['offset_CI_above_zero']=base.DeltaAge_CI_low>0
 rng=np.random.default_rng(SEED+200);draw=draws(p.tumour_type,10000,rng);chosen=rng.integers(0,2000,draw.shape);means=(a[draw,chosen,1]-p.chronological_age.to_numpy()[draw]).mean(1);within=(a[:,:,1]-p.chronological_age.to_numpy()[:,None]).mean(0)
 allrows=pd.concat([pd.read_csv(N/f'intermediate/{l}_triplet_metrics.tsv',sep='\t') for l in p.sample_id]);depth=allrows.pivot(index='CpGmarker',columns='sample_id',values='complete_molecules').reindex(columns=p.sample_id);common=depth.index[(depth>=12).all(1)].tolist();valid=len(common)>=10;table(depth.reset_index(),'qc/molecule_complete_read_coverage.tsv');table(allrows,'tables/molecule_locus_sample_metrics.tsv');mr=[]
 if valid:
  for label in p.sample_id:
   s=allrows[(allrows.sample_id==label)&allrows.CpGmarker.isin(common)];assert len(s)==len(common)
   r=dict(sample_id=label,common_triplets=len(common));r.update({k:float(s[k].mean()) for k in ['entropy3_bits','entropyBS_bits','hidden_entropy_bits','major3_fraction','majorBS_fraction','entropy3_normalised','entropyBS_common_normalisation','hidden_entropy8_bits','entropy3_8_bits','triplet_H_fraction']});mr.append(r)
  molecule=pd.DataFrame(mr);table(molecule,'tables/molecule_level_summary.tsv');base=base.merge(molecule,on='sample_id');assoc=[];ll=[]
  for feature in ['entropy3_bits','hidden_entropy_bits','major3_fraction']:
   for target in ['DNAm_age_acceleration','tumour_content','chronological_age']:
    st,loo=association(base[feature],base[target],base.tumour_type);assoc.append(dict(feature=feature,target=target,**st));ll.extend(dict(feature=feature,target=target,omitted_patient=l,rho=float(x)) for l,x in zip(p.sample_id,loo))
  st,loo=association(base.hidden_entropy_bits,base.DNAm_age_acceleration,base.tumour_type,adjust=base.triplet_H_fraction);assoc.append(dict(feature='hidden_entropy_bits',target='offset adjusted for triplet H fraction',**st));ll.extend(dict(feature='hidden_entropy_bits_adjusted_H',target='DNAm_age_acceleration',omitted_patient=l,rho=float(x)) for l,x in zip(p.sample_id,loo))
  ast=pd.DataFrame(assoc);ast['FDR']=bh(ast.p);table(ast,'tables/molecule_associations.tsv');table(pd.DataFrame(ll),'qc/molecule_leave_one_out.tsv')
  comparison=[]
  for feature in ['entropy3','hidden_entropy']:
   st,_=association(base[feature+'_bits'],base[feature+'8_bits'] if feature=='hidden_entropy' else base.entropy3_8_bits,base.tumour_type,nperm=999);comparison.append(dict(feature=feature,comparison='12 vs 8 reads',**st))
   st,_=association(base[feature+'8_bits'] if feature=='hidden_entropy' else base.entropy3_8_bits,base.DNAm_age_acceleration,base.tumour_type,nperm=9999);comparison.append(dict(feature=feature,comparison='8 reads vs offset',**st))
  # 16-read sensitivity uses its own common panel; compare 12/16 on those exact loci.
  common16=depth.index[(depth>=16).all(1)].tolist()
  if len(common16)>=10:
   for feature in ['entropy3','hidden_entropy']:
    pair=[]
    for label in p.sample_id:
     s=allrows[(allrows.sample_id==label)&allrows.CpGmarker.isin(common16)];pair.append([s[feature+'_bits'].mean(),s[feature+'16_bits' if feature=='hidden_entropy' else 'entropy3_16_bits'].mean()])
    pair=np.array(pair);st,_=association(pair[:,0],pair[:,1],p.tumour_type,nperm=999);comparison.append(dict(feature=feature,comparison='12 vs 16 reads; same restricted loci',**st));st,_=association(pair[:,1],p.DNAm_age_acceleration,p.tumour_type,nperm=9999);comparison.append(dict(feature=feature,comparison='16 reads vs offset',**st))
  table(pd.DataFrame(comparison),'qc/entropy_depth_sensitivity.tsv')
  exemplar=allrows[allrows.CpGmarker.isin(common)].groupby('CpGmarker').hidden_entropy_bits.median().sort_values(ascending=False,kind='stable').head(2).index.tolist();sample_id=p.iloc[np.argsort(p.DNAm_age_acceleration.to_numpy(),kind='stable')[9]].sample_id
 else:exemplar=[];sample_id=None;common16=[]
 table(base,'tables/age_patient_summary.tsv');qs=[json.loads((N/f'qc/{l}_sampling.json').read_text()) for l in p.sample_id];save(dict(primary_point_positive=int(sum(p.DNAm_age_acceleration>0)),sampling_CI_positive=int(sum(base.offset_CI_above_zero)),mean_offset=float(p.DNAm_age_acceleration.mean()),cohort_patient_and_molecule_CI=list(map(float,np.quantile(means,[.025,.975]))),conditional_fixed_cohort_CI=list(map(float,np.quantile(within,[.025,.975]))),median_individual_CI_width=float(np.median(base.Age_BS_CI_high-base.Age_BS_CI_low)),common_triplets=len(common),common_triplets16=len(common16),molecule_analysis_valid=valid,common_triplet_ids=common,exemplar_loci=exemplar,exemplar_patient=sample_id,median_downsample_age_shift=float(np.median(base.Age_BS_downsample_mean-base.Age_BS)),max_abs_downsample_mean_shift=float(max(abs(base.Age_BS_downsample_mean-base.Age_BS))),total_multiclock_molecules=sum(x['molecules_linking_multiple_clock_CpGs'] for x in qs),total_clock_molecules=sum(x['unique_clock_molecules'] for x in qs),bootstrap_zero_fraction=sum(x['bootstrap_zero_denominators'] for x in qs)/sum(x['bootstrap_observed_locus_replicates'] for x in qs),downsample_zero_fraction=sum(x['downsample_zero_denominators'] for x in qs)/sum(x['downsample_observed_locus_replicates'] for x in qs)),'qc/molecule_analysis_validation.json')
 print('MOLECULE ANALYSIS VALID:',valid,'common triplets',len(common),flush=True)

def run():
 p,*_=inputs()
 with concurrent.futures.ProcessPoolExecutor(max_workers=4) as pool:list(pool.map(sample,p.sample_id))
 summarize()
if __name__=='__main__':run()
