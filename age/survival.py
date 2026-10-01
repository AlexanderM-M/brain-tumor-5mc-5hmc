"""Reproduce GBM Cox estimates from pseudonyms and precomputed durations only."""
from pathlib import Path
import argparse
import csv
import json
import re
from collections import Counter
import numpy as np
from survival_core import fit,arrays,validate_math


def load(path):
    with path.open() as f:
        rows=list(csv.DictReader(f,delimiter='\t'))
    required={'sample_id','diagnosis','DNAm_age_acceleration','post_sampling_survival_days','event'}
    if not rows or not required.issubset(rows[0]):
        raise ValueError('Required columns: '+', '.join(sorted(required)))
    # Accept only the public input contract, never a clinical QC export.
    if set(rows[0])-required:
        raise ValueError('Unexpected fields: use only the five de-identified input columns.')
    assert len({r['sample_id'] for r in rows})==len(rows)
    for r in rows:
        assert re.fullmatch(r'GBM-\d{2}',r['sample_id'])
        assert r['diagnosis']=='Glioblastoma' and r['event'] in ('0','1')
        assert np.isfinite(float(r['DNAm_age_acceleration']))
        assert np.isfinite(float(r['post_sampling_survival_days'])) and float(r['post_sampling_survival_days'])>=0
    return rows


def save(path,rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n')
        w.writeheader();w.writerows({k:('NA' if v is None else v) for k,v in r.items()} for r in rows)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,required=True,help='Primary cohort: five columns, durations in days')
    p.add_argument('--sensitivity-earlier',type=Path,help='All-GBM durations under earlier candidate date; no exact dates')
    p.add_argument('--sensitivity-later',type=Path,help='All-GBM durations under later candidate date; no exact dates')
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();rows=load(a.input);a.output.mkdir(parents=True,exist_ok=True)
    t,e,x=arrays(rows);result=fit(t,e,x)
    if 'HR' not in result:raise RuntimeError('Standard Cox unstable; do not force an estimate: '+result['status'])
    validate_math(t,e,x,result)
    primary=[dict(analysis='primary',n=len(rows),deaths=int(e.sum()),censored=int((1-e).sum()),
                  median_followup_censored_days=float(np.median(t[e==0])) if (e==0).any() else None,**result)]
    for label,path in [('earlier_candidate',a.sensitivity_earlier),('later_candidate',a.sensitivity_later)]:
        if path:
            r=load(path);tt,ee,xx=arrays(r);fr=fit(tt,ee,xx)
            primary.append({k:dict(analysis=label,n=len(r),deaths=int(ee.sum()),censored=int((1-ee).sum()),
                                   median_followup_censored_days=float(np.median(tt[ee==0])) if (ee==0).any() else None,**fr).get(k) for k in primary[0]})
    save(a.output/'survival_primary_GBM.tsv',primary)
    omissions=[]
    for patient in rows:
        subset=[r for r in rows if r['sample_id']!=patient['sample_id']];tt,ee,xx=arrays(subset);r=fit(tt,ee,xx)
        omissions.append(dict(omitted_patient=patient['sample_id'],n=len(subset),deaths=int(ee.sum()),
                              HR=r.get('HR'),CI95_low=r.get('CI95_low'),CI95_high=r.get('CI95_high'),
                              status=r['status'],direction_consistent=bool(np.sign(r['log_HR'])==np.sign(result['log_HR'])) if 'log_HR' in r else None))
    save(a.output/'survival_leave_one_out.tsv',omissions)
    rng=np.random.default_rng(20261001);estimates=[];failures=Counter()
    for _ in range(2000):
        ix=rng.integers(0,len(rows),len(rows));r=fit(t[ix],e[ix],x[ix])
        if 'log_HR' in r:estimates.append(r['log_HR'])
        else:failures[r['status']]+=1
    success=len(estimates);ok=success>=500 and (2000-success)/2000<=.10
    diagnostic=dict(seed=20261001,resamples=2000,successful=success,failures=dict(failures),
                    failure_fraction=(2000-success)/2000,CI_released=ok,
                    HR_CI95=np.exp(np.quantile(estimates,[.025,.975])).tolist() if ok else None)
    (a.output/'bootstrap_diagnostics.json').write_text(json.dumps(diagnostic,indent=2)+'\n')


if __name__=='__main__':main()
