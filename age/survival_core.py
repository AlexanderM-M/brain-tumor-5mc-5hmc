"""Retained univariable Cox kernels; inputs are de-identified durations only."""
import math
import numpy as np
from scipy.optimize import brentq,minimize_scalar
from scipy.stats import chi2,norm

def likelihood(beta, times, events, x):
    """Efron log partial likelihood, score and observed information, 1 covariate."""
    eta = beta * x
    ll = score = information = 0.0
    for t in np.unique(times[events == 1]):
        risk = times >= t
        death = (times == t) & (events == 1)
        d = int(death.sum())
        z = eta[risk].max()
        w = np.exp(eta[risk] - z)
        wd = np.exp(eta[death] - z)
        xr, xd = x[risk], x[death]
        s0, s1, s2 = w.sum(), (w*xr).sum(), (w*xr*xr).sum()
        d0, d1, d2 = wd.sum(), (wd*xd).sum(), (wd*xd*xd).sum()
        ll += eta[death].sum()
        score += xd.sum()
        for j in range(d):
            frac = j / d
            a, b, c = s0-frac*d0, s1-frac*d1, s2-frac*d2
            ll -= z + np.log(a)
            score -= b/a
            information += c/a - (b/a)**2
    return float(ll), float(score), float(information)

def fit(times, events, x):
    times, events, x = np.asarray(times,float), np.asarray(events,int), np.asarray(x,float)
    if events.sum() == 0:
        return dict(status='FAILED_ZERO_EVENTS')
    if np.ptp(x) < 1e-10:
        return dict(status='FAILED_CONSTANT_PREDICTOR')
    x = x - x.mean()
    lo, hi = -.5, .5
    while likelihood(lo,times,events,x)[1] < 0 and abs(lo) < 64:
        lo *= 2
    while likelihood(hi,times,events,x)[1] > 0 and hi < 64:
        hi *= 2
    if likelihood(lo,times,events,x)[1] <= 1e-9 or likelihood(hi,times,events,x)[1] >= -1e-9:
        return dict(status='FAILED_MONOTONE_OR_BOUNDARY_LIKELIHOOD')
    beta = brentq(lambda b: likelihood(b,times,events,x)[1],lo,hi,xtol=1e-12)
    ll, score, info = likelihood(beta,times,events,x)
    if info <= 1e-8 or abs(score) > 1e-6:
        return dict(status='FAILED_SINGULAR_INFORMATION')
    se = 1/np.sqrt(info)
    if abs(beta) > 30 or se > 100:
        return dict(status='FAILED_EXTREME_ESTIMATE')
    lr = max(0,2*(ll-likelihood(0,times,events,x)[0]))
    # Harrell training concordance; only pairs with an observed earlier death.
    pairs = [(i,j) for i in range(len(x)) if events[i] for j in range(len(x)) if times[j]>times[i]]
    concordance = np.mean([float(beta*x[i]>beta*x[j])+.5*float(beta*x[i]==beta*x[j]) for i,j in pairs]) if pairs else np.nan
    return dict(status='ESTIMATED_SPARSE_EVENTS', log_HR=beta, log_HR_SE=float(se),
                HR=float(np.exp(beta)), CI95_low=float(np.exp(beta-1.959963984540054*se)),
                CI95_high=float(np.exp(beta+1.959963984540054*se)),
                p_value=float(chi2.sf(lr,1)), p_value_method='likelihood_ratio_1df',
                Wald_p_value=float(2*norm.sf(abs(beta/se))), LR_statistic=lr,
                concordance_index=float(concordance), comparable_pairs=len(pairs),
                log_partial_likelihood=ll, score=score, information=info)

def direct_ll(beta,t,e,x):
    """Independent event/risk loop, no score/information reuse."""
    out=0.0
    for time in sorted(set(t[e==1])):
        risks=[i for i in range(len(t)) if t[i]>=time]
        deaths=[i for i in range(len(t)) if t[i]==time and e[i]==1]
        z=max(beta*x[i] for i in risks)
        sw=sum(math.exp(beta*x[i]-z) for i in risks)
        sd=sum(math.exp(beta*x[i]-z) for i in deaths)
        out+=sum(beta*x[i] for i in deaths)
        for j in range(len(deaths)):
            out-=z+math.log(sw-j*sd/len(deaths))
    return out

def validate_math(t,e,x,result):
    b=result['log_HR'];h=1e-4
    independent=minimize_scalar(lambda v:-direct_ll(v,t,e,x),bounds=(b-3,b+3),method='bounded',options={'xatol':1e-12})
    assert independent.success and abs(independent.x-b)<1e-6
    num_score=(direct_ll(b+h,t,e,x)-direct_ll(b-h,t,e,x))/(2*h)
    num_info=-(direct_ll(b+h,t,e,x)-2*direct_ll(b,t,e,x)+direct_ll(b-h,t,e,x))/h**2
    assert abs(num_score)<1e-5
    assert abs(num_info-result['information'])<1e-4
    # Analytic two-person, one-death likelihood at beta 0: -log(2), score .5, I .25.
    simple=likelihood(0,np.array([1,2]),np.array([1,0]),np.array([1.,0.]))
    assert np.allclose(simple,[-math.log(2),.5,.25])
    # Explicit Efron tie test, including derivatives away from zero.
    tt=np.array([1.,1.,2.,3.]);ee=np.array([1,1,0,1]);xx=np.array([-1.,.5,2.,0.])
    for bb in [-.4,0.,.6]:
        ll,s,i=likelihood(bb,tt,ee,xx)
        assert abs(ll-direct_ll(bb,tt,ee,xx))<1e-12
        assert abs(s-(direct_ll(bb+h,tt,ee,xx)-direct_ll(bb-h,tt,ee,xx))/(2*h))<1e-6
        assert abs(i+(direct_ll(bb+h,tt,ee,xx)-2*ll+direct_ll(bb-h,tt,ee,xx))/h**2)<1e-5
    assert fit(np.array([1.,2.]),np.array([1,0]),np.array([1.,0.]))['status'].startswith('FAILED')
    return dict(independent_likelihood_optimum_difference=float(independent.x-b),
                score_finite_difference=float(num_score), information_finite_difference=float(num_info),
                analytic_and_Efron_tie_checks='PASS', monotone_likelihood_detection='PASS')

def arrays(rows):
    return (np.array([r['post_sampling_survival_days'] for r in rows],float),
            np.array([r['event'] for r in rows],int),
            np.array([float(r['DNAm_age_acceleration'])/10 for r in rows]))

def numeric(v):
    return None if v in ('NA','') else float(v)

def describe(rows,label):
    usable=[r for r in rows if r['endpoint_usable']=='True']
    censor=[float(r['post_sampling_survival_days']) for r in usable if r['event']=='0']
    alltimes=[float(r['post_sampling_survival_days']) for r in usable]
    age=[float(r['DNAm_age_acceleration']) for r in rows]
    return dict(group=label,n=len(rows),n_usable=len(usable),deaths=sum(r['event']=='1' for r in usable),
                censored=sum(r['event']=='0' for r in usable),
                median_followup_censored_days=float(np.median(censor)) if censor else None,
                median_followup_censored_years=float(np.median(censor))/365.25 if censor else None,
                censored_followup_min_days=min(censor) if censor else None,censored_followup_max_days=max(censor) if censor else None,
                median_observed_time_days=float(np.median(alltimes)) if alltimes else None,
                observed_time_min_days=min(alltimes) if alltimes else None,observed_time_max_days=max(alltimes) if alltimes else None,
                acceleration_median_years=float(np.median(age)) if age else None,
                acceleration_q25_years=float(np.quantile(age,.25)) if age else None,
                acceleration_q75_years=float(np.quantile(age,.75)) if age else None,
                acceleration_min_years=min(age) if age else None,acceleration_max_years=max(age) if age else None)
