"""Prespecified participant-level statistics for the external five-locus study."""
import os
from pathlib import Path
import numpy as np
from scipy import stats

B = 5000
P = 19999


def ci(v):
    return np.nanquantile(v, [.025, .975]).tolist()


def bh(pvalues):
    p = np.asarray(pvalues, float)
    safe = np.where(np.isfinite(p), p, 1.)
    order = np.argsort(safe)
    adjusted = np.minimum.accumulate((safe[order] * len(p) / np.arange(1, len(p)+1))[::-1])[::-1]
    out = np.empty(len(p))
    out[order] = np.minimum(1., adjusted)
    out[~np.isfinite(p)] = np.nan
    return out


def spearman(x, y, rng):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if len(x) < 8 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return dict(n=len(x), rho=np.nan, rho_CI_low=np.nan, rho_CI_high=np.nan, p=np.nan)
    rx, ry = stats.rankdata(x), stats.rankdata(y)
    rx -= rx.mean(); ry -= ry.mean()
    denom = np.linalg.norm(rx) * np.linalg.norm(ry)
    rho = float(rx @ ry / denom)
    draws, extreme = [], 0
    for start in range(0, B, 200):
        idx = rng.integers(0, len(x), (min(200, B-start), len(x)))
        a, b = stats.rankdata(x[idx], axis=1), stats.rankdata(y[idx], axis=1)
        a -= a.mean(axis=1, keepdims=True); b -= b.mean(axis=1, keepdims=True)
        with np.errstate(invalid='ignore', divide='ignore'):
            draws.extend(np.sum(a*b, axis=1) / np.sqrt(np.sum(a*a, axis=1)*np.sum(b*b, axis=1)))
    for start in range(0, P, 200):
        yp = rng.permuted(np.tile(ry, (min(200, P-start), 1)), axis=1)
        extreme += np.count_nonzero(np.abs(yp @ rx / denom) >= abs(rho) - 1e-14)
    lo, hi = ci(draws)
    return dict(n=len(x), rho=rho, rho_CI_low=lo, rho_CI_high=hi, p=(extreme+1)/(P+1))


def theil_bootstrap(x, y, rng):
    """Exact Theil-Sen refits under case bootstrap via weighted pair slopes.

    Identical-age pairs are excluded, as in scipy.stats.theilslopes.
    Case multiplicities weight each distinct original pair by c_i*c_j.
    """
    x, y = np.asarray(x, float), np.asarray(y, float)
    fit = stats.theilslopes(y, x, method='joint')
    i, j = np.triu_indices(len(x), 1)
    valid = x[i] != x[j]
    i, j = i[valid], j[valid]
    slopes = (y[j]-y[i])/(x[j]-x[i])
    order = np.argsort(slopes)
    slopes, i, j = slopes[order], i[order], j[order]
    draws = np.empty((B, 2))
    for b in range(B):
        idx = rng.integers(0, len(x), len(x))
        counts = np.bincount(idx, minlength=len(x))
        cum = np.cumsum(counts[i]*counts[j])
        if cum[-1] == 0:
            draws[b] = np.nan
            continue
        lower = np.searchsorted(cum, (cum[-1]+1)//2)
        upper = np.searchsorted(cum, (cum[-1]+2)//2)
        slope = (slopes[lower]+slopes[upper])/2
        draws[b] = (slope, np.median(y[idx]-slope*x[idx]))
    return float(fit.slope), float(fit.intercept), draws


def group_difference(tumour, control, rng):
    t, c = np.asarray(tumour, float), np.asarray(control, float)
    t, c = t[np.isfinite(t)], c[np.isfinite(c)]
    if len(t) < 2 or len(c) < 2:
        return dict(n_GBM=len(t), n_control=len(c), delta_beta=np.nan,
                    delta_CI_low=np.nan, delta_CI_high=np.nan, p=np.nan)
    delta = t.mean()-c.mean()
    boot = np.mean(rng.choice(t, (B, len(t))), axis=1)-np.mean(rng.choice(c, (B, len(c))), axis=1)
    pooled = np.r_[t, c]
    extreme = 0
    for start in range(0, P, 200):
        v = rng.permuted(np.tile(pooled, (min(200, P-start), 1)), axis=1)
        d = v[:, :len(t)].mean(axis=1)-v[:, len(t):].mean(axis=1)
        extreme += np.count_nonzero(abs(d) >= abs(delta)-1e-14)
    lo, hi = ci(boot)
    return dict(n_GBM=len(t), n_control=len(c), delta_beta=delta,
                delta_CI_low=lo, delta_CI_high=hi, p=(extreme+1)/(P+1))


def trajectory_deviation(age, beta, normal_age, slope, intercept, normal_boot, rng):
    age, beta = np.asarray(age, float), np.asarray(beta, float)
    ok = np.isfinite(age) & np.isfinite(beta) & (age >= min(normal_age)) & (age <= max(normal_age))
    age, beta = age[ok], beta[ok]
    if not len(age):
        return dict(n=0, delta=np.nan, low=np.nan, high=np.nan)
    idx = rng.integers(0, len(age), (B, len(age)))
    draws = (beta[idx] - normal_boot[:, 0, None]*age[idx] - normal_boot[:, 1, None]).mean(axis=1)
    lo, hi = ci(draws)
    return dict(n=len(age), delta=float(np.mean(beta-slope*age-intercept)), low=lo, high=hi)


def inv_age(eta):
    eta = np.asarray(eta)
    return np.where(eta < 0, 21*np.exp(np.minimum(eta, 0))-1, 21*eta+20)


def adjusted_rank_correlation(age, beta, covariates, rng):
    """Rank correlation after linear adjustment for supplied nuisance terms.
    CI: participant bootstrap with ranks and nuisance fits recomputed.
    P: approximate partial-correlation t test; sensitivity analysis only.
    """
    x, y, z = np.asarray(age, float), np.asarray(beta, float), np.asarray(covariates, float)
    good = np.isfinite(x) & np.isfinite(y) & np.isfinite(z).all(axis=1)
    x, y, z = x[good], y[good], z[good]
    def calc(x, y, z):
        a = np.column_stack([np.ones(len(x)), z])
        ranks = np.column_stack([stats.rankdata(x), stats.rankdata(y)])
        res = ranks-a @ np.linalg.lstsq(a, ranks, rcond=None)[0]
        return np.corrcoef(res.T)[0,1]
    rho = calc(x,y,z)
    draws = []
    for _ in range(B):
        idx = rng.integers(0,len(x),len(x))
        draws.append(calc(x[idx],y[idx],z[idx]))
    lo, hi = ci(draws)
    df = len(x)-np.linalg.matrix_rank(np.column_stack([np.ones(len(x)),z]))-1
    p = 2*stats.t.sf(abs(rho)*np.sqrt(df/(1-rho*rho)),df)
    return dict(n=len(x), rho=rho, rho_CI_low=lo, rho_CI_high=hi, p=p)


def robust_age_group(tage, tbeta, cage, cbeta, rng):
    """Huber IRLS beta~1+age/10+tumour, stratified participant bootstrap CI."""
    ta,tb,ca,cb = map(lambda x:np.asarray(x,float),(tage,tbeta,cage,cbeta))
    tg, cg = np.isfinite(ta)&np.isfinite(tb), np.isfinite(ca)&np.isfinite(cb)
    ta,tb,ca,cb = ta[tg],tb[tg],ca[cg],cb[cg]
    if len(ta)<8 or len(ca)<3:
        return dict(n_GBM=len(ta), n_control=len(ca), delta=np.nan, low=np.nan, high=np.nan)
    def fit(ta,tb,ca,cb):
        a = np.column_stack([np.ones(len(ta)+len(ca)),(np.r_[ta,ca]-50)/10,np.r_[np.ones(len(ta)),np.zeros(len(ca))]])
        y = np.r_[tb,cb]
        coef = np.linalg.lstsq(a,y,rcond=None)[0]
        for _ in range(50):
            resid = y-a@coef
            scale = np.median(abs(resid-np.median(resid)))/.6744897501960817
            if scale<1e-12:break
            w=np.sqrt(np.minimum(1,1.345*scale/np.maximum(abs(resid),1e-15)))
            new = np.linalg.lstsq(a*w[:,None],y*w,rcond=None)[0]
            if np.max(abs(new-coef))<1e-8:
                coef=new;break
            coef=new
        return coef[-1]
    delta=fit(ta,tb,ca,cb)
    draws=[]
    for _ in range(B):
        ti=rng.integers(0,len(ta),len(ta));ci_=rng.integers(0,len(ca),len(ca))
        draws.append(fit(ta[ti],tb[ti],ca[ci_],cb[ci_]))
    low,high=ci(draws)
    return dict(n_GBM=len(ta),n_control=len(ca),delta=delta,low=low,high=high)
