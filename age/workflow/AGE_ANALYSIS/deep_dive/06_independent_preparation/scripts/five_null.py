"""Prespecified conditional matched-five framework; no WGS dependencies."""
import os
from pathlib import Path

from resources import constrain

constrain()
import numpy as np

SEED = 20260917
N_NULL = 10000


def inverse_age(eta):
    eta = np.asarray(eta, dtype=float)
    return np.where(eta < 0, 21 * np.exp(np.minimum(eta, 0)) - 1, 21 * eta + 20)


def replacement_effect(beta, reference, weights, selected, intercept=0.696):
    beta, reference, weights = map(
        lambda x: np.asarray(x, float), (beta, reference, weights)
    )
    if (
        beta.ndim != 2
        or beta.shape[1:] != reference.shape
        or weights.shape != reference.shape
    ):
        raise ValueError("Incompatible feature dimensions")
    if not all(np.isfinite(x).all() for x in (beta, reference, weights)):
        raise ValueError("Impute using frozen reference BEFORE counterfactual scoring")
    if ((beta < 0) | (beta > 1)).any() or ((reference < 0) | (reference > 1)).any():
        raise ValueError("Beta outside [0,1]")
    ix = np.asarray(selected, int)
    if len(set(ix)) != len(ix) or (ix < 0).any() or (ix >= len(weights)).any():
        raise ValueError("Invalid selected features")
    full_eta = intercept + beta @ weights
    shift_eta = ((beta[:, ix] - reference[ix]) * weights[ix]).sum(axis=1)
    counter_eta = full_eta - shift_eta
    return dict(
        full_eta=full_eta,
        score_contribution=shift_eta,
        full_age=inverse_age(full_eta),
        counterfactual_age=inverse_age(counter_eta),
        contribution_years=inverse_age(full_eta) - inverse_age(counter_eta),
    )


def matching_pools(features, selected, context="strict"):
    """Fixed calipers; strict context. Insufficient matching is a result, not relaxed away.

    features: list of dictionaries (id, weight, median_depth, eligible,
    island_context, genic_context). No clock outcome or methylation level used.
    """
    if context not in ("strict", "island_only", "calipers_only"):
        raise ValueError("Unknown matching profile")
    by_id = {x["id"]: x for x in features}
    if len(by_id) != len(features) or len(set(selected)) != 5:
        raise ValueError("Require unique features and exactly five targets")
    pools = {}
    for target in selected:
        t = by_id[target]
        if not t["eligible"] or t["weight"] == 0 or t["median_depth"] <= 0:
            raise ValueError("Target fails fixed eligibility")
        pool = []
        for c in features:
            if not c["eligible"] or c["id"] in selected:
                continue
            wr = abs(c["weight"] / t["weight"])
            dr = c["median_depth"] / t["median_depth"]
            same_context = context == "calipers_only" or (
                c["island_context"] == t["island_context"]
                and (
                    context == "island_only" or c["genic_context"] == t["genic_context"]
                )
            )
            if (
                c["weight"] * t["weight"] > 0
                and 0.8 <= wr <= 1.25
                and 2 / 3 <= dr <= 1.5
                and same_context
            ):
                pool.append(c["id"])
        pools[target] = sorted(pool)
    return pools


def draw_matched_sets(pools, n=N_NULL, seed=SEED, max_attempts=1000000):
    """Independent target draws then reject collisions: uniform feasible assignments.

    Avoid sequential greedy draws, which bias later target choices. Replicate
    sets may recur; no within-set repeated locus. Empty/unmatchable pools fail.
    """
    if len(pools) != 5 or n < 1 or any(not p for p in pools.values()):
        raise ValueError("Five nonempty pools required")
    if any(len(set(p)) != len(p) for p in pools.values()):
        raise ValueError("Duplicate pool entry")
    from itertools import combinations

    # Hall's condition detects impossible collision-free assignments up front.
    lists = list(pools.values())
    for k in range(1, 6):
        for subset in combinations(lists, k):
            if len(set().union(*map(set, subset))) < k:
                raise ValueError("No feasible five-distinct-locus matching")
    rng = np.random.default_rng(seed)
    result = []
    attempts = 0
    while len(result) < n and attempts < max_attempts:
        candidate = tuple(rng.choice(p) for p in lists)
        attempts += 1
        if len(set(candidate)) == 5:
            result.append(candidate)
    if len(result) != n:
        raise RuntimeError(
            "Matching rejection budget exhausted; do not report partial null"
        )
    return result


def conditional_null(beta, reference, weights, ids, selected, sets):
    """Mean patient-level age change, same exact counterfactual for every set."""
    lookup = {x: i for i, x in enumerate(ids)}
    effect = replacement_effect(beta, reference, weights, [lookup[x] for x in selected])
    observed = float(effect["contribution_years"].mean())
    null = np.array(
        [
            replacement_effect(beta, reference, weights, [lookup[x] for x in s])[
                "contribution_years"
            ].mean()
            for s in sets
        ]
    )
    if not len(null):
        raise ValueError("Empty null")
    return dict(
        observed_mean_years=observed,
        null=null,
        conditional_tail=float((1 + np.sum(null >= observed)) / (1 + len(null))),
        selection_adjusted=False,
    )
