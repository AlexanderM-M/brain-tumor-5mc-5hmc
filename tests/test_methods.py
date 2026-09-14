"""Scientific invariants: weighting, selection, missingness and patient inference."""

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
from brain_5mc_5hmc import regional, promoters, molecules, association, robustness


def test_strand_aware_promoters():
    assert regional.regional_promoter(3000, 9000, "+", 10000) == (1000, 3500)
    assert regional.regional_promoter(3000, 9000, "-", 10000) == (8500, 10000)
    assert regional.regional_promoter(100, 9000, "+", 10000) == (0, 600)


def test_equal_cpg_weighting():
    # CpG 1 has 10 calls; CpG 2 has 100. Each CpG still receives equal weight.
    counts = np.array([[0, 8, 2], [90, 0, 10]])
    fractions = counts[:, 1:] / counts.sum(1)[:, None] * 100
    windows, features = regional.regional_aggregate(
        fractions, np.array([0, 0]), np.array([2]), np.array([0]), np.array([2])
    )
    assert_allclose(windows, [[40, 15]])
    assert_allclose(features, [[40, 15]])
    assert not np.allclose(windows[0], counts[:, 1:].sum(0) / counts.sum() * 100)


def test_opposing_change_and_bootstrap_reproducibility():
    group = np.array([True] * 3 + [False] * 3)
    beta = np.tile([60.0, 5.0], (6, 2, 1))
    beta[:3, 0] = [63.0, 2.0]
    rows = regional.regional_result_rows([{}, {}], beta, group)
    assert rows[0]["opposing_changes_with_small_combined_difference"]
    assert rows[0]["GBM_minus_meningioma_combined_pp"] == 0
    assert not rows[1]["exploratory_5hmC_candidate"]
    a = regional.regional_bootstrap_intervals(beta, group, replicates=100)
    b = regional.regional_bootstrap_intervals(beta, group, replicates=100)
    assert_array_equal(a, b)
    assert_allclose(a[0], [-3, -3])


def test_selection_excludes_heldout_values_and_uses_separate_controls():
    features = [
        dict(
            feature_id=f"f{i}", chromosome="chr1", start=i * 200000, common_CpGs=50, GC_fraction=0.5
        )
        for i in range(4)
    ]
    beta = np.tile([60.0, 5.0], (8, 4, 1))
    group = np.array([True] * 4 + [False] * 4)
    beta[:4, :2] = [63.0, 2.0]
    train = list(range(1, 8))
    pairs = promoters.selection_select(features, beta, beta, group, train)
    assert [(p["candidate_index"], p["comparison_index"]) for p in pairs] == [(0, 2), (1, 3)]
    changed = beta.copy()
    changed[0] = 999
    assert promoters.selection_select(features, changed, changed, group, train) == pairs
    assert len({p["comparison_index"] for p in pairs}) == len(pairs)


def test_block_metric_matches_optimized_calculation():
    x = np.random.default_rng(7).integers(0, 3, (25, 3))
    original = molecules.blocks_patterns_metrics(x, np.random.default_rng(42))[3] * 100
    optimized = molecules.followup_blocks_joint_excess(x, np.random.default_rng(42))
    assert_allclose(original, optimized, atol=1e-12)


def test_row_trades_preserve_both_margins_and_missingness():
    rng = np.random.default_rng(8)
    x = rng.integers(0, 2, (40, 30), dtype=np.int8)
    x[rng.random(x.shape) < 0.2] = -1
    before = x.copy()
    moves = association.row_null_trade(x, np.random.default_rng(42), 4000)
    assert moves > 0
    assert_array_equal(x < 0, before < 0)
    assert_array_equal((x == 1).sum(0), (before == 1).sum(0))
    assert_array_equal((x == 1).sum(1), (before == 1).sum(1))
    assert not np.array_equal(x, before)


def test_native_covariance_matches_direct_complete_pair_calculation():
    rng = np.random.default_rng(7)
    x = rng.integers(0, 3, (40, 30), dtype=np.int8)
    x[rng.random(x.shape) < 0.15] = -1
    ii, jj = np.array([0, 2, 4, 8]), np.array([1, 3, 9, 12])
    for encoding in range(3):
        out, info = association.native_binary_score(x, ii, jj, encoding, 42)
        duplicate, duplicate_info = association.native_binary_score(x, ii, jj, encoding, 42)
        assert_array_equal(out, duplicate)
        assert_array_equal(info, duplicate_info)
        assert info[2] == 1 and min(info[:2]) >= 64
        for k, (i, j) in enumerate(zip(ii, jj)):
            z = x[:, [i, j]]
            z = z[(z >= 0).all(1)]
            binary = z == 2 if encoding == 0 else z == 1 if encoding == 1 else z > 0
            expected = np.cov(binary.T, ddof=1)[0, 1] * 11 / 12 * 100
            assert_allclose(out[k, 0], expected, atol=1e-10)


def test_patient_scores_pair_identical_available_bins():
    scores = np.zeros((12, 3))
    pairs = [{"candidate_index": i, "comparison_index": i + 6} for i in range(6)]
    scores[:6] = [1, 2, 100]
    scores[6:, 2] = np.nan
    value, n = robustness.randomisation_patient(scores, pairs, 0, 2)
    assert value == 1.5 and n == 6
    value, n = robustness.randomisation_patient(scores, pairs[:5], 0, 2)
    assert np.isnan(value) and n == 5


def test_unsupported_label_allocations_remain_zero():
    values = np.arange(20, dtype=float)
    values[15:18] = np.nan
    statistic, _, n_men, _, eligible = robustness.randomisation_statistic(
        values, list(range(15, 20))
    )
    assert not eligible and n_men == 2 and statistic == 0
    assert len(robustness.randomisation_RANGES) == 36


def test_fork_worker_inherits_initialized_stage():
    from brain_5mc_5hmc import runtime

    with runtime.process_pool(max_workers=1) as pool:
        result = pool.submit(regional.regional_promoter, 3000, 9000, "+", 10000).result(timeout=10)
    assert result == (1000, 3500)
