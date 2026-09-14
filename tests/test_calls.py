"""Synthetic BAM checks for explicit MM/ML tags, strand and duplicate handling."""

from array import array
import json
import numpy as np
import pysam
from numpy.testing import assert_array_equal
from brain_5mc_5hmc import extraction, molecules


def read(name, m, h, reverse=False):
    r = pysam.AlignedSegment()
    r.query_name = name
    r.query_sequence = "A" * 150 + "CG" + "A" * 248
    r.flag = 16 if reverse else 0
    r.reference_id = 0
    r.reference_start = 0
    r.mapping_quality = 60
    r.cigar = [(0, 400)]
    r.set_tag("qs", 12.0)
    r.set_tag("MM", "C+mh,0;")
    r.set_tag("ML", array("B", [m, h]))
    return r


def test_probability_threshold_and_reverse_strand():
    for reverse in [False, True]:
        for m, h, state in [
            (0, 0, 0),
            (230, 0, 1),
            (0, 230, 2),
            (100, 100, -1),
            (204, 0, -1),
            (205, 0, 1),
        ]:
            r = read("synthetic", m, h, reverse)
            assert molecules.blocks_eligible(r)
            assert_array_equal(molecules.blocks_states_at(r, np.array([150])), [state])
    r = read("missing_h", 230, 0)
    r.set_tag("MM", "C+m,0;")
    r.set_tag("ML", array("B", [230]))
    assert_array_equal(molecules.blocks_states_at(r, np.array([150])), [-1])


def test_native_counts_against_scalar_oracle_and_duplicate_ingestion(workspace):
    path = workspace / "synthetic.bam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": f"chr{i}", "LN": 400} for i in range(1, 23)]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for i, (m, h) in enumerate([(0, 0), (230, 0), (0, 230), (100, 100)]):
            out.write(read(f"synthetic_{i}", m, h))
        out.write(read("synthetic_reverse", 0, 230, True))
        excluded = read("synthetic_low_mapq", 0, 230)
        excluded.mapping_quality = 10
        out.write(excluded)
    extraction.counts_validate(path, "synthetic_validation")
    result = json.loads((workspace / ".analysis/validation/synthetic_validation.json").read_text())
    assert result["native_scalar_oracle_exact_match"]
    assert result["four_threshold_region_match"]
    assert result["duplicate_reads_removed_on_repeat"] == 5
    assert result["totals"][:4] == [1, 1, 2, 1]
