"""All test data are generated synthetically in pytest's temporary directory."""

import pytest
import pysam
from brain_5mc_5hmc import runtime


@pytest.fixture(scope="session", autouse=True)
def workspace(tmp_path_factory):
    root = tmp_path_factory.mktemp("synthetic_workspace")
    fasta = root / "synthetic.fa"
    sequence = "A" * 150 + "CG" + "A" * 248
    fasta.write_text(
        "".join(f">chr{i}\n{sequence if i == 1 else 'A' * 400}\n" for i in range(1, 23))
    )
    pysam.faidx(str(fasta))
    runtime.configure(workspace=root, reference=fasta, workers=1)
    runtime.prepare_workspace()
    runtime.build_native()
    for stage in [
        "pooled_analysis",
        "story_molecules",
        "followup_molecules",
        "followup_select",
        "followup_distance",
        "followup_row_null",
        "strengthen_score",
        "strengthen_permutation",
        "deep_regional",
    ]:
        runtime.initialize(stage)
    return root
