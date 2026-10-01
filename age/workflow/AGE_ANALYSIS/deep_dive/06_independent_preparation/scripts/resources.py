"""Resource policy shared by every entry point; no workflow/process manipulation."""

import os
import sys


def constrain():
    allowed = os.sched_getaffinity(0) if hasattr(os, "sched_getaffinity") else set()
    for name in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
        "BLIS_NUM_THREADS",
        "OMP_THREAD_LIMIT",
    ):
        os.environ[name] = "1"
    os.environ["PYTHONDONTWRITEBYTECODE"] = "1"
    sys.dont_write_bytecode = True
    return sorted(allowed)


constrain()
