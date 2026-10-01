"""Processed-count runtime; raw-read extraction is intentionally not distributed here."""
import os
from pathlib import Path
import json
R=Path(os.environ['AGE_WORKSPACE']);D=R/'AGE_ANALYSIS/deep_dive';O=D/'06_PDR_EPY';I=D/'.intermediate/local_epialleles';SEED=20260922


def save(p,x):
    p=Path(p);t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(x,indent=2));t.replace(p)
