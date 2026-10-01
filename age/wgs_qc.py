"""Validate de-identified SAVANA processed outputs, without patient sequence files."""
from pathlib import Path
import argparse
import importlib.util
import json
import re


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sample',required=True)
    p.add_argument('--folder',type=Path,required=True)
    p.add_argument('--log',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--exit-code',type=int,required=True)
    a=p.parse_args()
    if not re.fullmatch(r'(GBM|MEN)-\d{2}',a.sample):p.error('Use a pseudonymised study ID.')
    source=Path(__file__).parent/'workflow/AGE_ANALYSIS/deep_dive/04_WGS_genomic_backbone/scripts/wgs_cna_outcome.py'
    spec=importlib.util.spec_from_file_location('wgs_cna_outcome',source);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    result,_=m.validate_cna(a.sample,a.folder,a.log,a.exit_code)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
