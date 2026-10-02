"""Install portable source into a separate processed-data workspace and run a stage.

This launcher never obtains clinical records or extracts raw sequencing files.
"""
from pathlib import Path
import argparse
import os
import shutil
import subprocess
import sys

BASE = Path(__file__).resolve().parent
STAGES = {
    'clock-chemistry': ('AGE_ANALYSIS/next_generation/scripts/chemistry.py', []),
    'clock-bootstrap': ('AGE_ANALYSIS/next_generation/scripts/molecules.py', []),
    'wgs-fit-stability': ('AGE_ANALYSIS/deep_dive/FINAL_INTEGRATION/scripts/post_wgs_analysis.py', ['fit_stability']),
    'purity-associations': ('AGE_ANALYSIS/deep_dive/FINAL_INTEGRATION/scripts/post_wgs_analysis.py', ['integrate']),
    'five-cpg': ('AGE_ANALYSIS/deep_dive/FINAL_INTEGRATION/scripts/post_wgs_analysis.py', ['five']),
    'joint-ploidy-qc': ('AGE_ANALYSIS/deep_dive/FINAL_INTEGRATION/scripts/post_wgs_analysis.py', ['ploidy_review']),
    'local-annotation': ('AGE_ANALYSIS/deep_dive/06_PDR_EPY/scripts/analyze_local_epialleles.py', ['annotate']),
    'local-metrics': ('AGE_ANALYSIS/deep_dive/06_PDR_EPY/scripts/analyze_local_epialleles.py', ['prepare']),
    'beta-setup': ('AGE_ANALYSIS/deep_dive/07_FIVE_LOCUS_BETA_VALIDATION/scripts/focused_beta.py', ['setup']),
    'beta-patient': ('AGE_ANALYSIS/deep_dive/07_FIVE_LOCUS_BETA_VALIDATION/scripts/focused_beta.py', ['sample']),
    'beta-chemistry': ('AGE_ANALYSIS/deep_dive/07_FIVE_LOCUS_BETA_VALIDATION/scripts/composition_focused.py', []),
    'beta-summary': ('AGE_ANALYSIS/deep_dive/07_FIVE_LOCUS_BETA_VALIDATION/scripts/report_focused_beta.py', []),
    'external-ingest': ('AGE_ANALYSIS/deep_dive/09_EXTERNAL_VALIDATION/scripts/read_archive.py', []),
    'external-targets': ('AGE_ANALYSIS/deep_dive/09_EXTERNAL_VALIDATION/scripts/extract_targeted_rows.py', []),
    'external-validation': ('AGE_ANALYSIS/deep_dive/09_EXTERNAL_VALIDATION/scripts/run_validation.py', []),
    'main-figures': ('AGE_ANALYSIS/deep_dive/FINAL_AGE_STORY/scripts/revise_figures.py', []),
}


def initialize(workspace):
    workspace.mkdir(parents=True, exist_ok=True)
    sources = list((BASE/'workflow').rglob('*.py'))
    for source in sources:
        target = workspace/source.relative_to(BASE/'workflow')
        if target.exists() and target.read_bytes() != source.read_bytes():
            raise RuntimeError('Workspace script differs; use a new workspace: '+str(target))
    for source in sources:
        target=workspace/source.relative_to(BASE/'workflow')
        target.parent.mkdir(parents=True,exist_ok=True)
        if not target.exists():
            shutil.copyfile(source,target)
    for folder in [
        'AGE_ANALYSIS/next_generation/tables','AGE_ANALYSIS/next_generation/qc',
        'AGE_ANALYSIS/next_generation/intermediate','AGE_ANALYSIS/next_generation/00_audit',
        'AGE_ANALYSIS/deep_dive/FINAL_INTEGRATION/tables','AGE_ANALYSIS/deep_dive/FINAL_INTEGRATION/provenance',
        'AGE_ANALYSIS/deep_dive/04_WGS_genomic_backbone/tables',
        'AGE_ANALYSIS/deep_dive/06_PDR_EPY/tables','AGE_ANALYSIS/deep_dive/06_PDR_EPY/validation',
        'AGE_ANALYSIS/deep_dive/07_FIVE_LOCUS_BETA_VALIDATION/tables',
        'AGE_ANALYSIS/deep_dive/07_FIVE_LOCUS_BETA_VALIDATION/validation',
        'AGE_ANALYSIS/deep_dive/09_EXTERNAL_VALIDATION/tables',
        'AGE_ANALYSIS/deep_dive/09_EXTERNAL_VALIDATION/validation',
        'AGE_ANALYSIS/deep_dive/FINAL_AGE_STORY/figures','figure-preview',
        'public-cohort-cache',
    ]:
        (workspace/folder).mkdir(parents=True,exist_ok=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace',type=Path,required=True)
    parser.add_argument('--wgs-processed',type=Path)
    parser.add_argument('--external-cache',type=Path)
    parser.add_argument('stage',choices=['init']+list(STAGES))
    parser.add_argument('arguments',nargs=argparse.REMAINDER)
    a=parser.parse_args();workspace=a.workspace.resolve()
    if workspace==BASE.parent or BASE.parent in workspace.parents:
        parser.error('Use a separate data workspace outside the repository.')
    if a.stage != 'init' and (workspace/'PUBLIC_INPUTS.json').exists():
        parser.error('Privacy-transformed public tables cannot be used for exact clinical-age analyses or manuscript figure regeneration. Use release/public_source_data.py for public outputs; see release/README.md.')
    initialize(workspace)
    if a.stage=='init':
        print('Source installed; no data analysis run. Supply processed inputs described in age/README.md.')
        return
    env=os.environ.copy()
    env.update(AGE_WORKSPACE=str(workspace),
               AGE_WGS_PROCESSED=str((a.wgs_processed or workspace/'inputs/wgs_processed').resolve()),
               AGE_EXTERNAL_CACHE=str((a.external_cache or workspace/'public-cohort-cache').resolve()),
               AGE_FIGURE_TEMP=str(workspace/'figure-preview'),
               PYTHONDONTWRITEBYTECODE='1')
    rel,args=STAGES[a.stage]
    extra=a.arguments[1:] if a.arguments[:1]==['--'] else a.arguments
    subprocess.run([sys.executable,str(workspace/rel)]+args+extra,env=env,cwd=workspace,check=True)


if __name__=='__main__':
    main()
