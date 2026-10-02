"""Display privacy-transformed public summaries; never refit study analyses."""
from pathlib import Path
import argparse
import csv
import re

ROOT = Path(__file__).resolve().parent
TABLES = ROOT / 'processed/AGE_ANALYSIS/deep_dive/FINAL_AGE_STORY'
COHORT_COLUMNS = ['Sample ID', 'Diagnosis', 'Age at sampling (privacy-minimized)',
                  'Sex', 'Sampling status', 'DNAm age (years)', 'WGS-derived purity',
                  'WGS purity status']
PHENOTYPE_COLUMNS = ['sample_id', 'tumour_type', 'chronological_age_public', 'Age_BS',
    'Age_BS_CI_low', 'Age_BS_CI_high', 'clock_CpGs_total', 'clock_CpGs_observed',
    'clock_CpG_coverage_fraction', 'unknown_beta_lower_years', 'unknown_beta_upper_years',
    'median_clock_CpG_depth', 'mean_clock_CpG_depth', 'Age_M', 'Age_M_CI_low', 'Age_M_CI_high',
    'hydroxymethylation_clock_effect', 'hydroxymethylation_clock_effect_CI_low',
    'hydroxymethylation_clock_effect_CI_high', 'absolute_measured_5mC_signal_fraction',
    'primary_purity', 'primary_ploidy', 'primary_purity_status',
    'relative_abs_log2_gt_0p2_fraction', 'purity_at_search_boundary',
    'five_reference_contrast_score', 'five_reference_contrast_years',
    'technical_status', 'scientific_fit_class', 'outcome_class', 'loh_status']
SURVIVAL_COLUMNS = ['Analysis', 'n', 'Deaths', 'Predictor', 'Scaling', 'HR',
                    '95% CI lower', '95% CI upper', 'Likelihood-ratio P', 'Notes']
ALLOWED = {f'GBM-{i:02d}' for i in range(1, 16)} | {f'MEN-{i:02d}' for i in range(1, 6)}


def read(path, columns):
    with path.open() as f:
        reader = csv.DictReader(f, delimiter='\t')
        if reader.fieldnames != columns:
            raise ValueError('Only the privacy-minimized public schema is accepted: ' + path.name)
        return list(reader)


def validate_local(rows, id_column, age_column):
    if len(rows) != 20 or {r[id_column] for r in rows} != ALLOWED:
        raise ValueError('Expected the fixed 20 study pseudonyms.')
    for row in rows:
        value = row[age_column]
        if value != '<18' and not (re.fullmatch(r'\d+', value) and 18 <= int(value) <= 120):
            raise ValueError('Public age must be adult completed years or <18.')


def load_public():
    cohort = read(TABLES/'supplementary/Supplementary_Table_1_cohort.tsv', COHORT_COLUMNS)
    phenotype = read(TABLES/'tables/patient_phenotype_WGS_chemistry.tsv', PHENOTYPE_COLUMNS)
    survival = read(TABLES/'supplementary/Supplementary_Table_2_survival_aggregate.tsv', SURVIVAL_COLUMNS)
    validate_local(cohort, 'Sample ID', 'Age at sampling (privacy-minimized)')
    validate_local(phenotype, 'sample_id', 'chronological_age_public')
    if len(survival) != 3 or any(re.search(r'\b(?:GBM|MEN)-\d+', v) for row in survival for v in row.values()):
        raise ValueError('Survival output must contain only three aggregate model summaries.')
    return cohort, phenotype, survival


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    cohort, phenotype, survival = load_public()
    args.output.mkdir(parents=True, exist_ok=True)
    from openpyxl import Workbook
    from openpyxl.styles import Font
    workbook = Workbook()
    note = workbook.active
    note.title = 'Privacy notice'
    note.append(['Privacy-transformed public summaries; not exact manuscript source data.'])
    note.append(['Adult ages are completed integer years; minors are <18.'])
    note.append(['Exact age and patient survival inputs require ethics approval and an executed DTA.'])
    note.append(['Do not rerun clinical-age analyses using these transformed ages.'])
    for name, rows in [('Public cohort', cohort), ('Molecular values', phenotype), ('Aggregate survival', survival)]:
        sheet = workbook.create_sheet(name)
        sheet.append(list(rows[0]))
        for row in rows:
            sheet.append(list(row.values()))  # Preserve category strings and supplied numeric text exactly.
        sheet.freeze_panes = 'A2'
        for cell in sheet[1]:
            cell.font = Font(bold=True)
    workbook.properties.creator = 'Study authors'
    workbook.properties.title = 'Privacy-transformed public source summaries'
    workbook.save(args.output/'privacy_transformed_public_summaries.xlsx')

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(10, 7))
    for i, row in enumerate(phenotype):
        age = float(row['Age_BS'])
        axes[0].errorbar(age, i, xerr=[[age-float(row['Age_BS_CI_low'])],
                                      [float(row['Age_BS_CI_high'])-age]], fmt='o', color='#236987')
        axes[1].barh(i, 100*float(row['absolute_measured_5mC_signal_fraction']), color='#236987')
    for ax in axes:
        ax.set_yticks(range(20))
        ax.set_yticklabels([r['sample_id'] for r in phenotype])
        ax.invert_yaxis()
    axes[0].set_xlabel('Projected DNAm age (model output, years)')
    axes[1].set_xlabel('Absolute measured 5mC clock signal (%)')
    axes[1].set_xlim(0, 100)
    fig.suptitle('Public molecular overview: supplied values only')
    fig.text(.05, .01, 'No chronological-age offsets or survival records; not a replacement manuscript figure.', fontsize=9)
    fig.tight_layout(rect=(0, .04, 1, .95))
    fig.savefig(args.output/'public_molecular_overview.pdf')
    plt.close(fig)
    print('Privacy-labelled public summaries generated; no scientific estimates recalculated.')


if __name__ == '__main__':
    main()
