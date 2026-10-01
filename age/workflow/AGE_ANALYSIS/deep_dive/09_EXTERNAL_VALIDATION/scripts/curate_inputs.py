"""Explicit sample linkage and eligibility; no biological model fitting."""
import os
from pathlib import Path
from pathlib import Path
import numpy as np
import pandas as pd
from read_archive import LOCI, PROJECT


def load(cache=Path(os.environ['AGE_EXTERNAL_CACHE'])):
    def meta(stem):
        return pd.read_csv(cache/(stem+'.metadata.tsv'), sep='\t', index_col=0, dtype=str)
    def mat(stem):
        d = pd.read_csv(cache/(stem+'.targets.tsv'), sep='\t', index_col=0)
        return d.apply(pd.to_numeric, errors='coerce')
    audit, cohorts = {}, {}
    # Normal brain: link SampleN titles to SampleN_Beta/ SampleN_detectP.
    nm = meta('GSE74193_series_matrix')
    normal_raw = mat('GSE74193_GEO_procData')
    def key(prefix):
        return next(k for k in nm if k.startswith(prefix))
    nm['age'] = pd.to_numeric(nm[key('age (')], errors='coerce')
    nm['donor'] = nm[key('brnum (')]
    nm['male'] = nm[key('sex (')].str.lower().map({'m':1., 'male':1., 'f':0., 'female':0.})
    nm['neuron'] = pd.to_numeric(nm[key('comp_neun_pos')], errors='coerce')
    nm['plate'] = nm[key('plate (')]
    nm['matrix_id'] = nm.title.str.extract(r'^(Sample\d+)_', expand=False)
    assert nm.matrix_id.notna().all() and nm.matrix_id.is_unique
    keep = ((nm['group'] == 'Control') & (nm['age'] >= 18)
            & (nm[key('bestqc (')].str.upper() == 'TRUE')
            & (nm[key('dropsample (')].str.upper() == 'FALSE'))
    nmeta = nm.loc[keep].copy()
    assert nmeta.donor.is_unique
    nval = normal_raw[[s+'_Beta' for s in nmeta.matrix_id]].T
    ndet = normal_raw[[s+'_detectP' for s in nmeta.matrix_id]].T
    nval.index = ndet.index = nmeta.index
    nval = nval.mask(ndet > .01).mask(ndet.isna())
    audit['GSE74193'] = {'arrays': len(nm), 'eligible_adult_controls': len(nmeta),
                       'age_range': [float(nmeta.age.min()), float(nmeta.age.max())],
                       'group_counts': nm['group'].value_counts().to_dict(),
                       'dropsample_counts': nm[key('dropsample (')].value_counts().to_dict(),
                       'target_valid_n': nval.notna().sum().to_dict()}
    # Primary resection GBM; explicitly exclude named recurrences and spheres.
    gm = meta('GSE60274_series_matrix')
    gv = mat('GSE60274_series_matrix').T
    assert set(gv.index) == set(gm.index)
    gm['age'] = pd.to_numeric(gm.age, errors='coerce')
    gm['donor'] = gm.title.str.extract(r'(GBM-\d+|NTB-\d+)', expand=False)
    gkeep = gm.title.str.startswith('genomic DNA from Surgical resection GBM-')
    gctrl = gm['disease state'].eq('Non-tumor brain')
    assert gm.loc[gkeep | gctrl, 'donor'].is_unique
    cohorts['GSE60274'] = dict(meta=gm.loc[gkeep], beta=gv.loc[gkeep],
                              control_meta=gm.loc[gctrl], control_beta=gv.loc[gctrl],
                              notes='64 primary GBM; historical diagnoses, IDH/G-CIMP not supplied; five epilepsy-surgery non-tumour controls; detection P unavailable in series matrix.')
    audit['GSE60274'] = {'arrays': len(gm), 'primary_GBM': int(gkeep.sum()),
                       'non_tumour_controls': int(gctrl.sum()),
                       'excluded_recurrences': int(gm.title.str.contains('recurrent', case=False).sum()),
                       'excluded_spheres': int(gm.title.str.contains('spheres', case=False).sum()),
                       'primary_age_available': int(gm.loc[gkeep, 'age'].notna().sum()),
                       'control_age_available': int(gm.loc[gctrl, 'age'].notna().sum())}
    # Nordic trial: deposited patient identifiers explicitly link beta columns.
    km = meta('GSE195684_series_matrix')
    kv = mat('GSE195684_Matrix_Processed_BADY_GBMnordic')
    km['age'] = pd.to_numeric(km.age, errors='coerce')
    km['donor'] = km.patient
    assert km.patient.is_unique
    cols = ['GBM-'+p for p in km.patient]
    kbeta = kv[cols].T
    kdet = kv[[s+' Detection Pval' for s in cols]].T
    kbeta.index = kdet.index = km.index
    kbeta = kbeta.mask(kdet > .01).mask(kdet.isna())
    cohorts['GSE195684'] = dict(meta=km, beta=kbeta, control_meta=None, control_beta=None,
                               notes='Nordic trial FFPE GBM; deposited metadata do not include individual IDH/G-CIMP status; per-probe detection P <=0.01.')
    audit['GSE195684'] = {'GBM': len(km), 'age_available': int(km.age.notna().sum()),
                         'target_valid_n': kbeta[LOCI].notna().sum().to_dict(),
                         'age_range': [float(km.age.min()), float(km.age.max())]}
    # TCGA sample IDs match exactly, not by approximate sample-name similarity.
    tm = pd.read_csv(cache/'TCGA_clinical.tsv', sep='\t', index_col=0, dtype=str)
    tv = mat('HumanMethylation450').T
    assert set(tv.index) <= set(tm.index) and tm.index.is_unique
    tm = tm.loc[tv.index].copy()
    tm['age'] = pd.to_numeric(tm.age_at_initial_pathologic_diagnosis, errors='coerce')
    tm['donor'] = tm['_PATIENT']
    tp = tm.sample_type.eq('Primary Tumor')
    tn = tm.sample_type.eq('Solid Tissue Normal')
    assert tm.loc[tp, 'donor'].is_unique and tm.loc[tn, 'donor'].is_unique
    # Exclude primary tumours paired with normal specimens from unpaired contrast.
    paired = set(tm.loc[tp, 'donor']) & set(tm.loc[tn, 'donor'])
    cohorts['TCGA_GBM'] = dict(meta=tm.loc[tp], beta=tv.loc[tp],
                              control_meta=tm.loc[tn], control_beta=tv.loc[tn],
                              exclude_from_unpaired=tm.index[tp & tm.donor.isin(paired)].tolist(),
                              notes='450K primary GBM; age at initial pathological diagnosis; G-CIMP sensitivity reported separately; supplied values rounded; detection P unavailable.')
    noncimp = tp & tm.G_CIMP_STATUS.eq('NON G-CIMP')
    cohorts['TCGA_non_GCIMP_sensitivity'] = dict(meta=tm.loc[noncimp], beta=tv.loc[noncimp],
                              control_meta=None, control_beta=None,
                              notes='Subset of TCGA_GBM; not independent replication; NON G-CIMP annotation does not itself verify IDH wildtype.')
    audit['TCGA_GBM'] = {'arrays': len(tm), 'sample_types': tm.sample_type.value_counts().to_dict(),
                        'primary_GBM': int(tp.sum()), 'age_available': int(tm.loc[tp,'age'].notna().sum()),
                        'G_CIMP_STATUS': tm.loc[tp,'G_CIMP_STATUS'].fillna('UNANNOTATED').value_counts().to_dict(),
                        'normal_tumour_paired_donors': len(paired), 'non_GCIMP_n': int(noncimp.sum())}
    for c in cohorts.values():
        vals = c['beta'].to_numpy()
        assert np.isfinite(vals).sum() and np.nanmin(vals) >= 0 and np.nanmax(vals) <= 1
        assert set(LOCI) <= set(c['beta'].columns)
    assert 0 <= np.nanmin(nval) <= np.nanmax(nval) <= 1
    patient = pd.read_csv(PROJECT/'AGE_ANALYSIS/deep_dive/FINAL_AGE_STORY/tables/patient_phenotype_WGS_chemistry.tsv', sep='\t').set_index('sample_id')
    local = pd.read_csv(PROJECT/'AGE_ANALYSIS/deep_dive/05_five_locus_signal/tables/five_locus_modification_values.tsv', sep='\t')
    ont = local.pivot(index='sample_id', columns='CpGmarker', values='total').loc[patient.index]
    return nmeta, nval, cohorts, patient, ont, audit
