"""Privacy boundaries for the public layer; no private input access."""
from pathlib import Path
import importlib.util
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('public_source_data', ROOT/'release/public_source_data.py')
public = importlib.util.module_from_spec(spec)
spec.loader.exec_module(public)


class PublicPrivacyTests(unittest.TestCase):
    def test_public_schema_and_cross_table_consistency(self):
        cohort, phenotype, survival = public.load_public()
        ages = {r['sample_id']: r['chronological_age_public'] for r in phenotype}
        self.assertEqual(ages, {r['Sample ID']: r['Age at sampling (privacy-minimized)'] for r in cohort})
        self.assertEqual(sum(a == '<18' for a in ages.values()), 1)
        self.assertTrue(all('DNAm_age_acceleration' not in r and 'remaining_DeltaAge' not in r for r in phenotype))
        self.assertTrue(all('Vital status' not in r and 'Post-sampling follow-up (days)' not in r for r in cohort))
        self.assertEqual([(r['HR'], r['95% CI lower'], r['95% CI upper'], r['Likelihood-ratio P']) for r in survival],
                         [('0.745', '0.342', '1.621', '0.414'), ('0.733', '0.343', '1.567', '0.375'), ('0.733', '0.343', '1.567', '0.375')])

    def test_rejects_precise_age_and_unapproved_local_keys(self):
        cohort, _, _ = public.load_public()
        for value in ['18.25', '17', '17.99']:
            changed = [dict(r) for r in cohort]
            changed[0]['Age at sampling (privacy-minimized)'] = value
            with self.assertRaises(ValueError):
                public.validate_local(changed, 'Sample ID', 'Age at sampling (privacy-minimized)')
        changed = [dict(r) for r in cohort]
        changed[0]['Sample ID'] = 'UNAPPROVED'
        with self.assertRaises(ValueError):
            public.validate_local(changed, 'Sample ID', 'Age at sampling (privacy-minimized)')

    def test_launcher_rejects_transformed_workspace_before_analysis(self):
        with tempfile.TemporaryDirectory() as folder:
            workspace = Path(folder)
            (workspace/'PUBLIC_INPUTS.json').write_text('{"privacy_transformed": true}')
            result = subprocess.run([sys.executable, str(ROOT/'age/run.py'), '--workspace', folder, 'main-figures'], capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('Privacy-transformed public tables cannot be used', result.stderr)
            self.assertFalse((workspace/'AGE_ANALYSIS').exists())


if __name__ == '__main__':
    unittest.main()
