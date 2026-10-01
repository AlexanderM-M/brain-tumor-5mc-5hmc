"""Synthetic validation only; never loads patient inputs or runs study analyses."""
from pathlib import Path
import importlib.util
import sys
import math
import numpy as np
import unittest
import tempfile
from unittest.mock import patch
import os

AGE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(AGE))
from survival_core import likelihood,fit,validate_math
from survival import load
from run import initialize


def module(name,relative):
    p=AGE/'workflow'/relative
    spec=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


class ReleaseTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.workspace=Path(self.temp.name)
        self.env=patch.dict(os.environ, {}, clear=False)
        self.env.start()
        self.addCleanup(self.env.stop)

    def test_cox_risk_likelihood_and_independent_optimum(self):
        assert np.allclose(likelihood(0,np.array([1.,2.]),np.array([1,0]),np.array([1.,0.])),[-math.log(2),.5,.25])
        t=np.array([1.,2.,3.,4.,5.]);e=np.array([1,1,1,0,0]);x=np.array([-1.,1.,0.,.4,-.3])
        r=fit(t,e,x);assert r['status']=='ESTIMATED_SPARSE_EVENTS';validate_math(t,e,x,r)
        assert fit(t,np.zeros(5),x)['status']=='FAILED_ZERO_EVENTS'


    def test_clock_transform_and_reference_substitution(self):
        tmp_path=self.workspace
        os.environ['AGE_WORKSPACE']=str(tmp_path)
        c=module('clock_common_synthetic','AGE_ANALYSIS/next_generation/scripts/common.py')
        assert np.allclose(c.inv(np.array([-1.,0.,1.])),[21/math.e-1,20,41])
        folder=AGE/'workflow/AGE_ANALYSIS/deep_dive/06_independent_preparation/scripts'
        sys.path.insert(0,str(folder))
        m=module('five_null_synthetic','AGE_ANALYSIS/deep_dive/06_independent_preparation/scripts/five_null.py')
        beta=np.array([[.4,.8],[.1,.2]]);ref=np.array([.5,.5]);weights=np.array([1.,-.2])
        r=m.replacement_effect(beta,ref,weights,[0]);expected=c.inv(.696+beta@weights)-c.inv(.696+np.column_stack([np.full(2,.5),beta[:,1]])@weights)
        assert np.allclose(r['contribution_years'],expected)


    def test_four_cpg_ordered_chemistry(self):
        m=module('epiallele_synthetic','AGE_ANALYSIS/deep_dive/06_PDR_EPY/scripts/epiallele_metrics.py')
        counts=np.zeros((1,81),dtype=np.uint32);counts[0,40]=10;counts[0,80]=10
        r=m.metrics(counts)
        assert r['PDR'][0]==0 and r['binary_entropy'][0]==0
        assert r['three_state_entropy'][0]>0 and r['hidden_information_bits'][0]>0


    def test_survival_input_rejects_extra_clinical_field(self):
        tmp_path=self.workspace
        p=tmp_path/'synthetic.tsv'
        p.write_text('sample_id\tdiagnosis\tDNAm_age_acceleration\tpost_sampling_survival_days\tevent\tdeath_date\nGBM-01\tGlioblastoma\t12\t100\t1\tNA\n')
        with self.assertRaises(ValueError):load(p)


    def test_workspace_install_is_non_destructive(self):
        tmp_path=self.workspace
        initialize(tmp_path)
        p=tmp_path/'AGE_ANALYSIS/next_generation/scripts/common.py';assert p.exists()
        initialize(tmp_path)
        p.write_text('# user modification\n')
        with self.assertRaises(RuntimeError):initialize(tmp_path)
