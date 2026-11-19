"""Fail-closed evidence and qualification guards, independent of solver runs."""
import copy
import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
from summarize_static_bc_sensitivity import validate_terminal_identity, cycle_metrics_qualified

class EvidenceContract(unittest.TestCase):
    def setUp(self):
        self.manifest = dict(status='COMPLETED', returncode=0, run_status=dict(
            status=0, termination_reason='COMPLETED', reached_time_s=3600., requested_time_s=3600.))
        self.cross = dict(provisional_running_snapshot=False, runs=[dict(
            manifest=copy.deepcopy(self.manifest), completed=True)])

    def test_terminal_and_failed_are_distinct(self):
        validate_terminal_identity(self.manifest, self.cross)
        failed=copy.deepcopy(self.manifest);failed['status']='FAILED';failed['returncode']=1
        failed['run_status'].update(status=13,termination_reason='SOLVER_FAILURE',reached_time_s=4.)
        validate_terminal_identity(failed, dict(runs=[dict(manifest=failed,completed=False)]))

    def test_provisional_cannot_be_terminal_evidence(self):
        self.cross['provisional_running_snapshot']=True
        with self.assertRaises(ValueError):validate_terminal_identity(self.manifest,self.cross)

    def test_copied_manifest_cannot_qualify(self):
        self.cross['runs'][0]['manifest']['run_status']['requested_time_s']=7200.
        with self.assertRaises(ValueError):validate_terminal_identity(self.manifest,self.cross)

    def test_completion_roundoff_matches_sealed_executor(self):
        self.manifest['run_status']['reached_time_s']=3599.9999999999995
        self.cross['runs'][0]['manifest']=copy.deepcopy(self.manifest)
        validate_terminal_identity(self.manifest,self.cross)
        self.manifest['run_status'].update(reached_time_s=float('inf'),requested_time_s=float('inf'))
        self.cross['runs'][0]['manifest']=copy.deepcopy(self.manifest)
        with self.assertRaises(ValueError):validate_terminal_identity(self.manifest,self.cross)

    def test_false_completion_cannot_qualify(self):
        self.manifest['run_status']['reached_time_s']=3599.
        self.cross['runs'][0]['manifest']=copy.deepcopy(self.manifest)
        with self.assertRaises(ValueError):validate_terminal_identity(self.manifest,self.cross)

    def test_inconsistent_summary_cannot_qualify(self):
        self.cross['runs'][0]['completed']=False
        with self.assertRaises(ValueError):validate_terminal_identity(self.manifest,self.cross)

    def test_qualification_requires_all_six_balances_and_accommodation(self):
        balances={k:dict(passes_inventory_1e6=True,passes_processed_1e6=True)
                  for k in ('l_system','l_pipe','l_riser','g_system','g_pipe','g_riser')}
        result=dict(completed=True,cycles_available=True,accommodated=True,
                    arithmetic_contract_passes=True,phase_region_balances=balances,
                    cycle_balance_audit=dict(intervals=10,phase_regions=copy.deepcopy(balances)))
        self.assertTrue(cycle_metrics_qualified(result))
        for flag in ('completed','cycles_available','accommodated','arithmetic_contract_passes'):
            changed=copy.deepcopy(result);changed[flag]=False
            self.assertFalse(cycle_metrics_qualified(changed))
        changed=copy.deepcopy(result);changed['phase_region_balances']['g_riser']['passes_inventory_1e6']=False
        self.assertFalse(cycle_metrics_qualified(changed))
        changed=copy.deepcopy(result);changed['cycle_balance_audit']['phase_regions'].pop('g_pipe')
        self.assertFalse(cycle_metrics_qualified(changed))

if __name__=='__main__':unittest.main()
