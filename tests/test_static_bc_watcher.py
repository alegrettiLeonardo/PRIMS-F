import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import watch_static_bc_sensitivity as watcher


class StaticBoundaryWatcherTest(unittest.TestCase):
    def check_unavailable_cycles(self,completed):
        # Synthetic already-analyzed files: no solver or analyzer is launched.
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);run=root/'run';run.mkdir()
            tools=root/'tools';tools.mkdir();(tools/'manifest.json').write_text('{}')
            plan=root/'plan.json';plan.write_text(json.dumps(dict(jobs=[dict(name='case',axis='ps',static_offset=5000.,output=str(run))])))
            seal=root/'seal.json';seal.write_text(json.dumps(dict(binary_sha256='binary',source_sha256={'source':'hash'})))
            state=dict(status='COMPLETED' if completed else 'FAILED',
                       plan_sha256=hashlib.sha256(plan.read_bytes()).hexdigest(),
                       seal_sha256=hashlib.sha256(seal.read_bytes()).hexdigest(),
                       binary_sha256='binary',source_sha256={'source':'hash'},
                       binary_unchanged=True,sources_unchanged=True)
            (run/'manifest.json').write_text(json.dumps(state))
            output=root/'analysis/case';output.mkdir(parents=True)
            (output/'accepted_conservation.json').write_text(json.dumps(dict(ledger=dict(arithmetic_contract_passes=True,balances={'l_system':{'passes_inventory_1e6':False}}))))
            (output/'cross_summary.json').write_text(json.dumps(dict(runs=[dict(completed=completed,analysis_available=False,reason='synthetic unavailable window')])) )
            with patch.object(sys,'argv',['watcher','--plan',str(plan),'--analysis-tools',str(tools),'--once']):
                watcher.main()
            result=json.loads((output/'summary.json').read_text())
            self.assertFalse(result['cycles_available'])
            self.assertFalse(result['accommodated'])
            self.assertEqual(result['needs_extension'],completed)
            self.assertEqual(result['analysis_reason'],'synthetic unavailable window')
            self.assertFalse(result['balances']['l_system']['passes_inventory_1e6'])

    def test_early_failure_preserves_unavailable_analysis_and_balance_failure(self):
        self.check_unavailable_cycles(False)

    def test_completed_without_resolved_cycles_is_flagged_for_review_extension(self):
        self.check_unavailable_cycles(True)


if __name__=='__main__':unittest.main()
