"""Prevent accidental mixing of historically distinct numerical campaigns."""
import importlib.util,json,tempfile,unittest
from pathlib import Path
P=Path(__file__).resolve().parents[1]/'validation/diagnostics/p0p1_20261119/compose_accepted_history.py'
spec=importlib.util.spec_from_file_location('history_versions',P);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
class HistoryVersions(unittest.TestCase):
 def folders(self,root,versions):
  result=[]
  for i,v in enumerate(versions):
   folder=root/str(i);folder.mkdir();(folder/'manifest.json').write_text(json.dumps(dict(operator=v)));result.append(dict(run=str(folder)))
  return result
 def test_historical_baselines_cannot_be_spliced(self):
  with tempfile.TemporaryDirectory()as folder:
   with self.assertRaises(ValueError):m.validate_operators(self.folders(Path(folder),['05ad799','490936a']))
 def test_short_and_full_references_to_same_sha_are_accepted(self):
  with tempfile.TemporaryDirectory()as folder:
   self.assertEqual(m.validate_operators(self.folders(Path(folder),['18aa7f7','18aa7f79f7b16442fa19ebbd9ef8819b0eff187a'])),'18aa7f7')
