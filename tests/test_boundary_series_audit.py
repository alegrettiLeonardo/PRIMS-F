import csv
from pathlib import Path
import tempfile
import unittest
from audit_boundary_series import audit


class BoundaryCoverageTests(unittest.TestCase):
    def test_empty_or_partial_coverage_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/"empty.csv"
            path.write_text("time_s,side,kind,status\n")
            with self.assertRaisesRegex(ValueError,"Incomplete"):
                audit(path,1)

    def test_domain_failures_are_counted_and_never_classified(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/"states.csv"
            with path.open("w",newline="") as stream:
                writer=csv.writer(stream)
                writer.writerow(["time_s","side","kind","status"])
                for side in (1,2):
                    for kind in (1,2):
                        writer.writerow([0,side,kind,7])
            result=audit(path,1)
            self.assertEqual(result["1_1"]["status_counts"],{7:1})
            self.assertEqual(result["2_2"]["incoming_counts"],{})
            self.assertIsNone(result["2_2"]["min_constraint_singular_value"])


if __name__=="__main__":
    unittest.main()
