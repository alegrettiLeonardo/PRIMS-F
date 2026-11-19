import csv
import importlib.util
import json
import math
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "verify_burr15_observations.py"
if not SCRIPT.exists():
    SCRIPT = Path(__file__).with_name("verify_burr15_observations.py")
SPEC = importlib.util.spec_from_file_location("burr15_verify", SCRIPT)
MOD = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MOD)


def field(t=0, ug=0.2, s=1.0, **changes):
    result = dict(time_s=t, s_m=s, theta_rad=0, alpha=0.5,
                  holdup_liquid=0.5, p_Pa=200000, ul_m_s=0.5, ug_m_s=ug,
                  rhol_kg_m3=998, rhog_kg_m3=2.3)
    result.update(changes)
    return result


def observed(f, **changes):
    # Closed half-full geometry, not the verifier's inversion or reference().
    threshold2 = ((998-2.3)/998)*9.81*0.05254/16
    ug2 = f["ug_m_s"]**2
    result = dict(time_s=f["time_s"], s_m=f["s_m"], closure_id=1, status=0,
                  hL_over_D=0.5, ug2_m2_s2=ug2, ug2_critical_m2_s2=threshold2,
                  margin_m2_s2=ug2-threshold2, unstable=int(ug2 >= threshold2))
    result.update(changes)
    return result


class BurrReferenceTests(unittest.TestCase):
    def test_independent_height_and_numeric_references(self):
        for alpha, height, critical in [
            (0.5, 0.5, 0.17927450395772561488),
            (0.804498890522114679, 0.25, 0.41776604691125290893),
            (0.195501109477885321, 0.75, 0.04854100642803791632)]:
            result = MOD.reference(field(alpha=alpha), 0.05254, 9.81)
            self.assertEqual(result["status"], 0)
            self.assertAlmostEqual(result["hL_over_D"], height, places=14)
            self.assertAlmostEqual(result["critical_ug_m_s"], critical, places=13)

    def test_height_tails_remain_interior(self):
        for alpha in [1e-100, 1e-14, 1e-6, 1-1e-6, 1-1e-14]:
            h, room, angle = MOD.height_geometry(alpha)
            self.assertGreater(room, 0)
            self.assertGreater(h, 0)
            self.assertGreater(angle, 0)
            cap = room if alpha < 0.5 else h
            self.assertTrue(math.isclose(MOD.cap_fraction(cap), min(alpha,1-alpha), rel_tol=2e-12))

    def test_domain_and_sign_are_not_squared_away(self):
        cases = [({"ug_m_s":-0.2},3), ({"ul_m_s":-0.2},3),
                 ({"theta_rad":0.01},4), ({"theta_rad":-math.pi/2},4),
                 ({"alpha":0},2), ({"alpha":1},2), ({"alpha":-0.1},1),
                 ({"rhog_kg_m3":998},1), ({"ug_m_s":math.nan},1)]
        for changes, status in cases:
            with self.subTest(changes=changes):
                self.assertEqual(MOD.reference(field(**changes),.05254,9.81)["status"],status)
        self.assertEqual(MOD.reference(field(),.05254,9.81,area=1)["status"],6)
        self.assertEqual(MOD.reference(field(),.05254,9.81,closure_id=3)["status"],5)
        self.assertEqual(MOD.reference(field(state_status=14),.05254,9.81)["status"],5)

    def test_density_and_diameter_scales(self):
        base = MOD.reference(field(),.05254,9.81)
        scaled = MOD.reference(field(rhol_kg_m3=9980,rhog_kg_m3=23),4*.05254,9.81)
        self.assertAlmostEqual(scaled["critical_ug_m_s"],2*base["critical_ug_m_s"],places=14)
        angled = MOD.reference(field(theta_rad=-2*math.pi/180),.05254,9.81)
        self.assertAlmostEqual(angled["critical_ug_m_s"],.17921989104729318697,places=14)


class BurrObservationTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def verify(self, fields, rows=None, **kwargs):
        rows = [observed(f) for f in fields] if rows is None else rows
        for name, values in [("fields.csv",fields),("burr15.csv",rows)]:
            with (self.root/name).open("w",newline="") as out:
                writer = csv.DictWriter(out,fieldnames=values[0].keys())
                writer.writeheader()
                writer.writerows(values)
        return MOD.verify(self.root/"fields.csv",self.root/"burr15.csv",.05254,
                          expected_cells=kwargs.pop("expected_cells",1),**kwargs)

    def test_loss_intervals_and_supplied_cycles(self):
        rows = [field(0,.1),field(.25,.2),field(.5,.3),field(.75,.1)]
        result = self.verify(rows,cycle_edges=[0,.5,1])
        self.assertTrue(result["consistent"])
        self.assertEqual(result["first_observed_loss"]["time_s"],.25)
        self.assertEqual(result["reference_bounds"]["ug_m_s"],{"min":.1,"max":.3})
        self.assertLess(result["reference_bounds"]["margin_m2_s2"]["min"],0)
        self.assertGreater(result["reference_bounds"]["margin_m2_s2"]["max"],0)
        interval, = result["observed_loss_intervals"]
        self.assertEqual(interval["onset_bracket_s"],[0,.25])
        self.assertEqual(interval["end_bracket_s"],[.5,.75])
        self.assertFalse(interval["left_censored"])
        self.assertFalse(interval["right_censored"])
        self.assertEqual([c["frames_with_observed_loss"] for c in result["supplied_cycles"]],[1,1])

    def test_identical_final_duplicates_are_not_extra_samples(self):
        result = self.verify([field(0),field(.25),field(.25)])
        self.assertEqual(result["counts"]["identical_duplicates"],2)
        self.assertEqual(result["counts"]["samples"],2)
        self.assertEqual(result["observed_loss_intervals"][0]["samples"],2)

    def test_conflicting_duplicate_refused(self):
        with self.assertRaisesRegex(ValueError,"Conflicting duplicate"):
            self.verify([field(0),field(0,.4)])

    def test_missing_burr_frame_breaks_continuity(self):
        fields = [field(0),field(.25),field(.5)]
        result = self.verify(fields,[observed(fields[0]),observed(fields[2])])
        self.assertFalse(result["consistent"])
        self.assertEqual(result["counts"]["field_only_frames"],1)
        intervals = result["observed_loss_intervals"]
        self.assertEqual(len(intervals),2)
        self.assertTrue(intervals[0]["right_censored"])
        self.assertTrue(intervals[1]["left_censored"])

    def test_common_sampling_gap_breaks_continuity(self):
        result = self.verify([field(0),field(1)],max_gap=.5)
        self.assertEqual(result["counts"]["time_gaps"],1)
        self.assertEqual(len(result["observed_loss_intervals"]),2)

    def test_expected_grid_detects_truncated_first_frame(self):
        result = self.verify([field(0)],expected_cells=2)
        self.assertFalse(result["consistent"])
        self.assertEqual(result["counts"]["incomplete_frames"],1)
        self.assertEqual(result["observed_loss_intervals"],[])

    def test_domain_failure_does_not_become_stable_observation(self):
        f0,f1,f2 = field(0),field(.25,-.2),field(.5)
        rows=[observed(f0),observed(f1,status=3,unstable=0),observed(f2)]
        result=self.verify([f0,f1,f2],rows)
        self.assertTrue(result["consistent"])
        self.assertEqual(result["domain_counts"]["DIRECTION"],1)
        self.assertEqual(len(result["observed_loss_intervals"]),2)

    def test_discriminates_density_perimeter_and_si_velocity_errors(self):
        f=field()
        for changes,key in [({"ug2_critical_m2_s2":3.734391812059664**2},"ug2_critical_m2_s2"),
                            ({"ug2_critical_m2_s2":4.680366052329626**2},"ug2_critical_m2_s2"),
                            ({"ug2_m2_s2":400},"ug2_m2_s2"),
                            ({"margin_m2_s2":math.nan},"margin_m2_s2"),
                            ({"status":3},"status"),({"unstable":0},"unstable")]:
            result=self.verify([f],[observed(f,**changes)])
            self.assertFalse(result["consistent"])
            self.assertEqual(result["mismatches"][key],1)
            json.dumps(result,allow_nan=False)

    def test_window_does_not_invent_unobserved_onset(self):
        result=self.verify([field(0,.1),field(.25),field(.5),field(.75,.1)],start=.25,end=.5)
        interval,=result["observed_loss_intervals"]
        self.assertTrue(interval["left_censored"])
        self.assertTrue(interval["right_censored"])
        self.assertEqual(result["counts"]["frames"],2)

    def test_cycle_summary_requires_matching_run_and_provenance(self):
        provenance={}
        for name in ["monitor.csv","manifest.json","run_metadata.txt","run_status.txt"]:
            (self.root/name).write_text(name)
            provenance[name]=MOD.file_sha256(self.root/name)
        run=dict(run=self.root.name,cycles_available=True,cycle_edges_s=[0,1],
                 conservation={"provenance":provenance})
        p=self.root/"cycles.json"
        p.write_text(json.dumps({"runs":[run]}))
        self.assertEqual(MOD.read_cycle_edges(p,self.root),[0,1])
        (self.root/"run_status.txt").write_text("other")
        with self.assertRaisesRegex(ValueError,"provenance"):
            MOD.read_cycle_edges(p,self.root)
        run["run"]="another_case"
        p.write_text(json.dumps({"runs":[run]}))
        with self.assertRaisesRegex(ValueError,"another run"):
            MOD.read_cycle_edges(p,self.root)

    def test_input_digests_are_reproducible(self):
        result=self.verify([field()])
        for item in result["inputs"]:
            self.assertEqual(item["sha256"],MOD.file_sha256(item["path"]))


if __name__ == "__main__":
    unittest.main()
