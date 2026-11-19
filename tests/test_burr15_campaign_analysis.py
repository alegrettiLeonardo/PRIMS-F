import unittest

from analyze_burr15_campaign import front_candidates, local_phase_transport, mass_balance, unique_centres


def profile(*points):
    return [{"s_m": position, "alpha": alpha} for position, alpha in points]


class CampaignObservationsTest(unittest.TestCase):
    def test_interpolates_a_translating_void_fraction_pulse(self):
        first = profile((0, .25), (1, .75), (2, .75), (3, .25))
        second = profile((1, .25), (2, .75), (3, .75), (4, .25))
        self.assertEqual(list(front_candidates(first, 5, [.5])),
                         [(.5, .5, "rising_alpha"), (.5, 2.5, "falling_alpha")])
        self.assertEqual(list(front_candidates(second, 5, [.5])),
                         [(.5, 1.5, "rising_alpha"), (.5, 3.5, "falling_alpha")])

    def test_ignores_a_matched_junction_holdup_jump(self):
        state = profile((0, .25), (1, .25), (2, .75), (3, .75))
        self.assertEqual(list(front_candidates(state, 1.5, [.5])), [])

    def test_transport_retains_velocity_sign_and_void_weighting(self):
        self.assertEqual(local_phase_transport(.25, 2.0, -4.0, 1000., 2., .01),
                         (1.5, -1., 15., -.02))

    def test_conflicting_terminal_duplicate_fails(self):
        base = {"time_s": "1.5", "s_m": "0.5", "alpha": ".25"}
        self.assertEqual(unique_centres([base, dict(base)]), [base])
        with self.assertRaisesRegex(ValueError, "Conflicting field records"):
            unique_centres([base, dict(base, alpha=".75")])

    def test_integral_mass_uses_cumulative_flux_differences(self):
        data = {"time_s": [0., 1., 2.]}
        for phase in ("g", "l"):
            data[f"M{phase}_kg"] = [10., 10.5, 11.]
            data[f"M{phase}_in_integral_kg"] = [0., 1.5, 3.]
            data[f"M{phase}_out_integral_kg"] = [0., 1., 2.]
            data[f"balance_{phase}_kg"] = [0., 0., 0.]
        report = mass_balance(data, [.25, 1.25])
        self.assertEqual(report["g"]["cycles"][0]["error_kg"], 0.)
        self.assertEqual(report["l"]["relative_to_initial_inventory"], 0.)


if __name__ == "__main__":
    unittest.main()
