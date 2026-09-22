import unittest

import numpy as np

from scope.classes_qc import VNM
from scope.vnm_tools import geom_sampling_from_vnm


class TestVNM(unittest.TestCase):

    def setUp(self):
        self.mode = VNM(1, 100.0)
        self.mode.set_mode([1, 2], [1, 8], [1.0, 1.0], [0.0, 0.0], [0.0, 0.0], is_mass_weighted=False)

    def test_weighting_and_unweighting_mutate_the_authoritative_mode(self):
        cartesian_mode = self.mode.mode.copy()
        expected_weighted = cartesian_mode * np.sqrt(np.asarray(self.mode.masses))[:, np.newaxis]

        returned_weighted = self.mode.mass_weight_mode()

        self.assertTrue(self.mode.is_mass_weighted)
        self.assertTrue(np.allclose(self.mode.mode, expected_weighted))
        self.assertIs(returned_weighted, self.mode.mode)
        self.assertNotIn("mode_format2", vars(self.mode))
        self.assertTrue(np.shares_memory(self.mode.mode_format2, self.mode.mode))

        returned_cartesian = self.mode.unweight_mode()

        self.assertFalse(self.mode.is_mass_weighted)
        self.assertTrue(np.allclose(self.mode.mode, cartesian_mode))
        self.assertIs(returned_cartesian, self.mode.mode)

    def test_atomic_participation_is_independent_of_mode_representation(self):
        cartesian_participation = self.mode.get_atomic_participation()
        self.mode.mass_weight_mode()
        weighted_participation = self.mode.get_atomic_participation()

        self.assertTrue(np.allclose(cartesian_participation, weighted_participation))
        self.assertAlmostEqual(float(np.sum(weighted_participation)), 1.0)
        self.assertGreater(weighted_participation[1], weighted_participation[0])

    def test_zero_norm_mode_has_no_atomic_participation(self):
        mode = VNM(1, 100.0)
        mode.set_mode([1], [1], [0.0], [0.0], [0.0], is_mass_weighted=True)

        with self.assertRaisesRegex(ValueError, "invalid norm"):
            mode.get_atomic_participation()

    def test_overlap_requires_the_same_weighting_convention(self):
        other = VNM(2, 100.0)
        other.set_mode([1, 2], [1, 8], [1.0, 1.0], [0.0, 0.0], [0.0, 0.0], is_mass_weighted=True)

        with self.assertRaisesRegex(ValueError, "same mass-weighting convention"):
            self.mode.overlap(other)

    def test_geometry_sampling_requires_mass_weighted_modes(self):
        with self.assertRaisesRegex(ValueError, "must be mass weighted"):
            geom_sampling_from_vnm(["H", "O"], np.zeros((2, 3)), [self.mode], n_samples=1, check_adjacencies=False)


if __name__ == "__main__":
    unittest.main()
