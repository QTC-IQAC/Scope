import unittest

import numpy as np

from scope.classes_qc import VNM
from scope.vnm_tools import displace_coords_with_vnm, geom_sampling_from_vnm, map_vnms


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

    def test_nonpermanent_conversions_return_independent_arrays(self):
        cartesian_mode = self.mode.mode.copy()
        weighted_mode = self.mode.mass_weight_mode(permanent=False)

        self.assertFalse(self.mode.is_mass_weighted)
        self.assertTrue(np.array_equal(self.mode.mode, cartesian_mode))
        self.assertFalse(np.shares_memory(weighted_mode, self.mode.mode))

        self.mode.mass_weight_mode()
        stored_weighted_mode = self.mode.mode.copy()
        unweighted_mode = self.mode.unweight_mode(permanent=False)

        self.assertTrue(self.mode.is_mass_weighted)
        self.assertTrue(np.array_equal(self.mode.mode, stored_weighted_mode))
        self.assertTrue(np.allclose(unweighted_mode, cartesian_mode))
        self.assertFalse(np.shares_memory(unweighted_mode, self.mode.mode))

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

    def test_overlap_locally_reconciles_weighting_conventions(self):
        expected_weighted = self.mode.mass_weight_mode(permanent=False)
        other = VNM(2, 100.0)
        other.set_mode([1, 2], [1, 8], expected_weighted[:, 0], expected_weighted[:, 1], expected_weighted[:, 2], is_mass_weighted=True)

        self.assertAlmostEqual(self.mode.overlap(other), 1.0)
        self.assertFalse(self.mode.is_mass_weighted)
        self.assertTrue(other.is_mass_weighted)

    def test_geometry_sampling_locally_weights_modes_without_mutating_them(self):
        original_mode = self.mode.mode.copy()

        geometries, q_coords, energies = geom_sampling_from_vnm(["H", "O"], np.zeros((2, 3)), [self.mode], n_samples=1, check_adjacencies=False)

        self.assertEqual(geometries.shape, (1, 2, 3))
        self.assertEqual(q_coords.shape, (1, 1))
        self.assertEqual(energies.shape, (1,))
        self.assertFalse(self.mode.is_mass_weighted)
        self.assertTrue(np.array_equal(self.mode.mode, original_mode))

    def test_mapping_and_displacement_use_temporary_weighted_modes(self):
        original_mode = self.mode.mode.copy()
        expected_weighted = self.mode.mass_weight_mode(permanent=False)
        weighted_mode = VNM(2, 100.0)
        weighted_mode.set_mode([1, 2], [1, 8], expected_weighted[:, 0], expected_weighted[:, 1], expected_weighted[:, 2], is_mass_weighted=True)

        mappings = map_vnms([self.mode], [weighted_mode])
        displaced = displace_coords_with_vnm([self.mode], np.zeros((2, 3)), amplitude=1)

        self.assertAlmostEqual(mappings[0]["overlap"], 1.0)
        self.assertTrue(np.allclose(displaced, expected_weighted * 0.1))
        self.assertFalse(self.mode.is_mass_weighted)
        self.assertTrue(np.array_equal(self.mode.mode, original_mode))


if __name__ == "__main__":
    unittest.main()
