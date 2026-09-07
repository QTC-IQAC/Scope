import unittest
from contextlib import redirect_stdout
from io import StringIO

import numpy as np

from scope.overlap import get_extended_info, overlap_molecules


class TestMoleculeOverlap(unittest.TestCase):

    def test_debug_reports_search_decisions(self):
        labels          = np.array(["O", "H"])
        coord1          = np.array([[0.0, 0.0, 0.0], [0.96, 0.0, 0.0]])
        coord2          = coord1 + np.array([2.0, -1.0, 0.5])
        adjmat          = np.array([[0, 1], [1, 0]])
        captured_output = StringIO()

        with redirect_stdout(captured_output):
            overlap_molecules(labels, coord1, labels, coord2, adjmat1=adjmat, adjmat2=adjmat, max_graph_mappings=1, debug=1)

        debug_output = captured_output.getvalue()
        self.assertIn("OVERLAP_MOLECULES: Starting overlap", debug_output)
        self.assertIn("HUNGARIAN_SEARCH: iteration=", debug_output)
        self.assertIn("GRAPH_MAPPING_SEARCH: Finished", debug_output)
        self.assertIn("OVERLAP_MOLECULES: Finished successfully", debug_output)

    def test_other_module_retains_the_public_overlap_imports(self):
        from scope.other import overlap_molecules as compatibility_import

        self.assertIs(compatibility_import, overlap_molecules)

    def test_extended_information_contains_element_labels(self):
        labels = np.array(["C", "O", "H", "H"])
        adjacency = np.array([
            [0, 1, 1, 0],
            [1, 0, 0, 1],
            [1, 0, 0, 0],
            [0, 1, 0, 0],
        ])
        adjacency_numbers = np.sum(adjacency, axis=1)

        information = get_extended_info(labels, adjacency, adjacency_numbers)

        self.assertFalse(any("None" in value for value in information))
        self.assertNotEqual(information[0], information[1])

    def test_graph_search_recovers_permuted_rotated_molecule(self):
        labels = np.array(["C", "C", "H", "H", "H", "H", "H", "H"])
        coordinates = np.array([
            [-0.75, 0.00, 0.00],
            [0.75, 0.00, 0.00],
            [-1.10, 1.00, 0.00],
            [-1.10, -0.50, np.sqrt(3) / 2],
            [-1.10, -0.50, -np.sqrt(3) / 2],
            [1.10, -1.00, 0.00],
            [1.10, 0.50, np.sqrt(3) / 2],
            [1.10, 0.50, -np.sqrt(3) / 2],
        ])
        adjacency = np.zeros((8, 8), dtype=int)
        adjacency[0, 1] = adjacency[1, 0] = 1
        for hydrogen in [2, 3, 4]:
            adjacency[0, hydrogen] = adjacency[hydrogen, 0] = 1
        for hydrogen in [5, 6, 7]:
            adjacency[1, hydrogen] = adjacency[hydrogen, 1] = 1

        permutation = np.array([1, 0, 6, 7, 5, 4, 2, 3])
        rotation = np.array([
            [0.0, -1.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
        ])
        mobile_labels = labels[permutation]
        mobile_coordinates = coordinates[permutation] @ rotation.T + np.array([3.0, -2.0, 1.0])
        mobile_adjacency = adjacency[np.ix_(permutation, permutation)]

        _, labels1, coord1, labels2, coord2, mapping = overlap_molecules(labels, coordinates, mobile_labels, mobile_coordinates, adjmat1=adjacency, adjmat2=mobile_adjacency)

        self.assertTrue(np.array_equal(labels1, labels2))
        self.assertTrue(np.allclose(coord1, coord2, atol=1e-10))
        self.assertEqual(sorted(mapping), list(range(len(labels))))

    def test_bond_orders_constrain_the_mapping(self):
        labels = np.array(["C", "O", "O", "H"])
        coordinates = np.array([
            [0.0, 0.0, 0.0],
            [1.2, 0.0, 0.0],
            [-1.0, 0.5, 0.0],
            [0.0, 1.1, 0.7],
        ])
        adjacency = np.array([
            [0, 1, 1, 1],
            [1, 0, 0, 0],
            [1, 0, 0, 0],
            [1, 0, 0, 0],
        ])
        bond_orders = adjacency.astype(float)
        bond_orders[0, 2] = bond_orders[2, 0] = 2.0

        mobile_coordinates = coordinates.copy()
        mobile_coordinates[[1, 2]] = mobile_coordinates[[2, 1]]

        _, _, _, _, _, mapping = overlap_molecules(labels, coordinates, labels, mobile_coordinates, adjmat1=adjacency, adjmat2=adjacency, bond_orders1=bond_orders, bond_orders2=bond_orders)

        self.assertEqual(mapping[1], 1)
        self.assertEqual(mapping[2], 2)

        resonance_adjacency  = np.array([
            [0, 1, 1, 0],
            [1, 0, 0, 1],
            [1, 0, 0, 0],
            [0, 1, 0, 0],
        ])
        resonance_orders1   = resonance_adjacency.astype(float)
        resonance_orders2   = resonance_adjacency.astype(float)
        resonance_orders1[0, 1] = resonance_orders1[1, 0] = 2.0
        resonance_orders2[0, 2] = resonance_orders2[2, 0] = 2.0

        with self.assertWarnsRegex(RuntimeWarning, "stored bond orders cannot be mapped exactly"):
            isgood, _, _, _, _, mapping = overlap_molecules(labels, coordinates, labels, coordinates, adjmat1=resonance_adjacency, adjmat2=resonance_adjacency, bond_orders1=resonance_orders1, bond_orders2=resonance_orders2)

        self.assertTrue(isgood)
        self.assertEqual(mapping, [0, 1, 2, 3])

    def test_terminal_metal_is_retained_for_metal_centering(self):
        labels = np.array(["Fe", "C", "H"])
        coordinates = np.array([
            [0.0, 0.0, 0.0],
            [1.8, 0.0, 0.0],
            [2.8, 0.5, 0.0],
        ])
        adjacency = np.array([
            [0, 1, 0],
            [1, 0, 1],
            [0, 1, 0],
        ])
        permutation = np.array([2, 1, 0])
        mobile_labels = labels[permutation]
        mobile_coordinates = coordinates[permutation] + np.array([2.0, -1.0, 0.5])
        mobile_adjacency = adjacency[np.ix_(permutation, permutation)]

        _, _, coord1, _, coord2, _ = overlap_molecules(labels, coordinates, mobile_labels, mobile_coordinates, center_method="metal", adjmat1=adjacency, adjmat2=mobile_adjacency)

        self.assertTrue(np.allclose(coord1, coord2, atol=1e-10))


if __name__ == "__main__":
    unittest.main()
