import unittest

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from scope.classes_specie import import_molecule, import_rdkit_molecule


class TestRDKitImport(unittest.TestCase):

    def test_implicit_hydrogens_are_materialized_for_a_3D_molecule(self):
        rdkit_molecule = Chem.MolFromSmiles('CC')
        conformer      = Chem.Conformer(rdkit_molecule.GetNumAtoms())
        conformer.SetAtomPosition(0, (0.0, 0.0, 0.0))
        conformer.SetAtomPosition(1, (1.54, 0.0, 0.0))
        conformer.Set3D(True)
        rdkit_molecule.AddConformer(conformer)

        self.assertEqual(rdkit_molecule.GetNumAtoms(), 2)
        self.assertEqual(sum(atom.GetNumImplicitHs() for atom in rdkit_molecule.GetAtoms()), 6)

        imported_molecule = import_molecule(rdkit_molecule)

        self.assertEqual(rdkit_molecule.GetNumAtoms(), 2)
        self.assertEqual(imported_molecule.natoms, 8)
        self.assertEqual(imported_molecule.labels.count('H'), 6)
        self.assertEqual(imported_molecule.rdkit_obj.GetNumAtoms(), 8)
        self.assertEqual(sum(atom.GetNumImplicitHs() for atom in imported_molecule.rdkit_obj.GetAtoms()), 0)

    def test_formal_atomic_charges_survive_the_complete_import(self):
        cases = [
            ('CCO', 0),
            ('[NH4+]', 1),
            ('CC(=O)[O-]', -1),
            ('[NH3+]CC(=O)[O-]', 0),
        ]

        for smiles, expected_charge in cases:
            with self.subTest(smiles=smiles):
                rdkit_molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
                embed_status   = AllChem.EmbedMolecule(rdkit_molecule, randomSeed=7)
                self.assertEqual(embed_status, 0)
                rdkit_atomic_charges = [atom.GetFormalCharge() for atom in rdkit_molecule.GetAtoms()]

                intermediate_molecule = import_rdkit_molecule(rdkit_molecule)
                imported_molecule     = import_molecule(rdkit_molecule)

                self.assertEqual(intermediate_molecule.atomic_charges, rdkit_atomic_charges)
                self.assertEqual(imported_molecule.atomic_charges, rdkit_atomic_charges)
                self.assertEqual(intermediate_molecule.charge, expected_charge)
                self.assertEqual(imported_molecule.charge, expected_charge)

    def test_rdkit_adjacency_is_preserved(self):
        rdkit_molecule = Chem.AddHs(Chem.MolFromSmiles('CCO'))
        embed_status   = AllChem.EmbedMolecule(rdkit_molecule, randomSeed=7)
        self.assertEqual(embed_status, 0)

        imported_molecule = import_molecule(rdkit_molecule)
        rdkit_adjacency    = Chem.GetAdjacencyMatrix(rdkit_molecule).astype(int)

        self.assertTrue(np.array_equal(imported_molecule.adjmat, rdkit_adjacency))
        self.assertTrue(np.array_equal(imported_molecule.adjnum, rdkit_adjacency.sum(axis=1)))

    def test_fragmented_rdkit_molecule_is_rejected(self):
        rdkit_molecule = Chem.AddHs(Chem.MolFromSmiles('CC.O=C=O'))
        embed_status   = AllChem.EmbedMolecule(rdkit_molecule, randomSeed=7)
        self.assertEqual(embed_status, 0)

        with self.assertRaisesRegex(ValueError, 'disconnected fragments'):
            import_molecule(rdkit_molecule)


if __name__ == '__main__':
    unittest.main()
