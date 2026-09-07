import unittest

import numpy as np

from scope.classes_data import Data
from scope_azo.azo_classes import State_azo


class TestAzoThermodynamics(unittest.TestCase):

    def setUp(self):
        source = type('Source', (), {'object_type': 'molecule', 'spin_multiplicity': 3, 'spin': 2, 'name': 'triplet'})()
        self.state = State_azo(source, 'thermal_test')
        self.state._z                = 1
        self.state.VNMs              = [object()]
        self.state.freqs_cm          = np.array([-100.0, 500.0])
        self.state.results['energy'] = Data('energy', -100.0, 'au')

    def test_default_thermal_settings_are_stored_in_gtot_eff(self):
        result = self.state.get_gtot_eff(temp=298.15)

        self.assertEqual(result.vib_options, {'typ': 'HO', 'FR_cutoff': None, 'FR_alpha': None, 'imaginary': 'ignore'})
        self.assertEqual(result.settings, ['svib_typ', 'imaginary', 'p_sh'])
        self.assertEqual(result.p_sh, 0.0002)
        self.assertIn('(with svib_typ=HO, imaginary=ignore and p_sh=0.0002)', repr(result))

    def test_matching_settings_can_be_extended_to_another_temperature(self):
        self.state.get_gtot_eff(temp=298.15)
        self.state.get_gtot_eff(temp=350.0)

        self.assertEqual(len(self.state.results['Gtot_eff']), 2)

    def test_changed_surface_hopping_probability_requires_overwrite(self):
        self.state.get_gtot_eff(temp=298.15, p_sh=0.0002)

        with self.assertRaises(ValueError): self.state.get_gtot_eff(temp=298.15, p_sh=0.001)

        result = self.state.get_gtot_eff(temp=298.15, p_sh=0.001, overwrite=True)
        self.assertEqual(result.p_sh, 0.001)
        self.assertEqual(len(self.state.results['Gtot_eff']), 1)


if __name__ == '__main__':
    unittest.main()
