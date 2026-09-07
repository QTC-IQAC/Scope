import unittest

import numpy as np

from scope.classes_data import Collection, Data
from scope.classes_state import State
from scope.thermodynamics import get_Hvib, get_Svib


class TestVibrationalEntropy(unittest.TestCase):

    def test_get_Svib_records_qrrho_settings(self):
        result = get_Svib([50.0, 100.0, 500.0], 298.15, freq_units='cm', typ='QRRHO', FR_cutoff=75.0, FR_alpha=3.0)

        self.assertEqual(result.svib_typ, 'QRRHO')
        self.assertEqual(result.fr_cutoff, 75.0)
        self.assertEqual(result.fr_alpha, 3.0)
        self.assertIn('(at temperature=298.15)', repr(result))
        self.assertIn('(with svib_typ=QRRHO, fr_cutoff=75.0, fr_alpha=3.0 and imaginary=ignore)', repr(result))
        self.assertNotIn('vib_options', repr(result))
        self.assertNotIn('imaginary_frequencies', repr(result))

    def test_collection_settings_are_propagated_to_data(self):
        collection = Collection('test', 'temperature')
        data       = Data('test', 1.0, 'au')
        data.add_property('temperature', 298.15)
        collection.add_data(data)
        collection.add_setting('method', 'HO')

        self.assertEqual(collection.settings, ['method'])
        self.assertEqual(data.settings, ['method'])
        self.assertIn('(with method=HO)', repr(data))

    def test_get_Svib_rejects_unknown_model(self):
        with self.assertRaises(ValueError):
            get_Svib([100.0], 298.15, freq_units='cm', typ='QRRHO_typo')

    def test_get_Svib_retains_the_previous_default_alias(self):
        result = get_Svib([100.0], 298.15, freq_units='cm', typ='default')

        self.assertEqual(result.svib_typ, 'HO')

    def test_state_requires_overwrite_when_entropy_model_changes(self):
        source = type('Source', (), {'object_type': 'molecule', 'spin_multiplicity': 1})()
        state = State(source, 'thermal_test')
        state._z                = 1
        state.VNMs              = [object()]
        state.freqs_cm          = np.array([50.0, 100.0, 500.0])
        state.results['energy'] = Data('energy', -100.0, 'au')

        state.get_thermal_data(temp=298.15)
        HO_Svib = state.results['Svib'].datas[0].value

        with self.assertRaises(ValueError):
            state.get_thermal_data(temp=298.15, vib_options={'typ': 'QRRHO'})

        state.get_thermal_data(temp=298.15, overwrite=True, vib_options={'typ': 'QRRHO', 'FR_cutoff': 75.0, 'FR_alpha': 3.0})
        QRRHO_Svib = state.results['Svib'].datas[0]
        QRRHO_Gtot = state.results['Gtot'].datas[0]

        self.assertNotEqual(HO_Svib, QRRHO_Svib.value)
        self.assertEqual(QRRHO_Svib.svib_typ, 'QRRHO')
        self.assertEqual(QRRHO_Gtot.svib_typ, 'QRRHO')
        self.assertEqual(QRRHO_Gtot.fr_cutoff, 75.0)
        self.assertEqual(QRRHO_Gtot.fr_alpha, 3.0)
        self.assertEqual(QRRHO_Gtot.vib_options, {'typ': 'QRRHO', 'FR_cutoff': 75.0, 'FR_alpha': 3.0, 'imaginary': 'ignore'})

        state.get_thermal_data(temp=[298.15, 350.0], vib_options={'typ': 'QRRHO', 'FR_cutoff': 75.0, 'FR_alpha': 3.0})

        self.assertEqual(len(state.results['Svib']), 2)
        self.assertEqual(len(state.results['Gtot']), 2)

    def test_imaginary_frequency_treatments_are_consistent(self):
        positive_frequencies = [100.0, 500.0]
        mixed_frequencies    = [-100.0, 500.0]

        ignored_Hvib   = get_Hvib(mixed_frequencies, 298.15, freq_units='cm')
        ignored_Svib   = get_Svib(mixed_frequencies, 298.15, freq_units='cm')
        positive_Hvib  = get_Hvib([500.0], 298.15, freq_units='cm')
        positive_Svib  = get_Svib([500.0], 298.15, freq_units='cm')
        absolute_Hvib  = get_Hvib(mixed_frequencies, 298.15, freq_units='cm', imaginary='absolute')
        absolute_Svib  = get_Svib(mixed_frequencies, 298.15, freq_units='cm', imaginary='absolute')
        reference_Hvib = get_Hvib(positive_frequencies, 298.15, freq_units='cm')
        reference_Svib = get_Svib(positive_frequencies, 298.15, freq_units='cm')

        self.assertAlmostEqual(ignored_Hvib.value, positive_Hvib.value)
        self.assertAlmostEqual(ignored_Svib.value, positive_Svib.value)
        self.assertAlmostEqual(absolute_Hvib.value, reference_Hvib.value)
        self.assertAlmostEqual(absolute_Svib.value, reference_Svib.value)
        self.assertEqual(ignored_Hvib.imaginary_frequencies, [-100.0])
        self.assertEqual(ignored_Svib.imaginary_frequencies, [-100.0])

        with self.assertRaises(ValueError): get_Hvib(mixed_frequencies, 298.15, freq_units='cm', imaginary='raise')
        with self.assertRaises(ValueError): get_Svib(mixed_frequencies, 298.15, freq_units='cm', imaginary='raise')

    def test_state_requires_overwrite_when_imaginary_treatment_changes(self):
        source = type('Source', (), {'object_type': 'molecule', 'spin_multiplicity': 1})()
        state = State(source, 'imaginary_test')
        state._z                = 1
        state.VNMs              = [object()]
        state.freqs_cm          = np.array([-100.0, 500.0])
        state.results['energy'] = Data('energy', -100.0, 'au')

        state.get_thermal_data(temp=298.15)
        self.assertEqual(state.results['Hvib'].datas[0].imaginary, 'ignore')
        self.assertEqual(state.results['Svib'].datas[0].imaginary, 'ignore')

        with self.assertRaises(ValueError): state.get_thermal_data(temp=298.15, vib_options={'imaginary': 'absolute'})

        state.get_thermal_data(temp=298.15, overwrite=True, vib_options={'imaginary': 'absolute'})
        self.assertEqual(state.results['Hvib'].datas[0].imaginary, 'absolute')
        self.assertEqual(state.results['Svib'].datas[0].imaginary, 'absolute')
        self.assertEqual(state.results['Gtot'].datas[0].imaginary, 'absolute')

    def test_state_rejects_unknown_vibrational_options(self):
        source = type('Source', (), {'object_type': 'molecule', 'spin_multiplicity': 1})()
        state = State(source, 'option_test')

        with self.assertRaises(ValueError): state.get_thermal_data(vib_options={'imaginery': 'ignore'})
        with self.assertRaises(TypeError): state.get_thermal_data(vib_options='ignore')


if __name__ == '__main__':
    unittest.main()
