"""Offline regression suite, including exact comparison with commit 78cf0a8.
Run with Python 3.6+: python -W error -m unittest discover -s tests -v
Synthetic fixtures only; no user calculations or external dependencies.
"""
import importlib.util
import math
import random
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import molden2inporb as current

spec = importlib.util.spec_from_file_location('legacy', str(ROOT/'tests/fixtures/legacy_molden2inporb.py'))
legacy = importlib.util.module_from_spec(spec)
spec.loader.exec_module(legacy)


def fixture(shells_by_atom, marker='', sparse=False, nmo=None, atoms=True):
    """Each shell is (label, component count); deterministic, nontrivial MOs."""
    lines = ['[Molden Format]']
    if atoms:
        lines += ['[Atoms] Angs']
        lines += ['H {} 1 {} 0 0'.format(i+1, i) for i in range(len(shells_by_atom))]
    if marker:
        lines += marker.splitlines()
    lines += ['[GTO]']
    nb = 0
    for i, shells in enumerate(shells_by_atom, 1):
        lines += ['{} 0'.format(i)]
        for label, size in shells:
            lines += ['{} 2 1.00'.format(label), '2.0 0.7' + (' 0.6' if label == 'sp' else ''),
                      '0.5 0.3' + (' 0.4' if label == 'sp' else '')]
            nb += size
        lines += ['']
    lines += ['[MO]']
    rng = random.Random(819)
    for mo in range(nb if nmo is None else nmo):
        lines += [' Sym= {}.1'.format(mo+1), ' Ene= {:.8f}'.format(-5+mo*.13),
                  ' Spin= Alpha', ' Occup= {}'.format(2 if mo < nb//2 else 0)]
        for i in range(nb):
            value = rng.uniform(-2, 2)
            if not sparse or i % 3 == 0:
                lines += [' {} {:.17g}'.format(i+1, value)]
    return '\n'.join(lines)+'\n'


def converted(module, path, order=None):
    atoms = module.parse_gto_shells(path)
    indices = module.build_index_map(atoms) if order is None else module.build_index_map(atoms, order)
    nb, orbs = module.parse_molden_mo(path, indices)
    if module is legacy:
        for orb in orbs:
            for i in module.compute_d_pure_positions(atoms):
                orb['coeffs'][i] /= math.sqrt(3)
    else:
        module.normalize_orbitals(orbs, atoms, order)
    return nb, orbs


def sections(path):
    result = {}
    key = None
    for line in Path(path).read_text().splitlines():
        if line.startswith('#'):
            key = line
            result[key] = []
        elif key is not None:
            result[key].append(line)
    return result


class ConverterTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def source(self, text):
        path = self.root/'input.molden'
        path.write_text(text)
        return path

    def run_cli(self, path, *args):
        return subprocess.run([sys.executable, '-W', 'error', str(ROOT/'molden2inporb.py'),
                               str(path)] + list(args), stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, universal_newlines=True)

    def test_exact_legacy_coefficients_across_spd_contractions(self):
        # 48 layouts, plus a second atom with different contractions each time.
        for ns in range(1, 4):
            for np in range(4):
                for nd in range(4):
                    with self.subTest(ns=ns, np=np, nd=nd):
                        path = self.source(fixture([[('s', 1)]*ns+[('p', 3)]*np+[('d', 6)]*nd,
                                                    [('s', 1), ('s', 1), ('p', 3)]]))
                        before = converted(legacy, path)
                        after = converted(current, path)
                        self.assertEqual(before, after)  # exact float equality, no tolerance
                        a = current.parse_gto_shells(path)
                        b = legacy.parse_gto_shells(path)
                        self.assertEqual(current.build_index_map(a), legacy.build_index_map(b))
                        self.assertEqual(current.compute_d_pure_positions(a), legacy.compute_d_pure_positions(b))

    def test_legacy_output_sections_unchanged(self):
        path = self.source(fixture([[('s', 1)]*3+[('p', 3)]*2+[('d', 6)]*2]))
        nb, old = converted(legacy, path)
        _, new = converted(current, path)
        legacy.write_gssorb(str(self.root/'old'), nb, old)
        current.write_gssorb(str(self.root/'new'), nb, new)
        before, after = sections(self.root/'old'), sections(self.root/'new')
        for key in ('#ORB', '#OCC', '#OCHR', '#INDEX'):
            self.assertEqual(before[key], after[key])
        self.assertEqual(before['#INFO'][:4], after['#INFO'][:4])

    def test_sparse_and_no_atoms_legacy_input(self):
        path = self.source(fixture([[('s', 1), ('p', 3), ('d', 6)]], sparse=True, atoms=False))
        self.assertEqual(converted(legacy, path), converted(current, path))

    def test_legacy_fractional_occupations_and_unused_scale(self):
        text = fixture([[('s', 1), ('p', 3), ('d', 6)]])
        text = text.replace('Occup= 2', 'Occup= 0.5001').replace('Occup= 0\n', 'Occup= 0.5\n')
        text = text.replace('2 1.00', '2 2.00')
        path = self.source(text)
        nb, old = converted(legacy, path)
        _, new = converted(current, path)
        self.assertEqual(old, new)
        legacy.write_gssorb(str(self.root/'old'), nb, old)
        current.write_gssorb(str(self.root/'new'), nb, new)
        self.assertEqual(sections(self.root/'old')['#INDEX'], sections(self.root/'new')['#INDEX'])

    def test_default_cli_only_needs_molden(self):
        path = self.source(fixture([[('s', 1), ('p', 3), ('d', 6)]]))
        result = self.run_cli(path)
        self.assertEqual(result.returncode, 0, result.stderr)
        expected = path.with_suffix('.INPORB')
        self.assertTrue(expected.is_file())
        data = sections(expected)
        self.assertIn('#INDEX', data)
        self.assertEqual(int(data['#INFO'][2]), 10)

    def test_explicit_same_atom_order_equals_default(self):
        path = self.source(fixture([[('s', 1)], [('s', 1), ('p', 3)]]))
        self.assertEqual(self.run_cli(path).returncode, 0)
        other = self.root/'explicit.INPORB'
        result = self.run_cli(path, '-o', str(other), '--same-atom-order')
        self.assertEqual(result.returncode, 0, result.stderr)
        a, b = sections(path.with_suffix('.INPORB')), sections(other)
        for key in ('#ORB', '#OCC', '#OCHR', '#ONE', '#INDEX'):
            self.assertEqual(a[key], b[key])

    def shell(self, label, size, marker=''):
        path = self.source(fixture([[(label, size)]], marker))
        return current.parse_gto_shells(path)

    def test_f_reference_order_and_divisors(self):
        a = self.shell('f', 10)
        self.assertEqual(current.molcas_cart_order(3), (0, 4, 5, 3, 9, 6, 1, 8, 7, 2))
        self.assertEqual(current.normalization_divisors(a), [math.sqrt(n) for n in (15, 3, 3, 3, 1, 3, 15, 3, 3, 15)])

    def test_g_reference_order_and_divisors(self):
        a = self.shell('g', 15)
        self.assertEqual(current.molcas_cart_order(4), (0, 3, 4, 9, 12, 10, 5, 13, 14, 7, 1, 6, 11, 8, 2))
        self.assertEqual(current.normalization_divisors(a), [math.sqrt(n) for n in (105, 15, 15, 9, 3, 9, 15, 3, 3, 15, 105, 15, 9, 15, 105)])

    def test_spherical_d_no_cartesian_scaling(self):
        a = self.shell('d', 5, '[5D]')
        self.assertEqual(current.build_index_map(a), [2, 3, 1, 4, 0])
        self.assertEqual(current.normalization_divisors(a), [1]*5)
        self.assertEqual(current.compute_d_pure_positions(a), set())

    def test_spherical_f_and_g(self):
        for label, size, marker, take in (
                ('f', 7, '[5D]', [6, 4, 2, 0, 1, 3, 5]),
                ('g', 9, '[5D7F]\n[9G]', [8, 6, 4, 2, 0, 1, 3, 5, 7])):
            a = self.shell(label, size, marker)
            idx = current.build_index_map(a)
            self.assertEqual([idx.index(i) for i in range(size)], take)
            self.assertEqual(current.normalization_divisors(a), [1]*size)

    def test_mixed_representation_marker(self):
        a = self.shell('f', 10, '[5D10F]')
        self.assertTrue(a[0][0]['cart'])
        self.assertEqual(a[0][0]['ncomp'], 10)
        a = self.shell('d', 6, '[7F]')
        self.assertTrue(a[0][0]['cart'])

    def test_sp_and_nonadjacent_contractions(self):
        path = self.source(fixture([[('sp', 4), ('sp', 4)]]))
        a = current.parse_gto_shells(path)
        idx = current.build_index_map(a)
        self.assertEqual([idx.index(i) for i in range(8)], [0, 4, 1, 5, 2, 6, 3, 7])

    def test_atom_reorder_is_opt_in(self):
        path = self.source(fixture([[('s', 1)], [('p', 3)]]))
        a = current.parse_gto_shells(path)
        self.assertEqual(current.build_index_map(a), [0, 1, 2, 3])
        self.assertEqual(current.build_index_map(a, [2, 1]), [3, 0, 1, 2])
        with self.assertRaises(ValueError):
            current.build_index_map(a, [2, 2])

    def test_target_log_automatically_matches_atoms(self):
        path = self.source(fixture([[('s', 1)], [('p', 3)]]))
        log = self.root/'target.log'
        log.write_text('Character Table for C1\nCartesian Coordinates / Bohr, Angstrom\n'
                       'Center Label x y z x y z\n1 H1 1.889726 0 0 1 0 0\n'
                       '2 H2 0 0 0 0 0 0\n\nNumber of basis functions 4\n')
        a = current.parse_gto_shells(path)
        self.assertEqual(current.target_order_from_log(log, path, a), [2, 1])
        output = self.root/'matched.INPORB'
        result = self.run_cli(path, '--target-log', str(log), '-o', str(output))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('2,1', result.stdout)
        log.write_text(log.read_text().replace('functions 4', 'functions 5'))
        with self.assertRaises(ValueError):
            current.target_order_from_log(log, path, a)

    def test_fixed_width_energy_reader_and_rectangular_dimensions(self):
        orbs = [dict(coeffs=[1.0]+[0.0]*12, occup=2., ene=-3.123456),
                dict(coeffs=[0.0, 1.0]+[0.0]*11, occup=0., ene=0.)]
        out = self.root/'rect.INPORB'
        current.write_gssorb(out, 13, orbs)
        data = sections(out)
        self.assertEqual([int(s) for s in data['#INFO'][2:4]], [13, 2])
        raw = data['#ONE'][1]
        self.assertEqual(len(raw), 24)
        self.assertEqual([float(raw[i:i+12]) for i in (0, 12)], [-3.1235, 0.])

    def test_zero_energies_and_single_orbital(self):
        out = self.root/'zero.INPORB'
        current.write_gssorb(out, 1, [dict(coeffs=[1.], occup=2., ene=0.)])
        self.assertEqual(float(sections(out)['#ONE'][1]), 0.)

    def test_index_labels_wrap_and_optional_omission(self):
        orbs = [dict(coeffs=[0.]*110, occup=2., ene=-1.) for _ in range(110)]
        out = self.root/'many.INPORB'
        current.write_gssorb(out, 110, orbs)
        index = sections(out)['#INDEX'][1:]
        self.assertEqual([s[0] for s in index], list('01234567890'))
        self.assertTrue(all(len(s) == 12 for s in index))
        current.write_gssorb(out, 110, orbs, write_index=False)
        self.assertNotIn('#INDEX', sections(out))

    def test_fortran_numbers_and_multiple_mos_python36(self):
        text = fixture([[('s', 1), ('p', 3)]], nmo=3).replace('Ene= -5.00000000', 'Ene= -5.0D+00')
        path = self.source(text)
        nb, orbs = converted(current, path)
        self.assertEqual((nb, len(orbs), orbs[0]['ene']), (4, 3, -5.))

    def test_bad_inputs_fail_before_output(self):
        text = fixture([[('s', 1)]])
        for bad in (text+'1 0.2\n', text+'2 0.2\n', text.replace('Spin= Alpha', 'Spin= Beta'),
                    text.replace('s 2 1.00', 'h 2 1.00'), text.replace('2.0 0.7', '-2.0 0.7')):
            with self.subTest(text=bad):
                path = self.source(bad)
                result = self.run_cli(path)
                self.assertNotEqual(result.returncode, 0)
                self.assertFalse(path.with_suffix('.INPORB').exists())
        path = self.source(text)
        self.assertNotEqual(self.run_cli(path, '-o', str(path)).returncode, 0)
        self.assertEqual(path.read_text(), text)


if __name__ == '__main__':
    unittest.main()
