#!/usr/bin/env python3
"""Standard Molden -> OpenMolcas INPORB 2.2, C1, Python 3.6+.

python molden2inporb.py input.molden [-o output.INPORB]
Default: preserve atom order; automatically reorder/normalize AO components.
Use --target-log to match a different target atom order without a manual list.
Source/target must have identical radial contractions and angular representation.
"""
import argparse
import math
import os
import re
import socket
from datetime import datetime
from pathlib import Path


# Standard Molden Cartesian order, checked against its format specification.
_MOLDEN_CART = {
    0: ((0, 0, 0),),
    1: ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
    2: ((2, 0, 0), (0, 2, 0), (0, 0, 2), (1, 1, 0), (1, 0, 1), (0, 1, 1)),
    3: ((3, 0, 0), (0, 3, 0), (0, 0, 3), (1, 2, 0), (2, 1, 0),
        (2, 0, 1), (1, 0, 2), (0, 1, 2), (0, 2, 1), (1, 1, 1)),
    4: ((4, 0, 0), (0, 4, 0), (0, 0, 4), (3, 1, 0), (3, 0, 1),
        (1, 3, 0), (0, 3, 1), (1, 0, 3), (0, 1, 3), (2, 2, 0),
        (2, 0, 2), (0, 2, 2), (2, 1, 1), (1, 2, 1), (1, 1, 2)),
}
_CARTESIAN_C = {
    l: tuple(powers.index(p) for p in sorted(powers, reverse=True))
    for l, powers in _MOLDEN_CART.items()
}


def _number(value):
    result = float(value.replace('D', 'E').replace('d', 'e'))
    if not math.isfinite(result):
        raise ValueError('Non-finite number: ' + value)
    return result


def _sections(filepath):
    text = Path(filepath).read_text(encoding='utf-8-sig')
    matches = list(re.finditer(r'^\s*\[([^\]]+)\]([^\r\n]*)', text, re.M))
    sections = {}
    for i, m in enumerate(matches):
        name = m[1].strip().upper()
        if name in sections:
            raise ValueError('Duplicate section [' + name + ']')
        end = matches[i+1].start() if i+1 < len(matches) else len(text)
        sections[name] = (m[2].strip(), text[m.end():end])
    return sections


def molcas_cart_order(l):
    if l not in _CARTESIAN_C:
        raise ValueError('Unsupported angular momentum; only s through g are supported')
    return _CARTESIAN_C[l]


def _component_order(shell):
    l = shell['l']
    if shell['cart'] or l < 2:
        return molcas_cart_order(l)
    molden_m = [0] + [m for n in range(1, l+1) for m in (n, -n)]
    return tuple(molden_m.index(m) for m in range(-l, l+1))


def parse_gto_shells(filepath):
    """Original list-of-shell-lists API, with atom_id metadata added."""
    sections = _sections(filepath)
    if 'GTO' not in sections:
        raise ValueError('No [GTO] section found')
    spherical = {2: False, 3: False, 4: False}
    specified = {}
    markers = {'5D': {2: True, 3: True}, '5D7F': {2: True, 3: True},
               '5D10F': {2: True, 3: False}, '7F': {3: True}, '9G': {4: True}}
    for name in sections:
        if name in markers:
            for l, flag in markers[name].items():
                if l in specified and specified[l] != flag:
                    raise ValueError('Conflicting angular markers')
                specified[l] = spherical[l] = flag
        elif re.fullmatch(r'(?:\d+[DFGHI])+', name):
            raise ValueError('Unsupported angular marker [' + name + ']')
    lines = [line.strip() for line in sections['GTO'][1].splitlines() if line.strip()]
    atoms, atom_ids = [], set()
    i = 0
    while i < len(lines):
        t = lines[i].split()
        if re.fullmatch(r'\d+', t[0]):
            aid = int(t[0])
            if len(t) > 2 or (len(t) == 2 and t[1] != '0') or aid in atom_ids or aid < 1:
                raise ValueError('Invalid/duplicate GTO atom header: ' + lines[i])
            atom_ids.add(aid)
            atoms.append([])
            i += 1
            continue
        label = t[0].lower()
        if not atoms or label not in ('s', 'p', 'd', 'f', 'g', 'sp'):
            raise ValueError('Unsupported/malformed shell: ' + lines[i])
        if len(t) not in (2, 3):
            raise ValueError('Unsupported shell header: ' + lines[i])
        if len(t) == 3:
            # The Molden specification says this legacy field is nonfunctional.
            # The original converter also ignored it; do not rescale primitives.
            _number(t[2])
        nprim = int(t[1])
        if nprim < 1 or i+nprim >= len(lines):
            raise ValueError('Empty/truncated primitive list')
        for line in lines[i+1:i+1+nprim]:
            p = [_number(v) for v in line.split()]
            if len(p) != (3 if label == 'sp' else 2) or p[0] <= 0:
                raise ValueError('Invalid primitive row: ' + line)
        for ang in ('s', 'p') if label == 'sp' else (label,):
            l = 'spdfg'.index(ang)
            cart = not spherical.get(l, False)
            atoms[-1].append({'l': l, 'nprim': nprim, 'ncontr': 1, 'cart': cart,
                             'ncomp': (l+1)*(l+2)//2 if cart else 2*l+1, 'atom_id': aid})
        i += nprim+1
    if not atoms or any(not shells for shells in atoms):
        raise ValueError('Missing basis shells')
    return atoms


def _target_rows(atoms, atom_order=None):
    """Yield (source AO index, shell, source component) in target AO order."""
    ids = [shells[0].get('atom_id', i+1) for i, shells in enumerate(atoms)]
    order = ids if atom_order is None else list(atom_order)
    if sorted(order) != sorted(ids):
        raise ValueError('Atom order must contain each GTO atom ID exactly once')
    data = {}
    base = 0
    for aid, shells in zip(ids, atoms):
        groups = {}
        for sh in shells:
            groups.setdefault(sh['l'], []).append((base, sh))
            base += sh['ncomp']
        data[aid] = groups
    for aid in order:
        for l, group in sorted(data[aid].items()):
            if len({(sh['cart'], sh['ncomp']) for _, sh in group}) != 1:
                raise ValueError('Mixed representations at the same angular momentum')
            for comp in _component_order(group[0][1]):
                for start, sh in group:
                    yield start+comp, sh, comp


def build_index_map(atoms, atom_order=None):
    """Source-to-target AO indices; original call signature retained."""
    idx_map = [0]*sum(sh['ncomp'] for shells in atoms for sh in shells)
    for target, (source, _, _) in enumerate(_target_rows(atoms, atom_order)):
        idx_map[source] = target
    return idx_map


def compute_d_pure_positions(atoms, atom_order=None):
    """Retained API: Cartesian xx/yy/zz target positions, never spherical d."""
    return {i for i, (_, sh, comp) in enumerate(_target_rows(atoms, atom_order))
            if sh['l'] == 2 and sh['cart'] and comp in (0, 1, 2)}


def normalization_divisors(atoms, atom_order=None):
    divisors = []
    for _, sh, comp in _target_rows(atoms, atom_order):
        factor = 1
        if sh['cart']:
            for n in _MOLDEN_CART[sh['l']][comp]:
                for k in range(2*n-1, 0, -2):
                    factor *= k
        divisors.append(math.sqrt(factor))
    return divisors


def normalize_orbitals(orbitals, atoms, atom_order=None):
    divisors = normalization_divisors(atoms, atom_order)
    for orb in orbitals:
        for i, divisor in enumerate(divisors):
            if divisor != 1.0:
                # Preserve the original d division's exact floating-point result.
                orb['coeffs'][i] /= divisor


def parse_molden_mo(filepath, idx_map):
    sections = _sections(filepath)
    if 'MO' not in sections:
        raise ValueError('No [MO] section found in Molden file')
    nbasis = len(idx_map)
    if not nbasis or sorted(idx_map) != list(range(nbasis)):
        raise ValueError('Invalid AO permutation')
    text = sections['MO'][1]
    starts = [m.start() for m in re.finditer(r'^\s*Sym\s*=', text, re.M | re.I)]
    if not starts or text[:starts[0]].strip():
        raise ValueError('MO must start with Sym=')
    starts.append(len(text))
    orbitals = []
    # No zero-width re.split: that is not supported consistently by Python 3.6.
    for start, end in zip(starts, starts[1:]):
        block = text[start:end]
        fields = {k.lower(): v.strip() for k, v in re.findall(
            r'^\s*(Sym|Ene|Spin|Occup)\s*=\s*([^\r\n]+)', block, re.M | re.I)}
        if fields.get('spin', 'Alpha').lower() != 'alpha':
            raise ValueError('Beta/UHF orbitals require an unrestricted writer')
        coeffs = [0.0]*nbasis
        seen = set()
        for line in block.splitlines():
            if not line.strip() or re.match(r'^\s*(Sym|Ene|Spin|Occup)\s*=', line, re.I):
                continue
            t = line.split()
            if len(t) != 2:
                raise ValueError('Invalid MO coefficient: ' + line)
            source = int(t[0])-1
            if not 0 <= source < nbasis or source in seen:
                raise ValueError('Duplicate/out-of-range MO coefficient index: ' + t[0])
            seen.add(source)
            coeffs[idx_map[source]] = _number(t[1])
        if not seen:
            raise ValueError('MO has no coefficients')
        # Preserve sparse-zero input behavior: omitted AO coefficients are zero.
        orbitals.append({'occup': _number(fields.get('occup', '0')),
                         'ene': _number(fields.get('ene', '0')), 'coeffs': coeffs})
    if not orbitals or len(orbitals) > nbasis:
        raise ValueError('Invalid number of restricted orbitals')
    return nbasis, orbitals


def target_order_from_log(filepath, molden_path, atoms):
    """Match element/coordinates to the first C1 log geometry; not a basis proof."""
    text = Path(filepath).read_text(errors='replace')
    if not re.search(r'Character Table for C1\b', text):
        raise ValueError('Target log must have C1 symmetry')
    section = _sections(molden_path).get('ATOMS')
    if section is None or section[0].upper() not in ('ANGS', 'AU'):
        raise ValueError('Atom matching requires [Atoms] Angs or AU')
    factor = 0.529177209 if section[0].upper() == 'AU' else 1.0
    coords = {}
    for line in section[1].splitlines():
        if not line.strip():
            continue
        t = line.split()
        if len(t) != 6 or not re.match(r'[A-Za-z]+', t[0]):
            raise ValueError('Invalid atom row: ' + line)
        aid = int(t[1])
        if aid in coords:
            raise ValueError('Duplicate atom ID in [Atoms]')
        coords[aid] = (re.match(r'[A-Za-z]+', t[0]).group().capitalize(),
                       [_number(v)*factor for v in t[3:]])
    match = re.search(r'Cartesian Coordinates / Bohr, Angstrom', text)
    if not match:
        raise ValueError('Target coordinate table not found')
    order = []
    for line in text[match.end():].splitlines():
        t = line.split()
        if len(t) == 8 and t[0].isdigit() and re.fullmatch(r'[A-Za-z]+\d*', t[1]):
            symbol = re.match(r'[A-Za-z]+', t[1]).group().capitalize()
            xyz = [_number(v) for v in t[-3:]]
            matches = [aid for aid, (element, pos) in coords.items() if element == symbol
                       and max(abs(x-y) for x, y in zip(pos, xyz)) < 2e-6]
            if len(matches) != 1 or matches[0] in order:
                raise ValueError('Ambiguous/mismatched target atom: ' + line)
            order.append(matches[0])
        elif order:
            break
    ids = [a[0]['atom_id'] for a in atoms]
    if sorted(order) != sorted(ids) or sorted(coords) != sorted(ids):
        raise ValueError('Target atoms do not match Molden/GTO atoms')
    nbasis = sum(s['ncomp'] for a in atoms for s in a)
    counts = [int(n) for n in re.findall(r'(?:Number of basis functions|Basis functions)\s+(\d+)', text)]
    if not counts or any(n != nbasis for n in counts):
        raise ValueError('Target/Molden AO counts differ or target count is absent')
    return order


def fmt_es(val, width=22):
    if abs(val) < 1e-99:
        return f'{0.0:{width}.14E}'
    return re.sub(r'E([+-])0(\d{2})', r'E\1\2', f'{val:{width}.14E}')


def write_coeffs_block(fh, coeffs, per_line=5):
    for i in range(0, len(coeffs), per_line):
        fh.write(''.join(fmt_es(c) for c in coeffs[i:i+per_line]).rstrip() + '\n')


def write_gssorb(filepath, nbasis, orbitals, write_index=True):
    if not orbitals or any(len(o['coeffs']) != nbasis for o in orbitals):
        raise ValueError('Output MO dimensions do not match NBAS')
    with open(filepath, 'w') as f:
        f.write('#INPORB 2.2\n#INFO\n* SCF orbitals\n       0       1       2\n')
        # NBAS precedes NORB; equal in the legacy square case.
        f.write(f'      {nbasis}\n      {len(orbitals)}\n')
        date_str = datetime.now().strftime('%c')
        f.write(f'*BC:HOST {socket.gethostname()} PID {os.getpid()} DATE {date_str}\n#ORB\n')
        for i, orb in enumerate(orbitals, 1):
            f.write(f'* ORBITAL    1     {i:>2}\n')
            write_coeffs_block(f, orb['coeffs'])
        occ_values = [orb['occup'] for orb in orbitals]
        f.write('#OCC\n* OCCUPATION NUMBERS\n')
        write_coeffs_block(f, occ_values)
        f.write('#OCHR\n* OCCUPATION NUMBERS (HUMAN-READABLE)\n')
        for i in range(0, len(occ_values), 10):
            f.write(''.join(f'{v:8.4f}' for v in occ_values[i:i+10]).rstrip() + '\n')
        f.write('#ONE\n* ONE ELECTRON ENERGIES\n')
        ene_values = [orb['ene'] for orb in orbitals]
        # INPORB 2.2 uses 10(1X,ES11.4) here, not the MO coefficient format.
        # Keep real zero energies; do not invent a replacement spectrum.
        for i in range(0, len(ene_values), 10):
            f.write(''.join(f' {e:11.4E}' for e in ene_values[i:i+10]) + '\n')
        if write_index:
            f.write('#INDEX\n* 1234567890\n')
            # Preserve the original default. Natural occupations do not uniquely
            # identify an active space; use --no-index and explicit RASSCF input.
            pattern = ''.join('i' if occ > 0.5 else 's' for occ in occ_values)
            for row, i in enumerate(range(0, len(pattern), 10)):
                f.write(f'{row % 10} {pattern[i:i+10]}\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input_file')
    parser.add_argument('-o', '--output', help='Default: input basename + .INPORB')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--atom-order', help='Optional comma-separated Molden atom IDs in target order')
    group.add_argument('--target-log', help='Automatically match atoms against a C1 Molcas log')
    group.add_argument('--same-atom-order', action='store_true', help='Explicitly select the default; no reordering of atoms')
    parser.add_argument('--no-index', action='store_true', help='Omit occupation-based #INDEX; set the active space explicitly')
    args = parser.parse_args()
    output = args.output or os.path.splitext(args.input_file)[0] + '.INPORB'
    try:
        if Path(args.input_file).resolve() == Path(output).resolve():
            raise ValueError('Output must not overwrite the input Molden file')
        atoms = parse_gto_shells(args.input_file)
        order = None
        if args.atom_order is not None:
            order = [int(s.strip()) for s in args.atom_order.split(',')]
        elif args.target_log:
            order = target_order_from_log(args.target_log, args.input_file, atoms)
        idx_map = build_index_map(atoms, order)
        nbasis, orbitals = parse_molden_mo(args.input_file, idx_map)
        normalize_orbitals(orbitals, atoms, order)
        print(f'Input:  {args.input_file}\nOutput: {output}')
        print(f'Parsed {len(atoms)} atoms with shell structure:')
        for iatom, shells in enumerate(atoms, 1):
            for sh in shells:
                kind = 'Cartesian' if sh['cart'] else 'Spherical'
                print(f"  atom {iatom}: l={sh['l']}, nprim={sh['nprim']}, "
                      f"ncontr={sh['ncontr']}, {kind} ({sh['ncomp']} comp)")
        if order is not None:
            print('Target atom order (Molden IDs):', ','.join(map(str, order)))
        print(f'Total basis functions: {nbasis}\nTotal orbitals: {len(orbitals)}')
        write_gssorb(output, nbasis, orbitals, write_index=not args.no_index)
        print('Done.')
    except (ValueError, OSError) as error:
        parser.exit(2, 'Error: ' + str(error) + '\n')


if __name__ == '__main__':
    main()
