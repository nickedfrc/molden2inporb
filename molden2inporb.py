#!/usr/bin/env python3
"""
Convert Molden orbital file to Molcas INPORB format.

Parses [GTO] section to determine actual shell structure,
then reorders MO coefficients to match Molcas SEWARD convention:
  atom → shell → angular component → contracted function

Usage:
    python molden2gssorb.py input.molden
    python molden2gssorb.py input.molden -o output.INPORB
"""

import re
import os
import sys
import socket
import math
from datetime import datetime

# ---------------------------------------------------------------------------
# Molcas-internal Cartesian ordering
# ---------------------------------------------------------------------------
_CARTESIAN_C = {
    1: (0, 1, 2),
    2: (0, 3, 4, 1, 5, 2),
    3: (0, 3, 4, 5, 9, 7, 1, 6, 8, 2),
}


def molcas_cart_order(l):
    return _CARTESIAN_C.get(l, tuple(range((l + 1) * (l + 2) // 2)))


# ---------------------------------------------------------------------------
# [GTO] parser
# ---------------------------------------------------------------------------

def parse_gto_shells(filepath):
    with open(filepath) as f:
        text = f.read()

    sph_markers = set()
    for m in ('[5D]', '[5d]', '[7F]', '[7f]', '[9G]', '[9g]'):
        if m in text:
            sph_markers.add(m)

    gto = text.split('[GTO]', 1)[1]
    gto = gto.split('[', 1)[0]

    atoms = []
    shells = []
    iline = 0
    lines = [l.strip() for l in gto.strip().splitlines() if l.strip()]
    nlines = len(lines)
    l_map = {'s': 0, 'p': 1, 'd': 2, 'f': 3, 'g': 4, 'h': 5, 'i': 6}

    while iline < nlines:
        line = lines[iline]
        m = re.match(r'^(\d+)\s+\d+$', line)
        if m:
            if shells:
                atoms.append(shells)
                shells = []
            iline += 1
            continue

        m_shell = re.match(r'^([spdfghi])\s+(\d+)\s+([\d.]+)$', line)
        if m_shell:
            ang = m_shell.group(1)
            nprim = int(m_shell.group(2))
            l = l_map[ang]
            ncontr = 1

            marker = ''
            if l == 2 and ('[5D]' in sph_markers or '[5d]' in sph_markers):
                marker = '5D'
            elif l == 3 and ('[7F]' in sph_markers or '[7f]' in sph_markers):
                marker = '7F'
            elif l == 4 and ('[9G]' in sph_markers or '[9g]' in sph_markers):
                marker = '9G'

            cart = (marker == '')
            ncomp = (l + 1) * (l + 2) // 2 if cart else 2 * l + 1

            shells.append({
                'l': l, 'nprim': nprim, 'ncontr': ncontr,
                'cart': cart, 'ncomp': ncomp,
            })
            iline += 1
            iline += nprim
        else:
            iline += 1

    if shells:
        atoms.append(shells)

    return atoms


def compute_d_pure_positions(atoms):
    positions = set()
    pos = 0
    for shells in atoms:
        groups = []
        for sh in shells:
            l = sh['l']
            if groups and groups[-1]['l'] == l:
                groups[-1]['ncontr'] += 1
            else:
                groups.append({'l': l, 'ncomp': sh['ncomp'], 'ncontr': 1})
        for grp in groups:
            ncontr = grp['ncontr']
            ncomp = grp['ncomp']
            n = ncomp * ncontr
            if grp['l'] == 2:
                for ic in (0, 3, 5):
                    for j in range(pos + ic * ncontr, pos + (ic + 1) * ncontr):
                        positions.add(j)
            pos += n
    return positions


def build_index_map(atoms):
    idx_map = []
    m_base = 0
    g_base = 0

    for shells in atoms:
        groups = []
        for sh in shells:
            l = sh['l']
            if groups and groups[-1]['l'] == l:
                groups[-1]['ncontr'] += 1
            else:
                groups.append({
                    'l': l, 'cart': sh['cart'],
                    'ncomp': sh['ncomp'], 'ncontr': 1,
                })

        for grp in groups:
            l = grp['l']
            cart = grp['cart']
            ncomp = grp['ncomp']
            ncontr = grp['ncontr']

            if cart:
                perm = molcas_cart_order(l)
            else:
                perm = tuple(range(ncomp))
            comp_map = {m: g for g, m in enumerate(perm)}

            for ic in range(ncontr):
                for mc in range(ncomp):
                    molden_pos = m_base + ic * ncomp + mc
                    gc = comp_map[mc]
                    gssorb_pos = g_base + gc * ncontr + ic
                    if molden_pos >= len(idx_map):
                        idx_map.extend([0] * (molden_pos + 1 - len(idx_map)))
                    idx_map[molden_pos] = gssorb_pos

            ntotal = ncontr * ncomp
            m_base += ntotal
            g_base += ntotal

    return idx_map


# ---------------------------------------------------------------------------
# Molden MO parser
# ---------------------------------------------------------------------------

def parse_molden_mo(filepath, idx_map):
    with open(filepath) as f:
        text = f.read()

    mo_section = text.split('[MO]')
    if len(mo_section) < 2:
        raise ValueError("No [MO] section found in Molden file")
    mo_text = mo_section[1]

    parts = re.split(r'(\n\s*)(?=Sym=)', mo_text)
    blocks = []
    for p in parts:
        if p.strip().startswith('Sym='):
            blocks.append(p)
        elif p.strip() and not p.strip().startswith('Sym='):
            if blocks:
                blocks[-1] += p

    orbitals = []
    for block in blocks:
        block = block.strip()
        if not block:
            continue

        occup_match = re.search(r'Occup=\s*([\d.\-Ee+]+)', block)
        ene_match = re.search(r'Ene=\s*([\d.\-Ee+]+)', block)
        occup = float(occup_match.group(1)) if occup_match else 0.0
        ene = float(ene_match.group(1)) if ene_match else 0.0

        coeffs_by_index = {}
        for line in block.splitlines():
            m = re.match(r'^\s*(\d+)\s+([\d.\-EeDd+]+)', line)
            if m:
                idx = int(m.group(1))
                val = float(m.group(2).replace('D', 'E'))
                coeffs_by_index[idx] = val

        nbasis = len(idx_map)
        coeffs = [coeffs_by_index.get(i, 0.0) for i in range(1, nbasis + 1)]
        reordered = [0.0] * nbasis
        for molden_i, gssorb_i in enumerate(idx_map):
            reordered[gssorb_i] = coeffs[molden_i]
        orbitals.append({'occup': occup, 'ene': ene, 'coeffs': reordered})

    return nbasis, orbitals


# ---------------------------------------------------------------------------
# InpOrb writer
# ---------------------------------------------------------------------------

def fmt_es(val, width=22):
    if abs(val) < 1e-99:
        return f"{0.0:{width}.14E}"
    s = f"{val:{width}.14E}"
    s = re.sub(r'E([+-])0(\d{2})', r'E\1\2', s)
    return s


def write_coeffs_block(fh, coeffs, per_line=5):
    for i in range(0, len(coeffs), per_line):
        chunk = coeffs[i:i + per_line]
        line = ''.join(fmt_es(c) for c in chunk)
        fh.write(line.rstrip() + '\n')


def write_gssorb(filepath, nbasis, orbitals):
    hostname = socket.gethostname()
    pid = os.getpid()
    date_str = datetime.now().strftime('%c')

    with open(filepath, 'w') as f:
        f.write('#INPORB 2.2\n')
        f.write('#INFO\n')
        f.write('* SCF orbitals\n')
        f.write('       0       1       2\n')
        f.write(f'      {len(orbitals)}\n')
        f.write(f'      {nbasis}\n')
        f.write(f'*BC:HOST {hostname} PID {pid} DATE {date_str}\n')
        f.write('#ORB\n')

        for i, orb in enumerate(orbitals, 1):
            f.write(f'* ORBITAL    1     {i:>2}\n')
            write_coeffs_block(f, orb['coeffs'])

        f.write('#OCC\n')
        f.write('* OCCUPATION NUMBERS\n')
        write_coeffs_block(f, [orb['occup'] for orb in orbitals])

        f.write('#OCHR\n')
        f.write('* OCCUPATION NUMBERS (HUMAN-READABLE)\n')
        occ_values = [orb['occup'] for orb in orbitals]
        for i in range(0, len(occ_values), 10):
            chunk = occ_values[i:i + 10]
            line = ''.join(f'{v:8.4f}' for v in chunk)
            f.write(line.rstrip() + '\n')

        f.write('#ONE\n')
        f.write('* ONE ELECTRON ENERGIES\n')
        ene_values = [orb['ene'] for orb in orbitals]
        if all(abs(e) < 1e-10 for e in ene_values):
            n = len(ene_values)
            ene_values = [round(-20.0 + i * (25.0 / (n - 1)), 4) for i in range(n)]
        write_coeffs_block(f, ene_values)

        f.write('#INDEX\n')
        f.write('* 1234567890\n')
        occ_pattern = ''.join('i' if occ > 0.5 else 's' for occ in occ_values)
        for line_num in range(0, (len(occ_pattern) + 9) // 10):
            chunk = occ_pattern[line_num * 10:(line_num + 1) * 10]
            f.write(f'{line_num} {chunk}\n')


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} input.molden [-o output.INPORB]")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = None
    args = sys.argv[2:]
    i = 0
    while i < len(args):
        if args[i] == '-o' and i + 1 < len(args):
            output_file = args[i + 1]
            i += 2
        else:
            i += 1

    if output_file is None:
        base = os.path.splitext(input_file)[0]
        output_file = base + '.INPORB'

    print(f"Input:  {input_file}")
    print(f"Output: {output_file}")

    atoms = parse_gto_shells(input_file)
    print(f"Parsed {len(atoms)} atoms with shell structure:")
    for iatom, shells in enumerate(atoms):
        for sh in shells:
            cart = "Cartesian" if sh['cart'] else "Spherical"
            print(f"  atom {iatom + 1}: l={sh['l']}, nprim={sh['nprim']}, "
                  f"ncontr={sh['ncontr']}, {cart} ({sh['ncomp']} comp)")

    idx_map = build_index_map(atoms)
    nbasis, orbitals = parse_molden_mo(input_file, idx_map)

    d_pure_positions = compute_d_pure_positions(atoms)
    sqrt3 = math.sqrt(3.0)
    for orb in orbitals:
        for i in d_pure_positions:
            orb['coeffs'][i] /= sqrt3

    print(f"Total basis functions: {nbasis}")
    print(f"Total orbitals: {len(orbitals)}")
    write_gssorb(output_file, nbasis, orbitals)
    print("Done.")


if __name__ == '__main__':
    main()
