# molden2inporb

Convert standard Molden molecular orbitals to OpenMolcas **INPORB 2.2**.
Runs with **Python 3.6 or newer**, using only the standard library.

**Only the Molden file is required by default.** The script automatically detects
the shells, converts AO component ordering and normalization, and writes an
`.INPORB` file with the same basename. No manual atom list or `--same-atom-order`
option is needed. Atom order is preserved by default; if the target Molcas input
uses a different order, `--target-log` matches atoms automatically by element
and coordinates.

## Quick start

```bash
python molden2inporb.py input.molden
```

This writes `input.INPORB` beside the input file. The existing explicit-output
command continues to work:

```bash
python molden2inporb.py input.molden -o output.INPORB
```

If you already installed the script as `molden2inporb` on your `PATH`, use:

```bash
molden2inporb input.molden
```

No packages, basis-set names, normalization flags, or atom lists are required
for the usual case with matching source/target atom order.

## Compatibility with the original script

The previous version is preserved as a regression reference from commit
[`78cf0a8`](https://github.com/nickedfrc/molden2inporb/commit/78cf0a80e365256de23ed6d7d4a3e3b8f313fa55)
in `tests/fixtures/legacy_molden2inporb.py`.

For the previously correct **Cartesian s/p/6D cases**, with matching atom order
and the same ascending angular-momentum/contraction layout:

- Both original command forms and the default output filename are unchanged.
- Atom order is preserved by default. There is no automatic element sorting.
- The same component-first AO permutation is used.
- Cartesian `xx`, `yy`, `zz` coefficients still use **division by `sqrt(3)`**,
  preserving the original floating-point results exactly, rather than replacing
  division with multiplication by an approximate reciprocal.
- MO order, occupations, sparse-zero coefficients, and the default
  occupation-based `#INDEX` assignment are retained.
- The public helper names `parse_gto_shells`, `build_index_map`,
  `compute_d_pure_positions`, `parse_molden_mo`, `fmt_es`, `write_coeffs_block`,
  and `write_gssorb` remain available with their original positional arguments.
- `[Atoms]` remains unnecessary when no coordinate-based atom matching is requested.

Tests compare **48 synthetic contraction layouts** against the old implementation
using exact equality of every converted coefficient and occupation. They also
check unchanged `#ORB`, `#OCC`, `#OCHR`, and previously valid `#INDEX` output text.
This is a tested compatibility scope, not a claim that every historical input
or every quantum-chemistry exporter has been verified.

The entire output file is **not byte-identical**: incorrect metadata/formatting
is fixed, and the birth-certificate timestamp/process ID naturally varies.
The intentional corrections are:

- `#INFO` writes **NBAS then NORB**; the old order failed for rectangular MO sets.
- `#ONE` uses the actual INPORB 2.2 energy format, **10 fields per line, 12
  characters per field**, instead of the MO-coefficient format.
- Zero orbital energies remain zero; no artificial energy spectrum is generated.
- `#INDEX` row labels wrap from 9 back to 0 after 100 orbitals.
- Incorrect f/g and spherical transformations are fixed, as detailed below.
- Malformed/duplicate/out-of-range coefficients, unsupported shells, and beta
  orbitals produce errors instead of silently writing a misleading file.

## Different atom orders: automatic matching

INPORB coefficients have no atom identities. Their rows must match the target
SEWARD AO order. If the target geometry lists atoms differently, provide its log:

```bash
python molden2inporb.py input.molden --target-log target.log
```

This matches atoms by **element and Cartesian coordinates**, checks the total
number of AOs and C1 symmetry, and prints the selected source atom IDs. You do
not need to construct a list yourself. Coordinates must have the same origin
and orientation; coordinate matching tolerates rounding in the printed log.
The first coordinate table is used, so use a single-geometry calculation log.

For advanced/manual use, `--atom-order 1,3,2` means that target atoms 1, 2, 3
correspond to Molden atom IDs 1, 3, 2. It is optional and mutually exclusive with
`--target-log` and `--same-atom-order`.

`--same-atom-order` explicitly selects the default and is accepted for
compatibility with the earlier development version; it is never required.

**A Molden file alone cannot reveal a different target input's atom order.**
The no-option command therefore preserves the input order, which protects
existing workflows. Atom matching checks coordinates and AO count, not the
complete target basis: matching radial contractions and angular representation
remain necessary.

## Basis and AO conventions

The conversion depends on the actual `[GTO]` shells and angular markers, **not
on whether the basis is named split-valence, cc-pVDZ, cc-pVTZ, or cc-pVQZ**.

### Ordering

Molden uses `atom -> contracted shell -> angular component`.
SEWARD uses `atom -> angular momentum -> component -> contracted function`.
The script groups contractions of each angular momentum on each atom; interleaved
`sp` shells are supported and expanded into their s and p contractions.

For a `(3s,2p,1d)` Cartesian atom:

| Format | AO order |
|---|---|
| Molden | s1, s2, s3, px1, py1, pz1, px2, py2, pz2, xx, yy, zz, xy, xz, yz |
| Molcas | s1, s2, s3, px1, px2, py1, py2, pz1, pz2, xx, xy, xz, yy, yz, zz |

Molcas Cartesian powers are ordered by decreasing x exponent, then decreasing
y exponent. The f order is:

```text
xxx, xxy, xxz, xyy, xyz, xzz, yyy, yyz, yzz, zzz
```

Its source indices in the standard Molden 10F order, starting at zero, are:

```text
0, 4, 5, 3, 9, 6, 1, 8, 7, 2
```

The old f permutation was incorrect, and the old g fallback did no reordering.

### Cartesian normalization

For a component with powers `(a,b,c)`, standard Molden coefficients transform as

```text
C_Molcas = C_Molden / sqrt[(2a-1)!! (2b-1)!! (2c-1)!!]
```

Here `(-1)!! = 1`. The factors extend the original correct 6D conversion:

| Shell | Component family | Divide coefficient by |
|---|---|---|
| s, p | all | 1 |
| d | xx, yy, zz | sqrt(3) |
| d | xy, xz, yz | 1 |
| f | xxx and permutations | sqrt(15) |
| f | xxy and permutations | sqrt(3) |
| f | xyz | 1 |
| g | xxxx and permutations | sqrt(105) |
| g | xxxy and permutations | sqrt(15) |
| g | xxyy and permutations | 3 |
| g | xxyz and permutations | sqrt(3) |

This converts AO conventions; it does not renormalize whole MO vectors.
The appropriate numerical check is `C.T @ S @ C = I`, using the AO overlap
matrix, not `C.T @ C = I`. Nonstandard exporter conventions must be checked
separately; the script cannot infer arbitrary source normalization conventions.

### Spherical shells

Standard real spherical Molden components are ordered
`0,+1,-1,+2,-2,...`; Molcas d/f/g components are ordered `-l,...,0,...,+l`.
They are reordered without the Cartesian normalization factors.

Supported markers follow the Molden specification:

| Marker | Meaning |
|---|---|
| none | Cartesian 6D/10F/15G |
| `[5D]` or `[5D7F]` | spherical 5D **and** 7F |
| `[5D10F]` | spherical 5D, Cartesian 10F |
| `[7F]` | spherical 7F; d remains Cartesian unless separately specified |
| `[9G]` | spherical 9G |

The old script treated spherical d as Cartesian during normalization, and did
not implement these spherical permutations correctly.

## Scope and limitations

- **Python 3.6+**, standard library only. No `math.prod` or zero-width `re.split`.
- Target **C1**, restricted orbitals; beta/UHF orbitals are rejected.
- Supported shells: **s, p, d, f, g**, and combined **sp**.
- Source and target must have the same AO space, radial contractions, contraction
  order, geometry orientation, and Cartesian/spherical representation.
- **No basis projection or 5D-to-6D conversion.** Changing representation changes
  the AO dimension and requires more than coefficient permutation/scaling.
- h/i and higher shells are rejected; do not assume every element or every
  member of the cc-pVXZ family falls within the supported angular range.
- Nonstandard normalization, special ECP/ghost conventions, symmetry-adapted
  target orbitals, and general-contraction multi-column GTO extensions are not
  automatically handled.
- Sparse MO indices retain the original interpretation: omitted entries are
  zero. The converter cannot distinguish intentionally omitted zeros from an
  accidentally omitted interior coefficient.

The default `#INDEX` still marks `Occup > 0.5` as inactive (`i`) and the rest as
secondary (`s`) to preserve original behavior. Natural occupations do **not**
uniquely define a CASSCF active space. To omit this heuristic:

```bash
python molden2inporb.py input.molden --no-index
```

Specify the intended `INACTIVE`, `RAS2`, and `NACTEL` in the target RASSCF input.
An orbital swap appropriate for SCF orbitals should not automatically be reused
after importing a differently ordered set of CASSCF orbitals.

## Verification

```bash
python -W error -m unittest discover -s tests -v
```

The offline suite covers legacy numerical compatibility and output sections,
default/explicit CLI usage, Python 3.6 MO parsing, f/g permutations and factors,
spherical markers, sp shells, automatic atom matching, sparse input, malformed
input, rectangular dimensions, fixed-width energy fields, and index labels.
GitHub Actions runs it on Windows/Linux with current Python and in a Python 3.6
container. All fixtures are synthetic; private molecular input/log files are not
part of this repository.

The revised converter was additionally checked locally on a 140-AO Cartesian
cc-pVDZ example: independently calculated `max|C.T @ S @ C - I|` was about
`2e-10`, and the converted first 35 orbitals matched the target log's printed
RASSCF coefficients within its four-decimal rounding after allowing a sign
change for each MO. This is a representation/format check, **not an OpenMolcas
rerun**. Real high-angular-momentum calculations have not been rerun here.

## Sources

- [Molden format specification](https://www.theochem.ru.nl/molden/molden_format.html)
- [OpenMolcas Molden export and normalization](https://github.com/Molcas/OpenMolcas/blob/master/src/property_util/molden_interface.F90)
- [OpenMolcas Cartesian and spherical labels](https://github.com/Molcas/OpenMolcas/blob/master/src/integral_util/make_labels.F90)
- [INPORB formatting parameters](https://github.com/Molcas/OpenMolcas/blob/master/src/property_util/inporbfmt.F90)
- [OpenMolcas INPORB writer](https://github.com/Molcas/OpenMolcas/blob/master/src/property_util/wrvec_.F90)
