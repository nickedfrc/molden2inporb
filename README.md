# molden2inporb
A simple python-based script that do translation from .molden file to OpenMolcas .InpOrb file, using Cartesian basis functions.
# molden2gssorb.py

Convert Molden orbital files to Molcas INPORB format.

## Features

1. **Basis function reordering**: Converts Molden ordering (contracted function → component)
   to SEWARD ordering (component → contracted function)
2. **Cartesian component reordering**: Corrects differences between the Molden standard
   component order and Molcas's internal `lx > ly > lz` descending convention
3. **Pure-d normalization compensation**: Divides the three pure d components
   (dxx, dyy, dzz) by √3 to compensate for differing Cartesian Gaussian
   normalization conventions between the Molden source and Molcas

## Usage

```bash
python molden2gssorb.py input.molden
python molden2gssorb.py input.molden -o output.INPORB
```

Default output filename: `<input_basename>.INPORB`

## SEWARD ordering conventions

SEWARD uses **component-first** ordering (`atom → shell → component → contracted function`),
in contrast to Molden's **contracted-function-first** ordering (`atom → shell → contracted function → component`).

All contracted functions sharing the same angular component are adjacent.
For example, an O atom with contracted basis (3s,2p,1d):

| | Order |
|---|------|
| SEWARD | s1, s2, s3, **px1, px2**, py1, py2, pz1, pz2, dxx, dxy, dxz, dyy, dyz, dzz |
| Molden | s1, s2, s3, px1, py1, pz1, px2, py2, pz2, dxx, dyy, dzz, dxy, dxz, dyz |

Cartesian component order within each shell follows `lx > ly > lz` descending:

| shell | component order |
|:---:|------|
| s | s |
| p | px, py, pz |
| d | dxx, dxy, dxz, dyy, dyz, dzz |
| f | fxxx, fxxy, fxxz, fxyy, fxyz, fxzz, fyyy, fyyz, fyzz, fzzz |

## Pure-d normalization compensation

For Cartesian d shells, the pure components (dxx, dyy, dzz) have a normalization
factor of 1/√3 arising from the double factorial (2lx−1)!!. Some programs embed
this factor in the basis function definition at the integral stage (e.g. Molcas),
while others carry it in the MO coefficients (common in Molden exports). This
script divides the pure d components by √3 so that the resulting INPORB is
consistent with Molcas's expectation, regardless of the Molden source.

The correction is applied automatically to all Cartesian d shells.

## Output format

- Format: INPORB 2.2
- 5 coefficients per line (`ES21.14`)
- Sections: `#INFO`, `#ORB`, `#OCC`, `#OCHR`, `#ONE`, `#INDEX`
