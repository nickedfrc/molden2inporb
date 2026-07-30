# molden2inporb
A simple python-based script that do translation from .molden file to OpenMolcas .InpOrb file, using Cartesian coordinates.
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

## Supported basis function types

| l | Type | Cartesian components | Spherical components | Notes |
|:---:|:---:|:---:|:---:|------|
| 0 | s | 1 | 1 | |
| 1 | p | 3 | — | |
| 2 | d | 6 | 5 | `[5D]` marker controls spherical |
| 3 | f | 10 | 7 | `[7F]` marker controls spherical |
| ≥4 | g/h/i | 15/21/28 | 9/11/13 | `[9G]` etc. markers |

## Pure-d normalization compensation

For Cartesian d shells, the pure components (dxx, dyy, dzz) have a normalization
factor of 1/√3 arising from the double factorial (2lx−1)!!. Some programs embed
this factor in the basis function definition at the integral stage (e.g. Molcas),
while others carry it in the MO coefficients (common in Molden exports). This
script divides the pure d components by √3 so that the resulting INPORB is
consistent with Molcas's expectation, regardless of the Molden source.

The correction is applied automatically to all Cartesian d shells.

## SEWARD ordering conventions

- Global: atom → shell → component → contracted function
- Cartesian components: `lx > ly > lz` descending priority
- d-shell order: dxx, dxy, dxz, dyy, dyz, dzz
- f-shell order: fxxx, fxxy, fxxz, fxyy, fxyz, fxzz, fyyy, fyyz, fyzz, fzzz

## Output format

- Format: INPORB 2.2
- 5 coefficients per line (`ES21.14`)
- Sections: `#INFO`, `#ORB`, `#OCC`, `#OCHR`, `#ONE`, `#INDEX`
