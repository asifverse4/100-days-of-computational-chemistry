# Day 7: Unit converter

**What:** Converts between hartree (Eh), eV, kcal/mol, kJ/mol, cm-1, nm, THz and kelvin, for one number or for a whole column of a CSV file.
**Why:** QM output is in hartree, papers quote kcal/mol or kJ/mol, spectra are in nm or cm-1, and Tauc plots need eV. Doing these by hand is where factor errors creep in.

## Usage
One value, several targets:
```bash
python tool.py 1 eh ev kcal/mol kj/mol cm-1 nm
```
```
1 Eh = 27.2114 eV
1 Eh = 627.509 kcal/mol
1 Eh = 2625.5 kJ/mol
1 Eh = 219475 cm-1
1 Eh = 45.5634 nm
```
An absorption edge to a band gap energy:
```bash
python tool.py 450 nm ev        # 450 nm = 2.7552 eV
```
No target unit shows every unit, and `-d` sets the significant digits:
```bash
python tool.py 2.5 ev -d 8
```
Convert a CSV column (for example `total_energy_eh` from Day 5's summary) and add the new columns:
```bash
python tool.py --csv example/energies.csv --column energy_eh --from eh --to ev kcal/mol
```
```
name,energy_eh,energy_eh_ev,energy_eh_kcal_mol
h_atom,-0.5,-13.60569312,-313.754737
he_atom,-2.9037,-79.01370224,-1822.09926
h2_molecule,-1.1745,-31.95977315,-737.0098773
```
Add `-o out.csv` to write a file. Empty cells stay empty. Unit names are forgiving: `eV`, `hartree`, `kcal mol-1`, `cm^-1`, `1/cm`, `Kelvin` all work. `--list` shows the supported units.

## Method
Every value is converted through joules per molecule using CODATA 2018 constants:

- 1 eV = 1.602176634e-19 J, and h, c, N_A and k_B are exact SI values.
- 1 Eh = 4.3597447222071e-18 J.
- 1 kcal = 4184 J (thermochemical calorie).
- Wavenumber: E = h c (100 x) for x in cm-1. Wavelength: E = h c / x for x in nm. Frequency: E = h x. Temperature: E = k_B x.

The tests check the results against reference values (for example 1 Eh = 27.211386 eV = 627.50947 kcal/mol = 219474.63 cm-1, and 1 eV = 1239.842 nm) and that converting there and back returns the input for every pair of units.

## Test
```bash
pytest days/day07-unit-converter
```

## Limitations
- These are per-molecule energy equivalents. Wavelength, wavenumber and temperature are not energies in themselves, they are treated as photon or thermal energy equivalents.
- nm needs a positive value, so a zero or negative energy shows `n/a` for nm. Negative cm-1 is allowed, so imaginary frequencies still convert.
- Vacuum wavelengths only. Air wavelengths differ slightly.
- Output is rounded to 6 significant digits by default. Use `-d` for more.
- The example energies in `example/energies.csv` are rounded non-relativistic reference values for H, He and H2.

## Next
Day 8: periodic-table property lookup CLI.
