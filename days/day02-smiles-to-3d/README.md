# Day 2: SMILES to 3D

**What:** Converts a SMILES string, or a file of SMILES, into optimized 3D geometries saved as `.xyz` files.
**Why:** Every xTB, ORCA or Psi4 calculation starts from a sensible 3D structure. This is the quick first step before any QM or docking work.

## Install
```bash
pip install rdkit
```

## Usage
Single molecule:
```bash
python tool.py "CCO" -o ethanol.xyz
```
Batch (one SMILES per line, optional name after a space, `#` for comments):
```bash
python tool.py example/input.smi --outdir example/output
```
Example output:
```
[ok] ethanol: 9 atoms, MMFF94 -1.34 kcal/mol -> example/output/ethanol.xyz
[ok] aspirin: 21 atoms, MMFF94 18.91 kcal/mol -> example/output/aspirin.xyz
```

## Method
1. Parse SMILES and add explicit hydrogens.
2. Embed in 3D with RDKit's ETKDG distance-geometry method (fixed seed for reproducibility).
3. Optimize with MMFF94, falling back to UFF if MMFF parameters are missing.
4. Write standard XYZ (atom count, comment line with name, SMILES, force field and energy, then coordinates in angstroms).

## Limitations
- Gives **one** conformer, not the global minimum. Use CREST or an RDKit conformer search for flexible molecules (planned for Day 25).
- Force-field energies are only comparable between conformers of the **same** molecule, not between different molecules.
- Metals and unusual bonding may fail to embed or fall back to UFF. Always inspect the geometry before a QM run.
- Not a substitute for a QM optimization. Use the output as the starting point for xTB or DFT.

## Next
Day 3: molecular descriptors table from a SMILES list.
