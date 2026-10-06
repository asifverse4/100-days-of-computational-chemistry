# Day 6: XYZ file toolkit

**What:** Six small commands for structure files: `info`, `split`, `merge`, `center`, `align` and `convert` (xyz, pdb and Turbomole `coord`). No dependencies beyond Python.
**Why:** Between RDKit, xTB, ORCA and a viewer you constantly need to split a multi-frame trajectory, glue conformers together, move a molecule to the origin or hand xTB a `coord` file. This replaces the throwaway scripts.

## Usage
Inspect a structure (atom count, Hill formula, centroid, center of mass, size):
```bash
python tool.py info example/ethanol.xyz
```
```
file: example/ethanol.xyz
frames: 1
[frame 1] atoms: 9  formula: C2H6O
  centroid (A):           0.0000     0.0000     0.0000
  center of mass (A):     0.3648     0.0228     0.1756
  bounding box (A):       3.1143     2.6109     2.0574
  comment: ethanol | CCO | MMFF94 energy = -1.337 kcal/mol
```
Merge files into one multi-structure file, then split it back:
```bash
python tool.py merge example/ethanol.xyz example/benzene.xyz -o all.xyz
python tool.py split all.xyz --outdir frames          # frames/all_001.xyz, all_002.xyz
```
Center on the centroid (default) or the center of mass:
```bash
python tool.py center example/ethanol.xyz --mode mass -o ethanol_com.xyz
```
Convert by file extension (`.xyz`, `.pdb`, `.coord`):
```bash
python tool.py convert example/ethanol.xyz -o ethanol.pdb
python tool.py convert example/ethanol.xyz -o coord      # for xtb: xtb coord --opt
```
Align a mobile structure onto a reference with a Kabsch fit:
```bash
python tool.py align example/ethanol.xyz shifted.xyz -o aligned.xyz --rmsd
python tool.py align ref_traj.xyz mob_traj.xyz --reference-frame 2 --mobile-frame 2 -o frame2_aligned.pdb
```

| Command | Options |
|---|---|
| `split` | `--outdir`, `--prefix`, `--format {xyz,pdb,coord}` |
| `merge` | `-o` (required), `--strict` (require identical atoms in every frame, for trajectories) |
| `center` | `-o` (default `<input>_centered.<ext>`), `--mode {centroid,mass}`; every structure in the file is centered separately |
| `align` | `-o` (default `<mobile>_aligned.<ext>`), `--reference-frame` and `--mobile-frame` (1-based), `--rmsd`; aligns one selected mobile frame onto one selected reference frame |
| `convert` | `-o` (required) |

Errors are reported with a line number (for example `line 3: expected 'symbol x y z'`) and the exit code is 2.

## Method
- **Centroid:** the mean of the atom positions. **Center of mass:** the mass-weighted mean, using standard atomic weights for H to Kr plus common heavier elements (Ag, Au, Pt, Pb and others). An unknown element gives a clear error.
- **Kabsch alignment:** both frames are centered on their centroids, a covariance matrix is built, and a 3x3 SVD-based rotation is solved with reflection correction (proper rotation, determinant +1). The fitted mobile coordinates are then shifted onto the reference centroid.
- **Formula:** Hill order (C, then H, then the rest alphabetically; plain alphabetical if there is no carbon).
- **Units:** `.xyz` and `.pdb` are in angstrom. Turbomole `coord` is in bohr (1 bohr = 0.529177 A, CODATA 2018) and is converted on reading and writing.
- **PDB:** written as `HETATM` records (residue `MOL`, chain `A`) with the element in columns 77-78. Several structures become `MODEL`/`ENDMDL` blocks. When reading, the element column is used, or the atom name if it is missing.
- **Alignment compatibility checks:** `align` requires identical atom counts and element ordering between the selected frames.
- **Multi-frame alignment:** `align` processes one frame pair at a time (`--reference-frame`, `--mobile-frame`, default 1). It writes one aligned structure.
- The examples in `example/` are the RDKit-built geometries from Day 2.

## Test
```bash
pytest days/day06-xyz-toolkit
```
The tests cover round trips (xyz to pdb to xyz, xyz to coord to xyz), PDB column positions, bad-input messages, Kabsch alignment behavior (translation/rotation recovery, reflection handling, compatibility errors) and each command.

## Limitations
- No bonds or connectivity. PDB output has no `CONECT` records, so viewers guess the bonds from distances.
- PDB writing is for viewing and docking inputs. Residue, chain and atom-name information is not kept.
- `coord` holds one structure and ignores periodic-cell and constraint blocks.
- Only the first three frames are summarized by `info`.
- Alignment is rigid-body only: no atom reordering, no atom mapping/permutation search and no scaling.

## Next
Day 7: unit converter (Hartree, eV, kcal/mol, kJ/mol, nm, cm-1).
