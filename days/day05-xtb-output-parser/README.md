# Day 5: xTB output parser

**What:** Reads xtb output files and collects the total energy, HOMO-LUMO gap, HOMO/LUMO energies, energy components and atomic charges into CSV tables.
**Why:** After a batch of xTB runs (Day 4) the numbers are scattered across dozens of log files. This turns them into one table you can sort, plot or paste into a thesis, and ranks conformers by relative energy.

## Usage
One file:
```bash
python tool.py path/to/xtb.out
```
A whole Day 4 results folder (scans subfolders for `*.out` and `*.log`), with per-atom charges:
```bash
python tool.py ../day04-xtb-batch/xtb_results -o summary.csv --charges-csv charges.csv
```
Rank conformers or isomers of one molecule by relative energy:
```bash
python tool.py conformers/ --relative -o ranked.csv
```
A tiny **synthetic** demo is included so you can try the command right away:
```bash
python tool.py example/synthetic
```
The numbers in `example/synthetic` are made up. They only show the file format. Use real xtb output for real work.

## Output columns
| Column | Meaning |
|---|---|
| status | `converged`, `not_converged`, `singlepoint` or `failed` |
| method, xtb_version | e.g. `GFN2-xTB`, `6.7.1` (empty if not found) |
| total_energy_eh | final total energy (hartree) |
| rel_energy_kcal | with `--relative`: energy above the lowest run (kcal/mol) |
| gradient_norm | final gradient norm (Eh/a0) |
| gap_ev, homo_ev, lumo_ev | HOMO-LUMO gap and orbital energies (eV) |
| scc_energy_eh, repulsion_eh, dispersion_eh | energy components from the summary block |
| total_charge | total charge of the system (e) |
| iterations | optimization cycles to convergence |
| charge_source | where the charges came from (see below) |

The charges file has one row per atom: `name, atom_index, element, charge, charge_source`.

## Method
1. Find xtb output files. Files that do not look like xtb output are ignored when scanning a folder.
2. Take the **last** value of each quantity from the summary block, so geometry optimizations report the final geometry.
3. Atomic charges come from the `charges` file xtb writes next to its output (final geometry). If it is missing, the charge table printed in the output is used instead, and `charge_source` says so (that table is for the starting geometry).
4. `--relative` converts hartree to kcal/mol (1 Eh = 627.5095 kcal/mol) relative to the lowest energy in the batch.

## Test
```bash
pytest days/day05-xtb-output-parser
```

## Limitations
- Written for the text format of xtb 6.x. The tests use excerpts in that format, not output from every xtb version. If a field comes back empty on a real run, open an issue and attach the `xtb.out`.
- Charges are Mulliken-type charges from xtb, which depend on the method. Do not compare them with DFT charges directly.
- `--relative` is only meaningful for isomers or conformers of one molecule. It warns when atom counts differ.
- The HOMO-LUMO gap is a GFN-xTB orbital gap, not an optical or DFT gap.
- Does not parse frequencies or thermochemistry yet.

## Next
Day 6: XYZ file toolkit (convert, merge, split, center).
