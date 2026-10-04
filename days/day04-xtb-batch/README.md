# Day 4: Batch xTB optimization

**What:** Runs GFN-xTB geometry optimizations on many molecules (a SMILES list or a folder of `.xyz` files) and collects the final energy, HOMO-LUMO gap and convergence status into one CSV.
**Why:** Before DFT, you usually pre-optimize a whole set of structures (ligands, carbon-dot fragments, catalysts) with xTB. Doing that by hand means running `xtb` dozens of times and copying numbers out of log files.

## Install
```bash
pip install rdkit            # only needed for .smi input
conda install -c conda-forge xtb
```
xTB is not on PyPI. On Windows, use WSL or conda. The tool stops with a clear message if `xtb` is not found.

## Usage
From a SMILES list (charge and radicals are read from the SMILES automatically):
```bash
python tool.py example/input.smi --outdir example/output
```
From a folder of `.xyz` files (set the charge yourself):
```bash
python tool.py geometries/ --charge -1 --solvent water --opt-level tight
```

| Option | Meaning |
|---|---|
| `--gfn {0,1,2}` | GFN-xTB level (default 2) |
| `--solvent NAME` | ALPB implicit solvent, e.g. `water`, `dmso` |
| `--opt-level` | crude, sloppy, loose, normal (default), tight, vtight |
| `--charge`, `--uhf` | total charge and unpaired electrons for `.xyz` input |
| `--timeout` | seconds allowed per molecule (default 3600) |
| `--xtb PATH` | path to the xtb executable |

## Output
Each molecule gets its own folder `<outdir>/<name>/` with `input.xyz`, the full `xtb.out` log and the optimized `xtbopt.xyz`. A summary goes to `<outdir>/results.csv`:

```
name,status,energy_eh,gap_ev,grad_norm,iterations,charge,uhf,runtime_s,opt_xyz
```
Status is one of `converged`, `not_converged`, `timeout`, `failed`, `unknown` (energy printed but no convergence message) or `invalid_input`. The exit code is 1 if any molecule is not `converged`.

Run it on `example/input.smi` and paste your own `results.csv` here, since values depend on your xtb version.

## Method
1. Build the 3D start geometry (RDKit ETKDG + MMFF94/UFF, same as Day 2) or read your `.xyz`.
2. Run `xtb input.xyz --opt <level> --gfn <n> --chrg <q> --uhf <u> [--alpb <solvent>]` in a separate folder per molecule.
3. Parse the **last** `TOTAL ENERGY`, `HOMO-LUMO GAP` and `GRADIENT NORM` lines, and the convergence message, from the xtb output.

Energies are in hartree (Eh) and gaps in eV. Reference: Bannwarth, Ehlert, Grimme, *J. Chem. Theory Comput.* 2019, 15, 1652 (GFN2-xTB).

## Test
```bash
pytest days/day04-xtb-batch
```
The tests use sample xtb output text and a mocked `subprocess`, so they run without xtb installed.

## Limitations
- Runs molecules one after another. No parallel jobs yet.
- Optimizes one conformer from a single start geometry. For flexible molecules, use a conformer search first (Day 25).
- The gap printed by xtb is a GFN-xTB orbital gap, not a DFT or optical gap.
- Unusual xtb versions may word their messages differently. If a good run shows `unknown`, check `xtb.out` and open an issue.

## Next
Day 5: xTB output parser for energies, gaps and atomic charges.
