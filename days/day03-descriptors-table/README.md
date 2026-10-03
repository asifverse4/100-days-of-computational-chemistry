# Day 3: Molecular descriptors table

**What:** Turns a list of SMILES into a CSV of common descriptors: formula, molecular weight, logP, TPSA, H-bond donors and acceptors, rotatable bonds, rings, QED and Lipinski rule-of-five violations.
**Why:** A descriptor table is the first screening step before docking, QSAR or choosing which molecules to run through xTB/DFT. It also gives you a ready table for papers and theses.

## Install
```bash
pip install rdkit
```

## Usage
Single molecule (prints a table):
```bash
python tool.py "CC(=O)Oc1ccccc1C(=O)O"
```
Batch to CSV (one SMILES per line, optional name after a space, `#` for comments):
```bash
python tool.py example/input.smi -o example/descriptors.csv
```

## Columns
| Column | Meaning |
|---|---|
| mw | Molecular weight (g/mol) |
| logp | Crippen logP (lipophilicity estimate) |
| tpsa | Topological polar surface area (Å²) |
| hbd / hba | H-bond donors / acceptors (Lipinski definitions) |
| rot_bonds | Rotatable bonds |
| qed | Quantitative estimate of drug-likeness (0 to 1) |
| lipinski_violations | Count of: MW > 500, logP > 5, HBD > 5, HBA > 10 |

## Test
```bash
pip install pytest
pytest days/day03-descriptors-table
```

## Limitations
- logP is a calculated estimate, not a measured value.
- Lipinski's rules suit oral drug-likeness only. They say nothing about materials such as carbon dots.
- Invalid SMILES are skipped with a message and the exit code is 1.

## Next
Day 4: batch xTB optimization wrapper.
