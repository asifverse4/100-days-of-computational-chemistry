"""Day 3: molecular descriptors table.

Compute common descriptors (formula, MW, logP, TPSA, H-bond donors/acceptors,
rotatable bonds, rings, QED) and Lipinski rule-of-five violations for a list
of SMILES, and save the result as CSV.

Usage:
    python tool.py "CC(=O)Oc1ccccc1C(=O)O"
    python tool.py example/input.smi -o example/descriptors.csv
"""
import argparse
import csv
import sys
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import QED, Crippen, Descriptors, Lipinski, rdMolDescriptors

COLUMNS = [
    "name", "smiles", "formula", "mw", "logp", "tpsa", "hbd", "hba",
    "rot_bonds", "rings", "aromatic_rings", "heavy_atoms", "qed", "lipinski_violations",
]


def describe(smiles: str, name: str) -> dict:
    """Return a dict of descriptors for one SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles!r}")
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
    return {
        "name": name,
        "smiles": smiles,
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "mw": round(mw, 2),
        "logp": round(logp, 2) + 0.0,
        "tpsa": round(rdMolDescriptors.CalcTPSA(mol), 2),
        "hbd": hbd,
        "hba": hba,
        "rot_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "rings": rdMolDescriptors.CalcNumRings(mol),
        "aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "heavy_atoms": mol.GetNumHeavyAtoms(),
        "qed": round(QED.qed(mol), 3),
        "lipinski_violations": violations,
    }


def read_entries(arg: str) -> list[tuple[str, str]]:
    """Read (smiles, name) pairs from a .smi file, or treat arg as one SMILES."""
    src = Path(arg)
    if not src.is_file():
        return [(arg, "molecule")]
    entries = []
    for i, line in enumerate(src.read_text().splitlines(), 1):
        parts = line.split()
        if parts and not parts[0].startswith("#"):
            entries.append((parts[0], parts[1] if len(parts) > 1 else f"mol{i}"))
    return entries


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", help="SMILES string or .smi file (one SMILES per line, optional name)")
    ap.add_argument("-o", "--output", help="write CSV to this path (default: print table)")
    args = ap.parse_args()

    rows, failed = [], 0
    for smi, name in read_entries(args.input):
        try:
            rows.append(describe(smi, name))
        except ValueError as err:
            print(f"[skip] {name}: {err}", file=sys.stderr)
            failed += 1

    if args.output:
        out = Path(args.output)
        out.parent.mkdir(parents=True, exist_ok=True)
        with out.open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=COLUMNS)
            writer.writeheader()
            writer.writerows(rows)
        print(f"Wrote {len(rows)} molecules to {out}")
    else:
        print("\t".join(COLUMNS))
        for r in rows:
            print("\t".join(str(r[c]) for c in COLUMNS))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
