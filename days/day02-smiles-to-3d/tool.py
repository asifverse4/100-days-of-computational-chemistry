"""Day 2: SMILES to 3D.

Convert a SMILES string (or a file with one SMILES per line) into optimized
3D geometries written as .xyz files, using RDKit (ETKDG embedding + MMFF94/UFF).

Usage:
    python tool.py "CCO" -o ethanol.xyz
    python tool.py example/input.smi --outdir example/output
"""
import argparse
import sys
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem


def smiles_to_mol(smiles: str, seed: int = 42, max_iters: int = 2000):
    """Return (mol with 3D coords, force field used, energy in kcal/mol)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles!r}")
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=seed) != 0:
        # fall back to random coordinates for difficult molecules
        if AllChem.EmbedMolecule(mol, randomSeed=seed, useRandomCoords=True) != 0:
            raise RuntimeError(f"3D embedding failed for {smiles!r}")
    if AllChem.MMFFHasAllMoleculeParams(mol):
        AllChem.MMFFOptimizeMolecule(mol, maxIters=max_iters)
        props = AllChem.MMFFGetMoleculeProperties(mol)
        energy = AllChem.MMFFGetMoleculeForceField(mol, props).CalcEnergy()
        return mol, "MMFF94", energy
    AllChem.UFFOptimizeMolecule(mol, maxIters=max_iters)
    energy = AllChem.UFFGetMoleculeForceField(mol).CalcEnergy()
    return mol, "UFF", energy


def write_xyz(mol, path: Path, comment: str) -> None:
    conf = mol.GetConformer()
    lines = [str(mol.GetNumAtoms()), comment]
    for atom in mol.GetAtoms():
        p = conf.GetAtomPosition(atom.GetIdx())
        lines.append(f"{atom.GetSymbol():<2} {p.x:12.6f} {p.y:12.6f} {p.z:12.6f}")
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", help="SMILES string or path to a .smi file (one SMILES per line, optional name after a space)")
    ap.add_argument("-o", "--output", help="output .xyz (single SMILES only)")
    ap.add_argument("--outdir", default=".", help="output directory for batch mode")
    ap.add_argument("--seed", type=int, default=42, help="random seed for embedding")
    args = ap.parse_args()

    src = Path(args.input)
    if src.is_file():
        entries = []
        for i, line in enumerate(src.read_text().splitlines(), 1):
            parts = line.split()
            if parts and not parts[0].startswith("#"):
                entries.append((parts[0], parts[1] if len(parts) > 1 else f"mol{i}"))
    else:
        entries = [(args.input, "molecule")]

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    failed = 0
    for smi, name in entries:
        try:
            mol, ff, e = smiles_to_mol(smi, seed=args.seed)
        except (ValueError, RuntimeError) as err:
            print(f"[skip] {name}: {err}", file=sys.stderr)
            failed += 1
            continue
        out = Path(args.output) if (args.output and len(entries) == 1) else outdir / f"{name}.xyz"
        write_xyz(mol, out, f"{name} | {smi} | {ff} energy = {e:.3f} kcal/mol")
        print(f"[ok] {name}: {mol.GetNumAtoms()} atoms, {ff} {e:.2f} kcal/mol -> {out}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
