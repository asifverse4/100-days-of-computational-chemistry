"""Day 4: batch xTB optimization wrapper.

Run GFN-xTB geometry optimizations on many molecules and collect the results
(final energy, HOMO-LUMO gap, convergence status) into one CSV.

Input can be a .smi file (needs RDKit), a folder of .xyz files, or one .xyz file.

Usage:
    python tool.py example/input.smi --outdir example/output
    python tool.py geometries/ --charge -1 --solvent water
    python tool.py mol.xyz --gfn 2 --opt-level tight
"""
import argparse
import csv
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

COLUMNS = [
    "name", "status", "energy_eh", "gap_ev", "grad_norm", "iterations",
    "charge", "uhf", "runtime_s", "opt_xyz",
]

ENERGY_RE = re.compile(r"TOTAL ENERGY\s+(-?\d+\.\d+)\s+Eh", re.IGNORECASE)
GAP_RE = re.compile(r"HOMO-LUMO GAP\s+(-?\d+\.\d+)\s+eV", re.IGNORECASE)
GRAD_RE = re.compile(r"GRADIENT NORM\s+(-?\d+\.\d+)\s+Eh", re.IGNORECASE)
ITER_RE = re.compile(r"CONVERGED AFTER\s+(\d+)\s+ITERATIONS", re.IGNORECASE)


def parse_xtb_output(text: str) -> dict:
    """Extract the final energy, gap, gradient norm and convergence status."""

    def last(rx):
        found = rx.findall(text)
        return float(found[-1]) if found else None

    upper = text.upper()
    energy = last(ENERGY_RE)
    if "GEOMETRY OPTIMIZATION CONVERGED" in upper:
        status = "converged"
    elif "FAILED TO CONVERGE" in upper:
        status = "not_converged"
    elif energy is not None:
        status = "unknown"  # energy printed but no convergence message
    else:
        status = "failed"
    iters = ITER_RE.findall(text)
    return {
        "status": status,
        "energy_eh": energy,
        "gap_ev": last(GAP_RE),
        "grad_norm": last(GRAD_RE),
        "iterations": int(iters[-1]) if iters else None,
    }


def clean_xyz(text: str, name: str) -> str:
    """Rewrite an XYZ block with a plain comment line (xtb can trip on odd comments)."""
    lines = text.strip().splitlines()
    try:
        n = int(lines[0].split()[0])
    except (ValueError, IndexError):
        raise ValueError("not an XYZ file (first line must be the atom count)") from None
    atoms = [ln for ln in lines[2:] if ln.strip()]
    if len(atoms) != n:
        raise ValueError(f"XYZ says {n} atoms but found {len(atoms)} coordinate lines")
    return f"{n}\n{name}\n" + "\n".join(atoms) + "\n"


def smiles_to_xyz(smiles: str, name: str) -> tuple[str, int, int]:
    """Return (xyz text, formal charge, unpaired electrons) from a SMILES string."""
    from rdkit import Chem  # imported here so .xyz input works without RDKit
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles!r}")
    charge = Chem.GetFormalCharge(mol)
    uhf = sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms())
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=42) != 0 and AllChem.EmbedMolecule(
        mol, randomSeed=42, useRandomCoords=True
    ) != 0:
        raise ValueError(f"3D embedding failed for {smiles!r}")
    if AllChem.MMFFHasAllMoleculeParams(mol):
        AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
    else:
        AllChem.UFFOptimizeMolecule(mol, maxIters=2000)
    conf = mol.GetConformer()
    rows = []
    for atom in mol.GetAtoms():
        p = conf.GetAtomPosition(atom.GetIdx())
        rows.append(f"{atom.GetSymbol():<2} {p.x:12.6f} {p.y:12.6f} {p.z:12.6f}")
    return f"{mol.GetNumAtoms()}\n{name}\n" + "\n".join(rows) + "\n", charge, uhf


def safe_name(name: str) -> str:
    return re.sub(r"[^\w.-]", "_", name)


def unique(name: str, seen: set) -> str:
    base, k = safe_name(name), 2
    out = base
    while out in seen:
        out, k = f"{base}_{k}", k + 1
    seen.add(out)
    return out


def collect_entries(arg: str, charge: int = 0, uhf: int = 0) -> list[dict]:
    """Return a list of {name, xyz, charge, uhf, error} dicts from a file or folder."""
    src, seen, entries = Path(arg), set(), []
    if src.is_dir():
        files = sorted(src.glob("*.xyz"))
        if not files:
            raise ValueError(f"no .xyz files found in {src}")
    elif src.is_file():
        files = [src]
    else:
        raise ValueError(f"input not found: {arg}")

    if files[0].suffix.lower() == ".smi":
        for i, line in enumerate(files[0].read_text().splitlines(), 1):
            parts = line.split()
            if not parts or parts[0].startswith("#"):
                continue
            name = unique(parts[1] if len(parts) > 1 else f"mol{i}", seen)
            try:
                xyz, chg, unp = smiles_to_xyz(parts[0], name)
                entries.append({"name": name, "xyz": xyz, "charge": chg, "uhf": unp, "error": None})
            except ValueError as err:
                entries.append({"name": name, "xyz": None, "charge": 0, "uhf": 0, "error": str(err)})
        return entries

    for f in files:
        name = unique(f.stem, seen)
        try:
            xyz = clean_xyz(f.read_text(), name)
            entries.append({"name": name, "xyz": xyz, "charge": charge, "uhf": uhf, "error": None})
        except ValueError as err:
            entries.append({"name": name, "xyz": None, "charge": charge, "uhf": uhf, "error": str(err)})
    return entries


def run_xtb(entry: dict, outdir: Path, xtb: str = "xtb", gfn: int = 2,
            solvent: str | None = None, opt_level: str = "normal",
            timeout: int = 3600) -> dict:
    """Optimize one molecule with xtb and return a result row."""
    wd = outdir / entry["name"]
    wd.mkdir(parents=True, exist_ok=True)
    (wd / "input.xyz").write_text(entry["xyz"])
    cmd = [xtb, "input.xyz", "--opt", opt_level, "--gfn", str(gfn),
           "--chrg", str(entry["charge"]), "--uhf", str(entry["uhf"])]
    if solvent:
        cmd += ["--alpb", solvent]

    row = {"name": entry["name"], "charge": entry["charge"], "uhf": entry["uhf"], "opt_xyz": ""}
    start = time.time()
    try:
        proc = subprocess.run(cmd, cwd=wd, capture_output=True, text=True,
                              timeout=timeout, check=False)
        text = proc.stdout + "\n" + proc.stderr
        (wd / "xtb.out").write_text(text)
        row.update(parse_xtb_output(text))
        if proc.returncode != 0 and row["status"] != "converged":
            row["status"] = "failed"
    except subprocess.TimeoutExpired:
        row.update({"status": "timeout", "energy_eh": None, "gap_ev": None,
                    "grad_norm": None, "iterations": None})
    row["runtime_s"] = round(time.time() - start, 1)
    if (wd / "xtbopt.xyz").is_file():
        row["opt_xyz"] = str(wd / "xtbopt.xyz")
    return row


def fmt(value, digits):
    return "n/a" if value is None else f"{value:.{digits}f}"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", help=".smi file, folder of .xyz files, or a single .xyz file")
    ap.add_argument("--outdir", default="xtb_results", help="output directory (default: xtb_results)")
    ap.add_argument("--gfn", type=int, choices=[0, 1, 2], default=2, help="GFN-xTB level (default: 2)")
    ap.add_argument("--solvent", help="ALPB implicit solvent, e.g. water, ethanol, dmso")
    ap.add_argument("--opt-level", default="normal",
                    choices=["crude", "sloppy", "loose", "normal", "tight", "vtight"],
                    help="optimization convergence level (default: normal)")
    ap.add_argument("--charge", type=int, default=0, help="total charge for .xyz input (SMILES charge is automatic)")
    ap.add_argument("--uhf", type=int, default=0, help="unpaired electrons for .xyz input")
    ap.add_argument("--timeout", type=int, default=3600, help="seconds per molecule (default: 3600)")
    ap.add_argument("--xtb", default="xtb", help="path to the xtb executable")
    args = ap.parse_args()

    if shutil.which(args.xtb) is None:
        print(f"xtb executable not found: {args.xtb!r}. Install xtb (e.g. conda install -c conda-forge xtb) "
              "or pass --xtb /path/to/xtb", file=sys.stderr)
        return 2
    try:
        entries = collect_entries(args.input, args.charge, args.uhf)
    except ValueError as err:
        print(err, file=sys.stderr)
        return 2

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rows, bad = [], 0
    for e in entries:
        if e["error"]:
            print(f"[skip] {e['name']}: {e['error']}", file=sys.stderr)
            rows.append({"name": e["name"], "status": "invalid_input", "charge": e["charge"], "uhf": e["uhf"]})
            bad += 1
            continue
        row = run_xtb(e, outdir, args.xtb, args.gfn, args.solvent, args.opt_level, args.timeout)
        rows.append(row)
        ok = row["status"] == "converged"
        bad += not ok
        print(f"[{'ok' if ok else 'warn'}] {row['name']}: {row['status']}, "
              f"E = {fmt(row['energy_eh'], 6)} Eh, gap = {fmt(row['gap_ev'], 3)} eV ({row['runtime_s']} s)")

    csv_path = outdir / "results.csv"
    with csv_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS, lineterminator="\n", restval="")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} rows to {csv_path}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
