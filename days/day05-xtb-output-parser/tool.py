"""Day 5: xTB output parser.

Read xtb output files (``xtb.out``, ``*.out`` or ``*.log``) and collect energies,
HOMO-LUMO gaps, orbital energies and atomic charges into CSV tables.
Works directly on the folders written by Day 4 (one ``xtb.out`` per molecule).

Usage:
    python tool.py path/to/xtb.out
    python tool.py ../day04-xtb-batch/xtb_results -o summary.csv --charges-csv charges.csv
    python tool.py conformers/ --relative -o ranked.csv
"""
import argparse
import csv
import re
import sys
from pathlib import Path

EH_TO_KCAL = 627.5094740631  # hartree to kcal/mol
NUM = r"(-?\d+\.\d+)"
FLAGS = re.IGNORECASE

# Values from the final SUMMARY block ("::" lines) or the "|" energy block.
SUMMARY_PATTERNS = {
    "total_energy_eh": rf"(?:::|\|)\s*total energy\s+{NUM}\s+Eh",
    "gradient_norm": rf"(?:::|\|)\s*gradient norm\s+{NUM}\s+Eh",
    "gap_ev": rf"(?:::|\|)\s*HOMO-LUMO gap\s+{NUM}\s+eV",
    "scc_energy_eh": rf"::\s*SCC energy\s+{NUM}\s+Eh",
    "dispersion_eh": rf"::\s*->\s*dispersion\s+{NUM}\s+Eh",
    "repulsion_eh": rf"::\s*repulsion energy\s+{NUM}\s+Eh",
    "total_charge": rf"::\s*total charge\s+{NUM}\s+e",
}
SUMMARY_RES = {k: re.compile(v, FLAGS) for k, v in SUMMARY_PATTERNS.items()}
HOMO_RE = re.compile(rf"{NUM}\s+\(HOMO\)\s+{NUM}", FLAGS)
LUMO_RE = re.compile(rf"{NUM}\s+\(LUMO\)\s+{NUM}", FLAGS)
ITER_RE = re.compile(r"CONVERGED AFTER\s+(\d+)\s+ITERATIONS", FLAGS)
VERSION_RE = re.compile(r"xtb version\s+([0-9][\w.]*)", FLAGS)
METHOD_RE = re.compile(r"Hamiltonian\s+(GFN\S*)", FLAGS)
NATOMS_RE = re.compile(r"(?:number of atoms|# atoms)\s+(\d+)", FLAGS)
TABLE_HEADER_RE = re.compile(r"^\s*#\s+Z\s+covCN\s+q\b")
TABLE_ROW_RE = re.compile(
    rf"^\s*(\d+)\s+(\d+)\s+([A-Za-z]{{1,2}})\s+{NUM}\s+{NUM}\s+{NUM}\s+{NUM}\s*$"
)

SUMMARY_COLUMNS = [
    "name", "status", "method", "xtb_version", "n_atoms", "total_energy_eh",
    "gradient_norm", "gap_ev", "homo_ev", "lumo_ev", "scc_energy_eh",
    "repulsion_eh", "dispersion_eh", "total_charge", "iterations",
    "charge_source", "source",
]
CHARGE_COLUMNS = ["name", "atom_index", "element", "charge", "charge_source"]
OK_STATUS = {"converged", "singlepoint"}


class NotXtbOutput(ValueError):
    """The file does not look like xtb output."""


def is_xtb_output(text: str) -> bool:
    return "xtb version" in text.lower() or "x T B" in text


def _last(rx, text: str, group: int = 1):
    found = rx.findall(text)
    if not found:
        return None
    hit = found[-1]
    value = hit[group - 1] if isinstance(hit, tuple) else hit
    return float(value)


def parse_charge_table(text: str) -> list[dict]:
    """Return the last 'Z covCN q C6AA alpha' table printed in the output."""
    table, current, reading = [], [], False
    for line in text.splitlines():
        if TABLE_HEADER_RE.match(line):
            reading, current = True, []
            continue
        if reading:
            m = TABLE_ROW_RE.match(line)
            if m:
                current.append({"index": int(m[1]), "element": m[3], "charge": float(m[5])})
            else:
                table, reading = current or table, False
    if reading and current:
        table = current
    return table


def parse_xtb_text(text: str) -> dict:
    """Extract energies, gaps, orbital energies and status from xtb output text."""
    info = {key: _last(rx, text) for key, rx in SUMMARY_RES.items()}
    info["homo_ev"] = _last(HOMO_RE, text, group=2)
    info["lumo_ev"] = _last(LUMO_RE, text, group=2)
    if info["gap_ev"] is None and None not in (info["homo_ev"], info["lumo_ev"]):
        info["gap_ev"] = round(info["lumo_ev"] - info["homo_ev"], 4)

    upper = text.upper()
    if "GEOMETRY OPTIMIZATION CONVERGED" in upper:
        info["status"] = "converged"
    elif "FAILED TO CONVERGE" in upper:
        info["status"] = "not_converged"
    elif info["total_energy_eh"] is not None and "ABNORMAL TERMINATION" not in upper:
        info["status"] = "singlepoint"
    else:
        info["status"] = "failed"

    iters = ITER_RE.findall(text)
    info["iterations"] = int(iters[-1]) if iters else None
    version, method, natoms = VERSION_RE.search(text), METHOD_RE.search(text), NATOMS_RE.search(text)
    info["xtb_version"] = version[1] if version else None
    info["method"] = method[1] if method else None
    info["n_atoms"] = int(natoms[1]) if natoms else None
    return info


def read_charges_file(path: Path) -> list[float]:
    """Read the plain 'charges' file xtb writes next to its output (one value per atom)."""
    return [float(tok) for tok in path.read_text().split()]


def read_xyz_elements(path: Path) -> list[str]:
    lines = path.read_text().splitlines()
    return [ln.split()[0] for ln in lines[2:] if ln.strip()]


def run_charges(folder: Path, table: list[dict]) -> tuple[list[dict], str]:
    """Pick the best available charges: the `charges` file (final geometry) or the output table."""
    cfile = folder / "charges"
    if cfile.is_file():
        try:
            values = read_charges_file(cfile)
        except ValueError:
            values = []
        elements: list[str] = []
        for xyz in ("xtbopt.xyz", "input.xyz"):
            if (folder / xyz).is_file():
                elements = read_xyz_elements(folder / xyz)
                break
        if not elements or len(elements) != len(values):
            same_size = len(table) == len(values)
            elements = [r["element"] for r in table] if same_size else ["X"] * len(values)
        if values:
            rows = [{"index": i, "element": el, "charge": q}
                    for i, (el, q) in enumerate(zip(elements, values), 1)]
            return rows, "charges_file (final geometry)"
    if table:
        return table, "output_table (initial geometry)"
    return [], ""


def run_name(path: Path) -> str:
    path = path.resolve()
    return path.parent.name if path.stem == "xtb" else path.stem


def parse_run(path: Path) -> tuple[dict, list[dict]]:
    """Parse one xtb output file. Returns (summary row, per-atom charge rows)."""
    text = path.read_text(errors="replace")
    if not is_xtb_output(text):
        raise NotXtbOutput(f"{path} does not look like xtb output")
    info = parse_xtb_text(text)
    charges, source = run_charges(path.resolve().parent, parse_charge_table(text))
    name = run_name(path)
    row = dict(info, name=name, source=str(path), charge_source=source)
    if row["n_atoms"] is None and charges:
        row["n_atoms"] = len(charges)
    atoms = [dict(c, name=name, atom_index=c["index"], charge_source=source) for c in charges]
    return row, atoms


def find_outputs(paths: list[str]) -> list[tuple[Path, bool]]:
    """Expand files and folders into (path, explicitly_given) pairs."""
    found, seen = [], set()
    for arg in paths:
        p = Path(arg)
        if p.is_dir():
            files = sorted(set(p.rglob("*.out")) | set(p.rglob("*.log")))
            items = [(f, False) for f in files]
        elif p.is_file():
            items = [(p, True)]
        else:
            raise ValueError(f"not found: {arg}")
        for f, explicit in items:
            key = f.resolve()
            if key not in seen:
                seen.add(key)
                found.append((f, explicit))
    return found


def add_relative_energies(rows: list[dict]) -> None:
    """Add rel_energy_kcal (relative to the lowest total energy in the batch)."""
    energies = [r["total_energy_eh"] for r in rows if r["total_energy_eh"] is not None]
    if not energies:
        return
    lowest = min(energies)
    for r in rows:
        e = r["total_energy_eh"]
        r["rel_energy_kcal"] = None if e is None else round((e - lowest) * EH_TO_KCAL, 3)
    sizes = {r["n_atoms"] for r in rows if r["n_atoms"] is not None}
    if len(sizes) > 1:
        print("[warn] --relative compares molecules with different atom counts; "
              "relative energies are only meaningful for isomers or conformers of one molecule.",
              file=sys.stderr)


def write_csv(path: Path, columns: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, lineterminator="\n",
                                restval="", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+", help="xtb output files and/or folders to scan (recursive)")
    ap.add_argument("-o", "--output", help="write the summary CSV here (default: print a table)")
    ap.add_argument("--charges-csv", help="also write per-atom charges (long format) to this CSV")
    ap.add_argument("--relative", action="store_true",
                    help="add rel_energy_kcal vs the lowest energy (isomers/conformers of ONE molecule)")
    args = ap.parse_args()

    try:
        outputs = find_outputs(args.inputs)
    except ValueError as err:
        print(err, file=sys.stderr)
        return 2

    rows, atoms, bad = [], [], 0
    for path, explicit in outputs:
        try:
            row, charges = parse_run(path)
        except NotXtbOutput as err:
            if explicit:
                print(f"[skip] {err}", file=sys.stderr)
                bad += 1
            continue
        rows.append(row)
        atoms.extend(charges)
        if row["status"] not in OK_STATUS:
            print(f"[warn] {row['name']}: {row['status']}", file=sys.stderr)
            bad += 1
    if not rows:
        print("no xtb output files found", file=sys.stderr)
        return 2

    columns = list(SUMMARY_COLUMNS)
    if args.relative:
        add_relative_energies(rows)
        columns.insert(columns.index("total_energy_eh") + 1, "rel_energy_kcal")

    if args.output:
        write_csv(Path(args.output), columns, rows)
        print(f"Wrote {len(rows)} runs to {args.output}")
    else:
        print("\t".join(columns))
        for r in rows:
            print("\t".join("" if r.get(c) is None else str(r.get(c, "")) for c in columns))
    if args.charges_csv:
        write_csv(Path(args.charges_csv), CHARGE_COLUMNS, atoms)
        print(f"Wrote {len(atoms)} atomic charges to {args.charges_csv}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
