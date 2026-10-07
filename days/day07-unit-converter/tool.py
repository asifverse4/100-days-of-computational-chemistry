"""Day 7: unit converter for energies, wavelengths and wavenumbers.

Convert between hartree, eV, kcal/mol, kJ/mol, cm-1, nm, THz and kelvin using
CODATA 2018 constants. Works on single values or on a CSV column.

Usage:
    python tool.py 1 eh ev kcal/mol kj/mol cm-1   # one value, several targets
    python tool.py 450 nm ev                      # optical gap from an absorption edge
    python tool.py 2.5 ev                         # no target: show every unit
    python tool.py --csv summary.csv --column total_energy_eh --from eh --to kcal/mol
"""
import argparse
import csv
import re
import sys
from pathlib import Path

# CODATA 2018 (exact SI values except the hartree energy)
H = 6.62607015e-34  # Planck constant, J s
C = 299792458.0  # speed of light, m/s
E_CHARGE = 1.602176634e-19  # elementary charge, C (1 eV in J)
N_A = 6.02214076e23  # Avogadro constant, 1/mol
K_B = 1.380649e-23  # Boltzmann constant, J/K
E_H = 4.3597447222071e-18  # hartree energy, J
KCAL = 4184.0  # thermochemical calorie in J

# Joules (per molecule) per one unit; nm is handled separately because it is inverse.
LINEAR = {
    "eh": E_H,
    "ev": E_CHARGE,
    "kcal/mol": KCAL / N_A,
    "kj/mol": 1000.0 / N_A,
    "cm-1": H * C * 100.0,
    "thz": H * 1e12,
    "k": K_B,
}
LABELS = {"eh": "Eh", "ev": "eV", "kcal/mol": "kcal/mol", "kj/mol": "kJ/mol",
          "cm-1": "cm-1", "nm": "nm", "thz": "THz", "k": "K"}
ORDER = ["eh", "ev", "kcal/mol", "kj/mol", "cm-1", "nm", "thz", "k"]

ALIASES = {
    "eh": "eh", "hartree": "eh", "hartrees": "eh", "ha": "eh", "e_h": "eh", "au": "eh",
    "ev": "ev", "electronvolt": "ev", "electronvolts": "ev",
    "kcal/mol": "kcal/mol", "kcal": "kcal/mol", "kcalmol": "kcal/mol",
    "kcalmol-1": "kcal/mol", "kcal/mole": "kcal/mol",
    "kj/mol": "kj/mol", "kj": "kj/mol", "kjmol": "kj/mol", "kjmol-1": "kj/mol",
    "cm-1": "cm-1", "1/cm": "cm-1", "wavenumber": "cm-1", "wavenumbers": "cm-1",
    "nm": "nm", "nanometer": "nm", "nanometers": "nm", "nanometre": "nm",
    "thz": "thz", "terahertz": "thz",
    "k": "k", "kelvin": "k",
}


def canonical(unit: str) -> str:
    """Map a user-typed unit name (eV, kcal/mol, cm^-1, ...) to its internal key."""
    key = unit.strip().lower().replace(" ", "").replace("\u207b\u00b9", "-1")
    key = key.replace("^", "").replace("\u00b7", "").replace(".", "")
    if key not in ALIASES:
        raise ValueError(f"unknown unit {unit!r}. Supported: {', '.join(LABELS[k] for k in ORDER)}")
    return ALIASES[key]


def to_joule(value: float, unit: str) -> float:
    """Energy in joules per molecule."""
    if unit == "nm":
        if value <= 0:
            raise ValueError("a wavelength must be positive")
        return H * C / (value * 1e-9)
    return value * LINEAR[unit]


def from_joule(energy: float, unit: str) -> float:
    if unit == "nm":
        if energy <= 0:
            raise ValueError("only a positive energy has a wavelength")
        return H * C / energy * 1e9
    return energy / LINEAR[unit]


def convert(value: float, src: str, dst: str) -> float:
    """Convert value from unit src to unit dst (names are resolved with canonical())."""
    return from_joule(to_joule(value, canonical(src)), canonical(dst))


def fmt(value: float, digits: int) -> str:
    return format(value, f".{digits}g")


def convert_csv(path: Path, column: str, src: str, targets: list[str]) -> tuple[list[str], list[dict]]:
    """Return (header, rows) with one extra column per target unit."""
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None or column not in reader.fieldnames:
            raise ValueError(f"column {column!r} not found in {path}")
        header = list(reader.fieldnames)
        keys = [canonical(t) for t in targets]
        new_cols = [f"{column}_{re.sub(r'[^a-z0-9]+', '_', LABELS[k].lower()).strip('_')}" for k in keys]
        rows = []
        for n, row in enumerate(reader, 2):
            raw = (row[column] or "").strip()
            for col, key in zip(new_cols, keys):
                if not raw:
                    row[col] = ""
                    continue
                try:
                    row[col] = format(convert(float(raw), src, key), ".10g")
                except ValueError as err:
                    raise ValueError(f"line {n}: {err}") from None
            rows.append(row)
    return header + new_cols, rows


def write_rows(fh, header: list[str], rows: list[dict]) -> None:
    writer = csv.DictWriter(fh, fieldnames=header, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("args", nargs="*", metavar="VALUE FROM [TO ...]",
                    help="a number, its unit, and the unit(s) to convert to (omit to show all units)")
    ap.add_argument("--csv", help="convert a column of this CSV file instead")
    ap.add_argument("--column", help="CSV column to convert")
    ap.add_argument("--from", dest="src", help="unit of the CSV column")
    ap.add_argument("--to", nargs="+", help="target unit(s) for the CSV column")
    ap.add_argument("-o", "--output", help="write the converted CSV here (default: print)")
    ap.add_argument("-d", "--digits", type=int, default=6, help="significant digits to print (default 6)")
    ap.add_argument("--list", action="store_true", help="list supported units and constants")
    args = ap.parse_args()

    try:
        if args.list:
            print("units:", ", ".join(LABELS[k] for k in ORDER))
            print("constants (CODATA 2018): h, c, e, N_A, k_B exact; hartree energy 4.3597447222071e-18 J; 1 cal = 4.184 J")
            return 0

        if args.csv:
            if not (args.column and args.src and args.to):
                raise ValueError("--csv needs --column, --from and --to")
            header, rows = convert_csv(Path(args.csv), args.column, args.src, args.to)
            if args.output:
                with open(args.output, "w", newline="") as fh:
                    write_rows(fh, header, rows)
                print(f"Wrote {len(rows)} rows to {args.output}")
            else:
                write_rows(sys.stdout, header, rows)
            return 0

        if len(args.args) < 2:
            ap.error("give a value and its unit, for example: 450 nm ev")
        try:
            value = float(args.args[0])
        except ValueError:
            raise ValueError(f"not a number: {args.args[0]!r}") from None
        src = canonical(args.args[1])
        targets = [canonical(t) for t in args.args[2:]] or [k for k in ORDER if k != src]
        shown = fmt(value, args.digits)
        for key in targets:
            try:
                result = fmt(convert(value, src, key), args.digits)
            except ValueError:
                result = "n/a (needs a positive value)"
            print(f"{shown} {LABELS[src]} = {result} {LABELS[key]}")
        return 0
    except ValueError as err:
        print(f"error: {err}", file=sys.stderr)
        return 2
    except OSError as err:
        print(f"error: {err}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
