"""Day 8: periodic-table property lookup.

Look up an element by symbol, name or atomic number, list elements by group,
period, block or category, or get the molar mass of a formula. Offline: the data
ships as elements.csv (see the README for sources and the license).

Usage:
    python tool.py Fe
    python tool.py 26 cu zinc --property electronegativity
    python tool.py --group 17
    python tool.py --block d --period 4 -o 3d_metals.csv
    python tool.py --formula "CuSO4.5H2O"
"""
import argparse
import csv
import difflib
import re
import sys
from functools import lru_cache
from pathlib import Path

DATA = Path(__file__).with_name("elements.csv")
INT_COLS = {"atomic_number", "group", "period"}
TEXT_COLS = {"symbol", "name", "block", "category", "electron_configuration"}

# British / alternative spellings people actually type
ALT_NAMES = {"aluminium": "Al", "caesium": "Cs", "sulphur": "S", "wolfram": "W"}
# short names for --property
PROPERTY_ALIASES = {
    "z": "atomic_number", "mass": "atomic_mass", "weight": "atomic_mass", "en": "electronegativity",
    "covalent": "covalent_radius_pm", "vdw": "vdw_radius_pm", "config": "electron_configuration",
    "mp": "melting_point_k", "bp": "boiling_point_k", "density": "density_g_cm3",
}
# label, column, unit for the single-element card
CARD = [
    ("atomic mass", "atomic_mass", "u"),
    ("electronegativity", "electronegativity", "(Pauling)"),
    ("electron config.", "electron_configuration", ""),
    ("covalent radius", "covalent_radius_pm", "pm"),
    ("van der Waals radius", "vdw_radius_pm", "pm"),
    ("melting point", "melting_point_k", "K"),
    ("boiling point", "boiling_point_k", "K"),
    ("density", "density_g_cm3", "g/cm3"),
]
TABLE_COLS = [("Z", "atomic_number"), ("Sym", "symbol"), ("Name", "name"), ("Mass", "atomic_mass"),
              ("Grp", "group"), ("Per", "period"), ("Blk", "block"), ("EN", "electronegativity"),
              ("Category", "category")]


@lru_cache(maxsize=1)
def load_elements() -> list[dict]:
    """Read elements.csv into a list of dicts (numbers converted, blanks become None)."""
    with DATA.open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        for key, val in r.items():
            if val == "":
                r[key] = None
            elif key in INT_COLS:
                r[key] = int(val)
            elif key not in TEXT_COLS:
                r[key] = float(val)
    return rows


def find_element(query: str) -> dict:
    """Find an element by atomic number, symbol or name (case-insensitive)."""
    elements, q = load_elements(), query.strip().lower()
    q = ALT_NAMES.get(q, q).lower()
    for e in elements:
        if q == str(e["atomic_number"]) or q == e["symbol"].lower() or q == e["name"].lower():
            return e
    names = [e["name"] for e in elements] + [e["symbol"] for e in elements]
    close = difflib.get_close_matches(query.strip().capitalize(), names, n=2, cutoff=0.7)
    hint = f" Did you mean {' or '.join(close)}?" if close else ""
    raise ValueError(f"unknown element {query!r}.{hint}")


def select(elements: list[dict], group=None, period=None, block=None, category=None) -> list[dict]:
    out = []
    for e in elements:
        if group is not None and e["group"] != group:
            continue
        if period is not None and e["period"] != period:
            continue
        if block is not None and e["block"] != block:
            continue
        if category is not None and category.lower() not in e["category"].lower():
            continue
        out.append(e)
    return out


def parse_formula(formula: str) -> dict[str, int]:
    """Count atoms in a formula such as C9H8O4, Ca(OH)2, K4[Fe(CN)6] or CuSO4.5H2O."""
    total: dict[str, int] = {}
    for part in formula.replace(" ", "").replace("\u00b7", ".").split("."):
        lead = re.match(r"(\d+)(.+)", part)
        mult, body = (int(lead[1]), lead[2]) if lead else (1, part)
        for symbol, n in _count_group(body).items():
            total[symbol] = total.get(symbol, 0) + n * mult
    if not total:
        raise ValueError("empty formula")
    return total


def _count_group(text: str) -> dict[str, int]:
    stack: list[dict[str, int]] = [{}]
    i = 0
    while i < len(text):
        ch = text[i]
        if ch in "([{":
            stack.append({})
            i += 1
        elif ch in ")]}":
            if len(stack) == 1:
                raise ValueError(f"unbalanced bracket at position {i + 1}")
            group = stack.pop()
            num = re.match(r"\d+", text[i + 1:])
            k = int(num[0]) if num else 1
            i += 1 + (len(num[0]) if num else 0)
            for symbol, n in group.items():
                stack[-1][symbol] = stack[-1].get(symbol, 0) + n * k
        else:
            sym = re.match(r"[A-Z][a-z]?", text[i:])
            if not sym:
                raise ValueError(f"unexpected {ch!r} at position {i + 1}")
            i += len(sym[0])
            num = re.match(r"\d+", text[i:])
            n = int(num[0]) if num else 1
            i += len(num[0]) if num else 0
            stack[-1][sym[0]] = stack[-1].get(sym[0], 0) + n
    if len(stack) != 1:
        raise ValueError("unbalanced bracket")
    return stack[0]


def molar_mass(formula: str) -> tuple[float, list[tuple[str, int, float]]]:
    """Return (molar mass in g/mol, [(symbol, count, mass contribution), ...])."""
    parts = []
    for symbol, n in parse_formula(formula).items():
        e = find_element(symbol)
        parts.append((e["symbol"], n, n * e["atomic_mass"]))
    return sum(p[2] for p in parts), parts


def show(value) -> str:
    return "n/a" if value is None else (f"{value:.10g}" if isinstance(value, float) else str(value))


def card(e: dict) -> str:
    where = f"{e['block']}-block, {e['category']}"
    group = show(e["group"]) if e["group"] else "none (f-block)"
    lines = [f"{e['symbol']}  {e['name']}  (Z = {e['atomic_number']})",
             f"  {'group / period':<22}{group} / {e['period']}   ({where})"]
    for label, col, unit in CARD:
        lines.append(f"  {label:<22}{show(e[col])} {unit}".rstrip() if e[col] is not None
                     else f"  {label:<22}n/a")
    return "\n".join(lines)


def table(elements: list[dict]) -> str:
    cells = [[h for h, _ in TABLE_COLS]] + [[show(e[c]) for _, c in TABLE_COLS] for e in elements]
    widths = [max(len(row[i]) for row in cells) for i in range(len(TABLE_COLS))]
    return "\n".join("  ".join(v.ljust(w) for v, w in zip(row, widths)).rstrip() for row in cells)


def write_csv(path: Path, elements: list[dict]) -> None:
    columns = list(load_elements()[0])
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        for e in elements:
            writer.writerow({k: ("" if v is None else v) for k, v in e.items()})


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("queries", nargs="*", help="symbol, name or atomic number (several allowed)")
    ap.add_argument("--formula", help="molar mass of a formula, e.g. C9H8O4 or CuSO4.5H2O")
    ap.add_argument("--group", type=int, help="list elements in this group (1-18)")
    ap.add_argument("--period", type=int, help="list elements in this period (1-7)")
    ap.add_argument("--block", choices=["s", "p", "d", "f"], help="list elements in this block")
    ap.add_argument("--category", help="list elements whose category contains this text, e.g. halogens")
    ap.add_argument("--property", help="print just one property (mass, en, mp, bp, density, covalent, vdw, ...)")
    ap.add_argument("-o", "--output", help="also save the selected elements, with every column, as CSV")
    args = ap.parse_args()

    try:
        if args.formula:
            total, parts = molar_mass(args.formula)
            print(f"{args.formula}: {total:.3f} g/mol")
            for symbol, n, mass in parts:
                print(f"  {symbol:<3}x {n:<4}{mass:9.3f} g/mol  {100 * mass / total:6.2f} %")
            return 0

        filters = (args.group, args.period, args.block, args.category)
        if not args.queries and all(f is None for f in filters):
            ap.error("give an element, a filter (--group, --period, --block, --category) or --formula")
        chosen = [find_element(q) for q in args.queries] if args.queries else load_elements()
        chosen = select(chosen, *filters)
        if not chosen:
            print("no elements match", file=sys.stderr)
            return 1

        if args.property:
            col = PROPERTY_ALIASES.get(args.property.lower(), args.property.lower())
            if col not in chosen[0]:
                raise ValueError(f"unknown property {args.property!r}. Choose from: {', '.join(chosen[0])}")
            for e in chosen:
                print(show(e[col]) if len(chosen) == 1 else f"{e['symbol']}\t{show(e[col])}")
        elif len(chosen) == 1 and args.queries:
            print(card(chosen[0]))
        else:
            print(table(chosen))
        if args.output:
            write_csv(Path(args.output), chosen)
            print(f"Wrote {len(chosen)} elements to {args.output}")
        return 0
    except ValueError as err:
        print(f"error: {err}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    try:
        sys.exit(main())
    except BrokenPipeError:  # e.g. piping into `head`
        sys.exit(0)
