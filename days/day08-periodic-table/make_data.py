"""Regenerate elements.csv from the `mendeleev` package (developer tool, not needed to use Day 8).

    pip install mendeleev
    python make_data.py

Rounds values to sensible precision and leaves missing values blank.
"""
import csv
import warnings
from pathlib import Path

from mendeleev import get_all_elements

warnings.filterwarnings("ignore")  # mendeleev warns about elements with several allotropes

COLUMNS = [
    "atomic_number", "symbol", "name", "atomic_mass", "group", "period", "block", "category",
    "electronegativity", "covalent_radius_pm", "vdw_radius_pm", "electron_configuration",
    "melting_point_k", "boiling_point_k", "density_g_cm3",
]


def num(value, digits=None, sig=None):
    if value is None:
        return ""
    if sig is not None:
        return f"{float(value):.{sig}g}"
    return round(float(value), digits) if digits else round(float(value))


def main() -> None:
    rows = []
    for e in get_all_elements():
        rows.append({
            "atomic_number": e.atomic_number,
            "symbol": e.symbol,
            "name": e.name,
            "atomic_mass": num(e.atomic_weight, 6),
            "group": "" if e.group_id is None else int(e.group_id),
            "period": e.period,
            "block": e.block,
            "category": e.series,
            "electronegativity": num(e.electronegativity_pauling(), 2),
            "covalent_radius_pm": num(e.covalent_radius_pyykko),
            "vdw_radius_pm": num(e.vdw_radius_alvarez),
            "electron_configuration": e.econf,
            "melting_point_k": num(e.melting_point, 2),
            "boiling_point_k": num(e.boiling_point, 2),
            "density_g_cm3": num(e.density, sig=6),
        })
    out = Path(__file__).with_name("elements.csv")
    with out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} elements to {out}")


if __name__ == "__main__":
    main()
