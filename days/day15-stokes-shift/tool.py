"""Day 15: Stokes shift calculator.

The Stokes shift is the gap between the absorption and emission maxima, given in nm, cm-1 and eV.
Enter the two wavelengths, pick the maxima from an absorption and an emission spectrum, or process a
table of samples.

    shift (nm)  = em - abs
    shift (cm-1) = 1e7/abs - 1e7/em
    shift (eV)  = 1239.842 * (1/abs - 1/em)

Usage:
    python tool.py 345 450
    python tool.py --abs example/synthetic_absorption.csv --em example/synthetic_emission.csv --abs-peak longest -o stokes.png
    python tool.py --table example/samples.csv --csv shifts.csv -o shifts.png
    python tool.py 3.59 2.76 --input-unit ev
"""
import argparse
import csv
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks, savgol_filter

HC_EV_NM = 1239.841984  # h*c in eV nm (CODATA 2018)
NM_TO_CM1 = 1.0e7
MIN_POINTS = 5
RESULT_COLUMNS = ["name", "abs_nm", "em_nm", "shift_nm", "shift_cm1", "shift_ev", "shift_cm1_unc", "shift_ev_unc"]
COLORS = {"abs": "#0072B2", "em": "#D55E00"}
EMISSION_RE = re.compile(r"^em|emi|\bpl\b|fluor|photolum")


@dataclass
class Spectrum:
    name: str
    x: np.ndarray
    y: np.ndarray


# ---------- the calculation ----------
def to_nm(value: float, unit: str) -> float:
    """Convert a wavelength (nm), photon energy (eV) or wavenumber (cm-1) to nm."""
    if value <= 0:
        raise ValueError(f"values must be positive, got {value:g}")
    return {"nm": value, "ev": HC_EV_NM / value, "cm-1": NM_TO_CM1 / value}[unit]


def stokes(abs_nm: float, em_nm: float, unc_abs: float = 0.0, unc_em: float = 0.0) -> dict:
    """Stokes shift in nm, cm-1 and eV, with propagated uncertainties for the given wavelength uncertainties."""
    if abs_nm <= 0 or em_nm <= 0:
        raise ValueError("wavelengths must be positive")
    cm1_unc = NM_TO_CM1 * np.sqrt((unc_abs / abs_nm**2) ** 2 + (unc_em / em_nm**2) ** 2)
    return {"abs_nm": abs_nm, "em_nm": em_nm, "shift_nm": em_nm - abs_nm,
            "shift_cm1": NM_TO_CM1 * (1 / abs_nm - 1 / em_nm), "shift_ev": HC_EV_NM * (1 / abs_nm - 1 / em_nm),
            "shift_nm_unc": float(np.hypot(unc_abs, unc_em)), "shift_cm1_unc": float(cm1_unc),
            "shift_ev_unc": float(cm1_unc * HC_EV_NM / NM_TO_CM1)}


def shift_warnings(name: str, r: dict) -> list[str]:
    notes = []
    if r["shift_nm"] < 0:
        notes.append(f"{name}: emission ({r['em_nm']:.1f} nm) is at a shorter wavelength than absorption "
                     f"({r['abs_nm']:.1f} nm): check the bands (anti-Stokes shift or a wrongly picked peak)")
    elif r["shift_nm"] < 5:
        notes.append(f"{name}: the shift is under 5 nm: the bands may overlap with scatter or share one transition")
    return notes


# ---------- reading ----------
def read_spectra(path: Path) -> list[Spectrum]:
    """Read wavelength (nm) in column 1 and one or more intensity columns."""
    try:
        lines = [ln for ln in path.read_text(errors="replace").splitlines() if ln.strip()]
    except OSError as err:
        raise ValueError(f"cannot read {path}: {err.strerror}") from None
    rows, header, ncols = [], None, 0
    for line in lines:
        delim = next((d for d in ("\t", ";", ",") if d in line), None)
        tokens = [t.strip() for t in (line.split(delim) if delim else line.split()) if t.strip()]
        try:
            values = [float(t.replace(",", ".") if delim == ";" else t) for t in tokens] if len(tokens) >= 2 else None
        except ValueError:
            values = None
        if values is not None:
            if not rows:
                ncols = len(values)
            if len(values) == ncols:
                rows.append(values)
        elif not rows:
            header = tokens
    if not rows:
        raise ValueError(f"{path}: no numeric data found")
    data = np.array(rows)
    data = data[np.argsort(data[:, 0], kind="stable")]
    _, keep = np.unique(data[:, 0], return_index=True)
    data = data[np.sort(keep)]
    if data.shape[0] < MIN_POINTS * 2:
        raise ValueError(f"{path}: only {data.shape[0]} data points, need at least {MIN_POINTS * 2}")
    names = header[1:] if header and len(header) == ncols else []
    return [Spectrum(names[i - 1] if i - 1 < len(names) else f"{path.stem}_{i}", data[:, 0].copy(), data[:, i].copy())
            for i in range(1, ncols)]


def pick_column(spectra: list[Spectrum], column: str | None, path: str) -> Spectrum:
    """Choose an intensity column by header name or 1-based number (default: the first)."""
    if column is None:
        return spectra[0]
    for s in spectra:
        if s.name == column:
            return s
    if column.isdigit() and 1 <= int(column) <= len(spectra):
        return spectra[int(column) - 1]
    raise ValueError(f"{path}: column {column!r} not found. Available: {', '.join(s.name for s in spectra)}")


def read_table(path: Path, unit: str) -> list[tuple[str, float, float]]:
    """Read (name, absorption maximum, emission maximum) rows from a CSV, converting the unit to nm."""
    try:
        lines = [ln for ln in path.read_text(errors="replace").splitlines() if ln.strip()]
    except OSError as err:
        raise ValueError(f"cannot read {path}: {err.strerror}") from None
    out, ia, ie, in_ = [], 1, 2, 0
    for k, line in enumerate(lines, 1):
        delim = next((d for d in ("\t", ";", ",") if d in line), None)
        tokens = [t.strip() for t in (line.split(delim) if delim else line.split())]
        try:
            float(tokens[1])
        except (ValueError, IndexError):
            if not out:  # header row: find the columns by name
                low = [t.lower() for t in tokens]
                ia = next((i for i, t in enumerate(low) if "abs" in t), ia)
                ie = next((i for i, t in enumerate(low) if i != ia and EMISSION_RE.search(t)), ie)
                in_ = next((i for i, t in enumerate(low) if i not in (ia, ie) and ("name" in t or "sample" in t)), in_)
                continue
            raise ValueError(f"{path}: line {k}: expected a name and two numbers") from None
        try:
            name, a, e = tokens[in_], float(tokens[ia]), float(tokens[ie])
        except (ValueError, IndexError):
            raise ValueError(f"{path}: line {k}: expected a name, an absorption and an emission value") from None
        out.append((name, to_nm(a, unit), to_nm(e, unit)))
    if not out:
        raise ValueError(f"{path}: no rows found (expected columns name, absorption maximum, emission maximum)")
    return out


# ---------- peak picking ----------
def refine_peak(x: np.ndarray, y: np.ndarray, i: int, top_fraction: float = 0.9) -> float:
    """Peak position from a parabola fitted to the points above 90% of the maximum (same method as Day 14)."""
    threshold = top_fraction * y[i]
    lo = hi = i
    while lo > 0 and y[lo - 1] >= threshold:
        lo -= 1
    while hi < y.size - 1 and y[hi + 1] >= threshold:
        hi += 1
    lo, hi = min(lo, max(i - 2, 0)), max(hi, min(i + 2, y.size - 1))  # always use at least 5 points
    if hi - lo + 1 >= 5:
        a, b, _ = np.polyfit(x[lo:hi + 1] - x[i], y[lo:hi + 1], 2)
        if a < 0 and x[lo] - x[i] <= -b / (2 * a) <= x[hi] - x[i]:
            return float(x[i] - b / (2 * a))
    return float(x[i])


def find_band(spec: Spectrum, mode: str = "max", rng: tuple[float, float] | None = None, smooth: int = 0,
              prominence: float | None = None) -> float:
    """Wavelength of the highest band (mode 'max') or of the longest-wavelength band (mode 'longest')."""
    x, y = spec.x, spec.y
    if rng:
        mask = (x >= rng[0]) & (x <= rng[1])
        if mask.sum() < MIN_POINTS:
            raise ValueError(f"{spec.name}: fewer than {MIN_POINTS} points in {rng[0]:g}-{rng[1]:g} nm")
        x, y = x[mask], y[mask]
    if smooth:
        if smooth % 2 == 0 or smooth < MIN_POINTS or smooth > y.size:
            raise ValueError("--smooth needs an odd window of at least 5 points, shorter than the data")
        y = savgol_filter(y, smooth, 3)
    if mode == "max":
        return refine_peak(x, y, int(np.argmax(y)))
    prom = prominence if prominence is not None else 0.05 * float(y.max() - y.min())
    idx, _ = find_peaks(y, prominence=max(prom, 1e-12))
    if idx.size == 0:
        raise ValueError(f"{spec.name}: no band found; try a smaller --prominence")
    return refine_peak(x, y, int(idx[-1]))


# ---------- output ----------
def plot_spectra(abs_spec: Spectrum, em_spec: Spectrum, r: dict, path: Path, dpi: int) -> None:
    fig, ax = plt.subplots(figsize=(6.0, 4.0))
    for key, spec, label in (("abs", abs_spec, "absorption"), ("em", em_spec, "emission")):
        top = float(np.nanmax(spec.y))
        ax.plot(spec.x, spec.y / top, color=COLORS[key], lw=1.5, label=label)
        ax.axvline(r[f"{key}_nm"], color=COLORS[key], ls=":", lw=1.1)
    ax.annotate("", xy=(r["em_nm"], 1.08), xytext=(r["abs_nm"], 1.08),
                arrowprops={"arrowstyle": "<->", "color": "black", "lw": 1.1})
    ax.text((r["abs_nm"] + r["em_nm"]) / 2, 1.11, f"{r['shift_nm']:.0f} nm = {r['shift_cm1']:.0f} cm$^{{-1}}$",
            ha="center", va="bottom", fontsize=9, bbox={"facecolor": "white", "edgecolor": "none", "pad": 1.5})
    ax.set_ylim(-0.03, 1.28)
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Normalized intensity")
    ax.legend(frameon=False, fontsize=8, loc="upper right", bbox_to_anchor=(1.0, 0.9))
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi)
    plt.close(fig)


def plot_table(rows: list[dict], path: Path, dpi: int) -> None:
    fig, ax = plt.subplots(figsize=(5.5, 0.6 * len(rows) + 1.6))
    names = [r["name"] for r in rows][::-1]
    values = [r["shift_cm1"] for r in rows][::-1]
    ax.barh(names, values, color="#0072B2", height=0.6)
    for k, v in enumerate(values):
        ax.text(v, k, f" {v:.0f}", va="center", fontsize=8)
    ax.set_xlabel("Stokes shift (cm$^{-1}$)")
    ax.set_xlim(0, max(max(values), 1) * 1.15)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi)
    plt.close(fig)


def result_table(rows: list[dict]) -> str:
    width = max([len("sample")] + [len(r["name"]) for r in rows])
    lines = [f"{'sample':<{width}}  {'abs (nm)':>8}  {'em (nm)':>8}  {'shift (nm)':>10}  {'shift (cm-1)':>12}  {'shift (eV)':>10}"]
    for r in rows:
        lines.append(f"{r['name']:<{width}}  {r['abs_nm']:>8.1f}  {r['em_nm']:>8.1f}  {r['shift_nm']:>10.1f}  "
                     f"{r['shift_cm1']:>12.0f}  {r['shift_ev']:>10.3f}")
    return "\n".join(lines)


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=RESULT_COLUMNS, lineterminator="\n")
        writer.writeheader()
        for r in rows:
            has_unc = r["shift_cm1_unc"] > 0
            writer.writerow({"name": r["name"], "abs_nm": round(r["abs_nm"], 3), "em_nm": round(r["em_nm"], 3),
                             "shift_nm": round(r["shift_nm"], 3), "shift_cm1": round(r["shift_cm1"], 1),
                             "shift_ev": round(r["shift_ev"], 4),
                             "shift_cm1_unc": round(r["shift_cm1_unc"], 1) if has_unc else "",
                             "shift_ev_unc": round(r["shift_ev_unc"], 4) if has_unc else ""})


def parse_range(values: list[float] | None, label: str) -> tuple[float, float] | None:
    if values is None:
        return None
    lo, hi = values
    if lo >= hi:
        raise ValueError(f"{label} needs LOW below HIGH")
    return lo, hi


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("values", nargs="*", type=float, metavar="ABS EM", help="absorption and emission maxima (nm by default)")
    ap.add_argument("--input-unit", choices=["nm", "ev", "cm-1"], default="nm", help="unit of the values and of a --table")
    ap.add_argument("--abs", dest="abs_file", help="absorption spectrum file (wavelength in nm, intensity)")
    ap.add_argument("--em", dest="em_file", help="emission spectrum file (wavelength in nm, intensity)")
    ap.add_argument("--abs-column", help="intensity column of --abs (header name or number)")
    ap.add_argument("--em-column", help="intensity column of --em (header name or number)")
    ap.add_argument("--abs-peak", choices=["max", "longest"], default="max",
                    help="highest band (default) or the longest-wavelength band of the absorption spectrum")
    ap.add_argument("--abs-range", type=float, nargs=2, metavar=("LOW", "HIGH"), help="search the absorption band here (nm)")
    ap.add_argument("--em-range", type=float, nargs=2, metavar=("LOW", "HIGH"), help="search the emission band here (nm)")
    ap.add_argument("--smooth", type=int, default=0, metavar="N", help="smooth spectra before picking (odd window)")
    ap.add_argument("--prominence", type=float, help="minimum band prominence for --abs-peak longest (default 5%% of range)")
    ap.add_argument("--table", help="CSV with columns name, absorption maximum, emission maximum")
    ap.add_argument("--unc", type=float, default=0.0, help="uncertainty of each wavelength in nm, for error bars on the shift")
    ap.add_argument("-o", "--output", help="save a plot (.png, .svg, .pdf): spectra with the shift, or bars for a table")
    ap.add_argument("--csv", help="save the results to a CSV")
    ap.add_argument("--dpi", type=int, default=300, help="PNG resolution (default 300)")
    args = ap.parse_args()

    modes = [bool(args.values), bool(args.abs_file or args.em_file), bool(args.table)]
    try:
        if sum(modes) != 1:
            raise ValueError("choose one input: two values (ABS EM), --abs with --em, or --table")
        if args.unc < 0:
            raise ValueError("--unc must not be negative")
        if args.output and Path(args.output).suffix.lower() not in (".png", ".svg", ".pdf"):
            raise ValueError("output must end in .png, .svg or .pdf")
        rows, spectra = [], None
        if args.values:
            if len(args.values) != 2:
                raise ValueError("give exactly two values: the absorption maximum and the emission maximum")
            if args.output:
                raise ValueError("-o needs spectra (--abs and --em) or a --table")
            a, e = (to_nm(v, args.input_unit) for v in args.values)
            rows.append({"name": "sample", **stokes(a, e, args.unc, args.unc)})
        elif args.table:
            for name, a, e in read_table(Path(args.table), args.input_unit):
                rows.append({"name": name, **stokes(a, e, args.unc, args.unc)})
        else:
            if not (args.abs_file and args.em_file):
                raise ValueError("--abs and --em must be given together")
            abs_spec = pick_column(read_spectra(Path(args.abs_file)), args.abs_column, args.abs_file)
            em_spec = pick_column(read_spectra(Path(args.em_file)), args.em_column, args.em_file)
            a = find_band(abs_spec, args.abs_peak, parse_range(args.abs_range, "--abs-range"), args.smooth, args.prominence)
            e = find_band(em_spec, "max", parse_range(args.em_range, "--em-range"), args.smooth)
            rows.append({"name": "sample", **stokes(a, e, args.unc, args.unc)})
            spectra = (abs_spec, em_spec)
        if args.output:
            if spectra:
                plot_spectra(spectra[0], spectra[1], rows[0], Path(args.output), args.dpi)
            else:
                plot_table(rows, Path(args.output), args.dpi)
    except ValueError as err:
        print(f"error: {err}", file=sys.stderr)
        return 2

    print(result_table(rows))
    if args.unc > 0:
        print(f"uncertainty from +/- {args.unc:g} nm on each wavelength: "
              + ", ".join(f"{r['name']} +/- {r['shift_cm1_unc']:.0f} cm-1" for r in rows))
    for r in rows:
        for note in shift_warnings(r["name"], r):
            print(f"[warn] {note}", file=sys.stderr)
    if args.output:
        print(f"[ok] plot -> {args.output}")
    if args.csv:
        write_csv(Path(args.csv), rows)
        print(f"Wrote {len(rows)} results to {args.csv}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except BrokenPipeError:  # e.g. piping into `head`
        sys.exit(0)
