"""Day 13: relative fluorescence quantum yield calculator.

Compares a sample with a reference standard of known quantum yield (IUPAC eq. 1, Brouwer 2011):

    QY_x = QY_ref * (F_x / F_ref) * (f_ref / f_x) * (n_x / n_ref)^2,   f = 1 - 10^(-A)

Gradient method (3 or more concentrations each): F is plotted against f and the slopes are compared.
Single-point method (one measurement each): the ratio is taken directly.

Usage:
    python tool.py --ref example/synthetic_reference.csv --sample example/synthetic_sample.csv \\
        --standard fluorescein --sample-solvent ethanol
    python tool.py --ref-point 0.05 1.2e6 --sample-point 0.05 4.0e5 --standard quinine-hclo4 --sample-solvent water
    python tool.py --list-standards
"""
import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

MIN_GRADIENT_POINTS = 3
PALETTE = {"ref": "#0072B2", "sample": "#D55E00"}
SOURCE = "Brouwer, Pure Appl. Chem. 2011, 83, 2213, Table 1"


@dataclass(frozen=True)
class Standard:
    name: str
    qy: float
    qy_unc: float | None
    solvent: str
    excitation: str  # wavelength used for the reported value
    usable_excitation: str
    emission: str
    reference: str


# Values and conditions from Table 1 of the IUPAC Technical Report (Brouwer 2011). Always cite the primary source.
STANDARDS = {
    "quinine-sulfate": Standard("quinine sulfate in 0.5 M H2SO4", 0.546, None, "water", "366 nm", "280-380 nm",
                                "380-580 nm", f"Melhuish, J. Phys. Chem. 1961, 65, 229 (25 C); {SOURCE}"),
    "quinine-hclo4": Standard("quinine in 0.1 M HClO4", 0.60, 0.02, "water", "347.5 nm", "280-380 nm",
                              "380-580 nm", f"NIST SRM 936, Velapoldi and Mielenz 1980; {SOURCE}"),
    "fluorescein": Standard("fluorescein in 0.1 M NaOH", 0.925, 0.015, "water", "488 or 514 nm", "470-490 nm",
                            "490-620 nm", f"Magde, Wong, Seybold, Photochem. Photobiol. 2002, 75, 327; {SOURCE}"),
    "rhodamine-6g": Standard("rhodamine 6G in ethanol", 0.95, 0.015, "ethanol", "530 nm", "470-510 nm",
                             "510-700 nm", f"Magde, Wong, Seybold, Photochem. Photobiol. 2002, 75, 327; {SOURCE}"),
    "rhodamine-6g-water": Standard("rhodamine 6G in water", 0.92, 0.02, "water", "530 nm", "470-510 nm",
                                   "510-700 nm", f"Magde, Wong, Seybold, Photochem. Photobiol. 2002, 75, 327; {SOURCE}"),
}

# Refractive indices near 20 C at the sodium D line. Dilute aqueous solutions are close to water.
SOLVENTS = {
    "water": 1.333, "methanol": 1.329, "ethanol": 1.361, "acetonitrile": 1.344, "acetone": 1.359,
    "dmso": 1.479, "dmf": 1.430, "thf": 1.407, "chloroform": 1.446, "dichloromethane": 1.424,
    "cyclohexane": 1.426, "hexane": 1.375, "toluene": 1.496,
}


@dataclass
class Measurement:
    a: np.ndarray  # absorbance at the excitation wavelength
    intensity: np.ndarray  # integrated fluorescence intensity


@dataclass
class Gradient:
    slope: float
    sigma: float | None
    intercept: float
    r2: float | None
    n: int


# ---------- reading ----------
def _split(line: str) -> list[str]:
    for d in ("\t", ";", ","):
        if d in line:
            return [t.strip() for t in line.split(d) if t.strip()]
    return line.split()


def read_measurement(path: Path) -> Measurement:
    """Read absorbance and integrated intensity columns (by header name, else the first two columns)."""
    try:
        lines = [ln for ln in path.read_text(errors="replace").splitlines() if ln.strip()]
    except OSError as err:
        raise ValueError(f"cannot read {path}: {err.strerror}") from None
    rows, header = [], None
    for line in lines:
        tokens = _split(line)
        decimal_comma = ";" in line
        try:
            values = [float(t.replace(",", ".") if decimal_comma else t) for t in tokens]
        except ValueError:
            values = None
        if values is not None and len(values) >= 2:
            rows.append(values)
        elif not rows:
            header = [t.lower() for t in tokens]
    if not rows:
        raise ValueError(f"{path}: no numeric data (expected absorbance and integrated intensity columns)")
    ncols = len(rows[0])
    rows = [r for r in rows if len(r) == ncols]
    ia, ii = 0, 1
    if header and len(header) == ncols:
        found_a = [k for k, h in enumerate(header) if "abs" in h]
        found_i = [k for k, h in enumerate(header) if any(w in h for w in ("int", "area", "fluor", "emission", "counts"))]
        if found_a and found_i:
            ia, ii = found_a[0], found_i[0]
    data = np.array(rows)
    return validate(Measurement(data[:, ia], data[:, ii]), str(path))


def validate(m: Measurement, label: str) -> Measurement:
    if np.any(m.a <= 0):
        raise ValueError(f"{label}: absorbance values must be above 0")
    if np.any(m.intensity < 0):
        raise ValueError(f"{label}: integrated intensities must not be negative")
    return m


# ---------- the calculation ----------
def absorption_factor(a: np.ndarray, use_absorbance: bool) -> np.ndarray:
    """f = 1 - 10^-A (IUPAC), or A itself when --absorbance is used."""
    return a.copy() if use_absorbance else 1.0 - 10.0 ** (-a)


def fit_gradient(x: np.ndarray, y: np.ndarray, through_origin: bool = False) -> Gradient:
    """Slope of intensity against absorption factor, with its standard error."""
    n = x.size
    if through_origin:
        slope = float(np.sum(x * y) / np.sum(x * x))
        resid = y - slope * x
        sigma = float(np.sqrt(np.sum(resid**2) / (n - 1) / np.sum(x * x)))
        intercept = 0.0
    else:
        (slope, intercept), cov = np.polyfit(x, y, 1, cov=True)
        sigma, slope, intercept = float(np.sqrt(cov[0, 0])), float(slope), float(intercept)
        resid = y - (slope * x + intercept)
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - float(np.sum(resid**2)) / ss_tot if ss_tot > 0 else None
    if slope <= 0:
        raise ValueError("the fitted slope is not positive; check the intensities")
    return Gradient(slope, sigma, intercept, r2, n)


def resolve_n(solvent: str | None, override: float | None, label: str) -> float:
    if override is not None:
        if override < 1.0:
            raise ValueError(f"{label} refractive index must be at least 1")
        return override
    if solvent is None:
        raise ValueError(f"give the {label} solvent (--{label}-solvent) or its refractive index (--{label}-n)")
    key = solvent.strip().lower()
    if key not in SOLVENTS:
        raise ValueError(f"unknown solvent {solvent!r} for the {label}. Known: {', '.join(sorted(SOLVENTS))}. "
                         f"Or give --{label}-n")
    return SOLVENTS[key]


def compute(ref: Measurement, sample: Measurement, *, ref_qy: float, ref_qy_unc: float | None, n_ref: float,
            n_sample: float, use_absorbance: bool = False, through_origin: bool = False) -> dict:
    """Relative quantum yield by the gradient method (>= 3 points each) or the single-point method (1 each)."""
    if ref_qy <= 0 or ref_qy > 1:
        raise ValueError("the reference quantum yield must be between 0 and 1")
    nr, ns = ref.a.size, sample.a.size
    index = (n_sample / n_ref) ** 2
    fr, fs = absorption_factor(ref.a, use_absorbance), absorption_factor(sample.a, use_absorbance)
    if nr == 1 and ns == 1:
        qy = ref_qy * (sample.intensity[0] / ref.intensity[0]) * (fr[0] / fs[0]) * index
        rel = (ref_qy_unc / ref_qy) if ref_qy_unc else 0.0
        return {"method": "single-point", "qy": float(qy), "uncertainty": float(qy * rel) if ref_qy_unc else None,
                "ref_fit": None, "sample_fit": None, "index_ratio": index, "fr": fr, "fs": fs}
    if nr >= MIN_GRADIENT_POINTS and ns >= MIN_GRADIENT_POINTS:
        gr = fit_gradient(fr, ref.intensity, through_origin)
        gs = fit_gradient(fs, sample.intensity, through_origin)
        qy = ref_qy * (gs.slope / gr.slope) * index
        rel2 = (gs.sigma / gs.slope) ** 2 + (gr.sigma / gr.slope) ** 2
        if ref_qy_unc:
            rel2 += (ref_qy_unc / ref_qy) ** 2
        return {"method": "gradient", "qy": float(qy), "uncertainty": float(qy * np.sqrt(rel2)),
                "ref_fit": gr, "sample_fit": gs, "index_ratio": index, "fr": fr, "fs": fs}
    raise ValueError(f"use either one point each (single-point method) or at least {MIN_GRADIENT_POINTS} points "
                     f"each (gradient method). Got {nr} reference and {ns} sample points")


def warnings_for(ref: Measurement, sample: Measurement, result: dict) -> list[str]:
    notes = []
    if max(ref.a.max(), sample.a.max()) > 0.1:
        notes.append("absorbance above 0.1 at the excitation wavelength: inner-filter effects can bias the "
                     "result (IUPAC recommends A < 0.1)")
    for label, fit, meas in (("reference", result["ref_fit"], ref), ("sample", result["sample_fit"], sample)):
        if fit is None:
            continue
        if fit.r2 is not None and fit.r2 < 0.99:
            notes.append(f"{label} fit is not very linear (R2 = {fit.r2:.3f}): check for inner-filter effects, "
                         "aggregation or outliers")
        if abs(fit.intercept) > 0.05 * meas.intensity.max():
            notes.append(f"{label} fit has a large intercept: check the background subtraction")
        if meas.a.max() / meas.a.min() < 2:
            notes.append(f"{label} absorbances cover a narrow range, so the slope is poorly determined")
    if result["qy"] > 1.0:
        notes.append("quantum yield above 1 is not physical: check the refractive indices, absorbances, "
                     "reference value and that the intensities are integrated areas")
    return notes


# ---------- output ----------
def plot_fit(ref: Measurement, sample: Measurement, result: dict, path: Path, x_label: str, dpi: int) -> None:
    fig, ax = plt.subplots(figsize=(5.0, 3.8))
    for key, label, meas, f, fit in (("ref", "reference", ref, result["fr"], result["ref_fit"]),
                                     ("sample", "sample", sample, result["fs"], result["sample_fit"])):
        ax.plot(f, meas.intensity, "o", color=PALETTE[key], ms=5, label=label)
        slope, intercept = (fit.slope, fit.intercept) if fit else (meas.intensity[0] / f[0], 0.0)
        xs = np.array([0.0, f.max() * 1.05])
        ax.plot(xs, slope * xs + intercept, "-" if fit else "--", color=PALETTE[key], lw=1.2,
                label=f"{label} slope {slope:.3g}")
    unc = f" +/- {result['uncertainty']:.3f}" if result["uncertainty"] else ""
    ax.set_title(f"Relative QY = {result['qy']:.3f}{unc} ({result['method']})", fontsize=10)
    ax.set_xlabel(x_label)
    ax.set_ylabel("Integrated fluorescence intensity (a.u.)")
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False, fontsize=8)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi)
    plt.close(fig)


def standards_table() -> str:
    intro = ("Reference standards (Brouwer, Pure Appl. Chem. 2011, 83, 2213, Table 1). "
             "Cite the primary source and check the conditions before using one.")
    lines = [intro, ""]
    for key, s in STANDARDS.items():
        unc = f" +/- {s.qy_unc}" if s.qy_unc else ""
        detail = f"    QY {s.qy}{unc} at {s.excitation} excitation; usable excitation {s.usable_excitation}; emission {s.emission}"
        lines += [f"{key}: {s.name}", detail, f"    {s.reference}"]
    return "\n".join(lines)


def write_csv(path: Path, row: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(row), lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)


def result_row(result: dict, ref_qy: float, n_ref: float, n_sample: float, standard: str) -> dict:
    gr, gs = result["ref_fit"], result["sample_fit"]
    return {"method": result["method"], "quantum_yield": round(result["qy"], 4),
            "uncertainty": "" if result["uncertainty"] is None else round(result["uncertainty"], 4),
            "reference": standard, "reference_qy": ref_qy, "n_reference": n_ref, "n_sample": n_sample,
            "ref_slope": "" if gr is None else float(f"{gr.slope:.6g}"),
            "sample_slope": "" if gs is None else float(f"{gs.slope:.6g}"),
            "ref_r2": "" if gr is None or gr.r2 is None else round(gr.r2, 5),
            "sample_r2": "" if gs is None or gs.r2 is None else round(gs.r2, 5)}


def load(file_arg: str | None, point: list[float] | None, label: str) -> Measurement:
    if file_arg and point:
        raise ValueError(f"give the {label} as a file or as a point, not both")
    if file_arg:
        return read_measurement(Path(file_arg))
    if point:
        return validate(Measurement(np.array([point[0]]), np.array([point[1]])), f"{label} point")
    raise ValueError(f"give the {label} data with --{label} FILE or --{label}-point ABSORBANCE INTENSITY")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", help="reference data: absorbance and integrated intensity columns")
    ap.add_argument("--sample", help="sample data: absorbance and integrated intensity columns")
    ap.add_argument("--ref-point", type=float, nargs=2, metavar=("A", "I"), help="single reference measurement")
    ap.add_argument("--sample-point", type=float, nargs=2, metavar=("A", "I"), help="single sample measurement")
    ap.add_argument("--standard", choices=sorted(STANDARDS), help="built-in reference standard (see --list-standards)")
    ap.add_argument("--ref-qy", type=float, help="reference quantum yield (overrides --standard)")
    ap.add_argument("--ref-qy-unc", type=float, help="uncertainty of the reference quantum yield")
    ap.add_argument("--ref-solvent", help="reference solvent (default: the standard's solvent)")
    ap.add_argument("--sample-solvent", help="sample solvent")
    ap.add_argument("--ref-n", type=float, help="reference refractive index (overrides --ref-solvent)")
    ap.add_argument("--sample-n", type=float, help="sample refractive index (overrides --sample-solvent)")
    ap.add_argument("--absorbance", action="store_true", help="use absorbance A instead of the absorption factor 1-10^-A")
    ap.add_argument("--through-origin", action="store_true", help="force the gradient fits through zero")
    ap.add_argument("--plot", help="save the intensity-vs-absorption-factor plot (.png, .svg or .pdf)")
    ap.add_argument("--csv", help="save the result to a one-row CSV")
    ap.add_argument("--dpi", type=int, default=300, help="PNG resolution (default 300)")
    ap.add_argument("--list-standards", action="store_true", help="show the built-in standards and exit")
    ap.add_argument("--list-solvents", action="store_true", help="show the built-in refractive indices and exit")
    args = ap.parse_args()

    if args.list_standards:
        print(standards_table())
        return 0
    if args.list_solvents:
        print("\n".join(f"{k:<16}{v}" for k, v in sorted(SOLVENTS.items())))
        return 0
    try:
        std = STANDARDS.get(args.standard) if args.standard else None
        ref_qy = args.ref_qy if args.ref_qy is not None else (std.qy if std else None)
        if ref_qy is None:
            raise ValueError("give the reference quantum yield with --standard or --ref-qy")
        ref_unc = args.ref_qy_unc if args.ref_qy_unc is not None else (std.qy_unc if std and args.ref_qy is None else None)
        n_ref = resolve_n(args.ref_solvent or (std.solvent if std else None), args.ref_n, "ref")
        n_sample = resolve_n(args.sample_solvent, args.sample_n, "sample")
        ref, sample = load(args.ref, args.ref_point, "ref"), load(args.sample, args.sample_point, "sample")
        result = compute(ref, sample, ref_qy=ref_qy, ref_qy_unc=ref_unc, n_ref=n_ref, n_sample=n_sample,
                         use_absorbance=args.absorbance, through_origin=args.through_origin)
        if args.plot:
            if Path(args.plot).suffix.lower() not in (".png", ".svg", ".pdf"):
                raise ValueError("--plot must end in .png, .svg or .pdf")
            plot_fit(ref, sample, result, Path(args.plot),
                     "Absorbance A" if args.absorbance else "Absorption factor 1 - 10^-A", args.dpi)
    except ValueError as err:
        print(f"error: {err}", file=sys.stderr)
        return 2

    label = std.name if std and args.ref_qy is None else "custom reference"
    print(f"Relative quantum yield ({result['method']} method)")
    print(f"  reference:        {label}, QY {ref_qy}")
    print(f"  refractive index: reference {n_ref}, sample {n_sample}, (n_s/n_ref)^2 = {result['index_ratio']:.4f}")
    for name, fit in (("reference", result["ref_fit"]), ("sample", result["sample_fit"])):
        if fit:
            print(f"  gradient {name + ':':<10} {fit.slope:.4g} +/- {fit.sigma:.2g} (R2 {fit.r2:.4f}, {fit.n} points)")
    unc = f" +/- {result['uncertainty']:.3f}" if result["uncertainty"] is not None else ""
    print(f"  QY sample = {result['qy']:.3f}{unc}")
    for note in warnings_for(ref, sample, result):
        print(f"[warn] {note}", file=sys.stderr)
    if args.plot:
        print(f"[ok] plot -> {args.plot}")
    if args.csv:
        write_csv(Path(args.csv), result_row(result, ref_qy, n_ref, n_sample, label))
        print(f"Wrote result to {args.csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
