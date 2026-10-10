"""Day 12: Tauc plot and optical band gap calculator.

Turn a UV-vis absorbance (or diffuse-reflectance) spectrum into a Tauc plot, fit the linear
region and extrapolate it to the energy axis to get the optical band gap Eg.

    (alpha * h * nu)^r = B * (h * nu - Eg)
    r = 2 direct allowed, 1/2 indirect allowed, 2/3 direct forbidden, 1/3 indirect forbidden

Usage:
    python tool.py example/synthetic_tauc.csv -o tauc.png
    python tool.py spectrum.csv --transition indirect --zero-at 700 --csv gaps.csv
    python tool.py spectrum.csv --fit-range 3.3 3.7          # choose the linear region yourself
    python tool.py drs.csv --signal reflectance-percent       # Kubelka-Munk for powders
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

HC_EV_NM = 1239.841984  # h*c in eV nm (CODATA 2018)
LN10 = 2.302585093
EXPONENTS = {"direct": 2.0, "indirect": 0.5, "direct-forbidden": 2.0 / 3.0, "indirect-forbidden": 1.0 / 3.0}
RESULT_COLUMNS = ["sample", "transition", "exponent", "band_gap_ev", "uncertainty_ev", "fit_emin_ev",
                  "fit_emax_ev", "r_squared", "n_points", "method", "slope", "intercept", "warning"]
MIN_POINTS = 5
PALETTE = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"]


@dataclass
class Spectrum:
    name: str
    x: np.ndarray
    y: np.ndarray


@dataclass
class Fit:
    slope: float
    intercept: float
    r2: float
    n: int
    lo: float
    hi: float
    gap: float
    sigma: float


# ---------- reading (same reader as Day 11) ----------
def _delimiter(line: str) -> str | None:
    for d in ("\t", ";", ","):
        if d in line:
            return d
    return None


def _numbers(tokens: list[str], decimal_comma: bool) -> list[float] | None:
    try:
        return [float(t.replace(",", ".") if decimal_comma else t) for t in tokens]
    except ValueError:
        return None


def read_spectra(path: Path) -> list[Spectrum]:
    """Read a text file with x in column 1 and one or more signal columns."""
    try:
        lines = [ln for ln in path.read_text(errors="replace").splitlines() if ln.strip()]
    except OSError as err:
        raise ValueError(f"cannot read {path}: {err.strerror}") from None
    rows, header, ncols = [], None, 0
    for line in lines:
        delim = _delimiter(line)
        tokens = [t.strip() for t in (line.split(delim) if delim else line.split())]
        tokens = [t for t in tokens if t != ""]
        values = _numbers(tokens, decimal_comma=(delim == ";")) if len(tokens) >= 2 else None
        if values is not None:
            if not rows:
                ncols = len(values)
            if len(values) == ncols:
                rows.append(values)
        elif not rows:
            header = tokens
    if not rows:
        raise ValueError(f"{path}: no numeric data found (expected x and signal columns)")
    data = np.array(rows)
    data = data[np.argsort(data[:, 0], kind="stable")]
    _, keep = np.unique(data[:, 0], return_index=True)
    data = data[np.sort(keep)]
    if data.shape[0] < MIN_POINTS * 2:
        raise ValueError(f"{path}: only {data.shape[0]} data points, need at least {MIN_POINTS * 2}")
    names = header[1:] if header and len(header) == ncols else []
    return [Spectrum(names[i - 1] if i - 1 < len(names) else (path.stem if ncols == 2 else f"{path.stem}_{i}"),
                     data[:, 0].copy(), data[:, i].copy()) for i in range(1, ncols)]


# ---------- physics ----------
def to_energy(spec: Spectrum, x_unit: str) -> Spectrum:
    """Return the spectrum on an ascending photon-energy axis in eV."""
    if np.any(spec.x <= 0):
        raise ValueError(f"{spec.name}: wavelengths and energies must be positive")
    x = HC_EV_NM / spec.x if x_unit == "nm" else spec.x.copy()
    order = np.argsort(x)
    return Spectrum(spec.name, x[order], spec.y[order])


def absorption_proxy(y: np.ndarray, signal: str, path_length_cm: float) -> np.ndarray:
    """alpha (cm^-1) from absorbance, or the Kubelka-Munk function F(R) from reflectance."""
    if signal == "absorbance":
        if path_length_cm <= 0:
            raise ValueError("--path-length-cm must be positive")
        return LN10 * y / path_length_cm
    r = y / 100.0 if signal == "reflectance-percent" else y
    if np.any(r <= 0) or np.any(r >= 1):
        raise ValueError("reflectance must be between 0 and 1 (or 0 and 100 for reflectance-percent), exclusive")
    return (1.0 - r) ** 2 / (2.0 * r)


def tauc_y(energy: np.ndarray, alpha: np.ndarray, transition: str) -> np.ndarray:
    """(alpha*h*nu)^r. Negative values (noise below zero) are clipped to zero."""
    return np.clip(alpha * energy, 0.0, None) ** EXPONENTS[transition]


def linear_fit(e: np.ndarray, y: np.ndarray) -> Fit:
    """Least-squares line through (e, y); the band gap is the x-intercept with propagated uncertainty."""
    (m, b), cov = np.polyfit(e, y, 1, cov=True)
    if m <= 0:
        raise ValueError("the fitted slope is not positive")
    resid = y - (m * e + b)
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - float(np.sum(resid**2)) / ss_tot if ss_tot > 0 else 0.0
    var = cov[1, 1] / m**2 + (b**2 / m**4) * cov[0, 0] - 2 * (b / m**3) * cov[0, 1]
    return Fit(float(m), float(b), r2, int(e.size), float(e.min()), float(e.max()), float(-b / m),
               float(np.sqrt(max(var, 0.0))))


def manual_fit(e: np.ndarray, y: np.ndarray, lo: float, hi: float) -> Fit:
    mask = (e >= lo) & (e <= hi)
    if mask.sum() < MIN_POINTS:
        raise ValueError(f"fewer than {MIN_POINTS} points between {lo} and {hi} eV")
    return linear_fit(e[mask], y[mask])


def auto_fit(e: np.ndarray, y: np.ndarray, width_ev: float = 0.25, r2_min: float = 0.99,
             emin: float | None = None, emax: float | None = None) -> tuple[Fit, bool]:
    """Slide a full-width window along the curve and pick the linear region.

    Among windows with R^2 >= r2_min, keep the steep ones (slope >= 80% of the steepest) and take the one
    at the lowest energy, i.e. the linear region closest to the absorption edge. If no window reaches
    r2_min, the best-R^2 window is used.
    Returns (fit, met_r2_threshold)."""
    mask = np.ones_like(e, dtype=bool)
    if emin is not None:
        mask &= e >= emin
    if emax is not None:
        mask &= e <= emax
    e, y = e[mask], y[mask]
    if e.size < MIN_POINTS:
        raise ValueError(f"fewer than {MIN_POINTS} points in the search range")
    if e[-1] - e[0] < width_ev:
        raise ValueError(f"the energy range ({e[-1] - e[0]:.2f} eV) is narrower than --window {width_ev}; use a smaller --window or --fit-range")
    fits = []
    for i in range(e.size):
        if e[i] + width_ev > e[-1] + 1e-9:
            break  # only full-width windows: short windows at the edge can look steep by chance
        j = int(np.searchsorted(e, e[i] + width_ev, side="right"))
        if j - i < MIN_POINTS:
            continue
        try:
            fits.append(linear_fit(e[i:j], y[i:j]))
        except ValueError:
            continue
    if not fits:
        raise ValueError("no linear region with a positive slope found; try --fit-range or a different --transition")
    good = [f for f in fits if f.r2 >= r2_min]
    if not good:
        return max(fits, key=lambda f: f.r2), False
    steepest = max(f.slope for f in good)
    steep = [f for f in good if f.slope >= 0.8 * steepest]
    return min(steep, key=lambda f: f.lo), True


def warning_for(fit: Fit, met_threshold: bool, r2_min: float) -> str:
    notes = []
    if not met_threshold:
        notes.append(f"no window reached R2 {r2_min}; best window used")
    if fit.gap <= 0 or fit.gap > fit.hi:
        notes.append("intercept is outside the sensible range, check the fit range and transition type")
    elif fit.lo - fit.gap > 1.0:
        notes.append("intercept is more than 1 eV below the fit range")
    return "; ".join(notes)


def analyze(spec: Spectrum, *, transition: str, signal: str, path_length_cm: float, x_unit: str,
            zero_at_nm: float | None, fit_range: tuple[float, float] | None, emin, emax, width_ev, r2_min):
    """Full pipeline for one spectrum. Returns (result row, energy, tauc y, fit)."""
    values = spec.y.copy()
    if zero_at_nm is not None:
        if signal != "absorbance":
            raise ValueError("--zero-at only applies to absorbance spectra")
        if x_unit != "nm":
            raise ValueError("--zero-at needs wavelengths in nm")
        values = values - values[int(np.argmin(np.abs(spec.x - zero_at_nm)))]
    prepared = to_energy(Spectrum(spec.name, spec.x, values), x_unit)
    alpha = absorption_proxy(prepared.y, signal, path_length_cm)
    y = tauc_y(prepared.x, alpha, transition)
    try:
        if fit_range:
            fit, met, method = manual_fit(prepared.x, y, *fit_range), True, "manual"
        else:
            (fit, met), method = auto_fit(prepared.x, y, width_ev, r2_min, emin, emax), "auto"
    except ValueError as err:
        raise ValueError(f"{spec.name}: {err}") from None
    row = {"sample": spec.name, "transition": transition, "exponent": round(EXPONENTS[transition], 4),
           "band_gap_ev": round(fit.gap, 4), "uncertainty_ev": round(fit.sigma, 4),
           "fit_emin_ev": round(fit.lo, 4), "fit_emax_ev": round(fit.hi, 4), "r_squared": round(fit.r2, 5),
           "n_points": fit.n, "method": method, "slope": float(f"{fit.slope:.6g}"),
           "intercept": float(f"{fit.intercept:.6g}"),
           "warning": warning_for(fit, met, r2_min) if method == "auto" else warning_for(fit, True, r2_min)}
    return row, prepared.x, y, fit


# ---------- output ----------
def plot_tauc(items, path: Path, *, title: str | None, dpi: int, figsize: tuple[float, float] | None) -> None:
    """items: list of (row, energy, y, fit). One panel per spectrum."""
    n = len(items)
    cols = 1 if n == 1 else 2
    rows = (n + cols - 1) // cols
    size = figsize or (5.0 * cols, 3.8 * rows)
    fig, axes = plt.subplots(rows, cols, figsize=size, squeeze=False)
    for ax in axes.flat[n:]:
        ax.set_visible(False)
    for k, (row, e, y, fit) in enumerate(items):
        ax = axes.flat[k]
        color = PALETTE[k % len(PALETTE)]
        x_lo, x_hi = max(e.min(), fit.gap - 0.7), min(e.max(), fit.hi + 0.7)
        view = (e >= x_lo) & (e <= x_hi)
        top = 1.2 * float(y[view].max()) if view.any() else 1.0
        ax.plot(e, y, color=color, lw=1.4, label="data")
        ax.axvspan(fit.lo, fit.hi, color=color, alpha=0.12, lw=0)
        line_x = np.array([max(fit.gap, x_lo), fit.hi + 0.15])
        ax.plot(line_x, fit.slope * line_x + fit.intercept, "--", color="black", lw=1.2, label="linear fit")
        ax.plot([fit.gap], [0], "o", color="black", ms=5)
        ax.annotate(f"Eg = {fit.gap:.2f} eV", (fit.gap, 0), textcoords="offset points", xytext=(-95, 45), fontsize=9,
                    arrowprops={"arrowstyle": "-", "lw": 0.8, "color": "black"})
        ax.set_xlim(x_lo, x_hi)
        ax.set_ylim(-0.04 * top, top)
        ax.set_xlabel("Photon energy (eV)")
        r = EXPONENTS[row["transition"]]
        ax.set_ylabel(rf"$(\alpha h\nu)^{{{r:.3g}}}$ (arb. units)")
        ax.set_title(f"{row['sample']} ({row['transition']})", fontsize=10)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        if k == 0:
            ax.legend(frameon=False, fontsize=8, loc="upper left")
    if title:
        fig.suptitle(title)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi)
    plt.close(fig)


def result_table(rows: list[dict]) -> str:
    width = max([len("sample")] + [len(r["sample"]) for r in rows])
    head = f"{'sample':<{width}}  {'transition':<10}  {'Eg (eV)':>8}  {'+/- (eV)':>8}  {'fit range (eV)':>15}  {'R2':>7}  {'n':>4}"
    lines = [head]
    for r in rows:
        rng = f"{r['fit_emin_ev']:.2f}-{r['fit_emax_ev']:.2f}"
        lines.append(f"{r['sample']:<{width}}  {r['transition']:<10}  {r['band_gap_ev']:>8.3f}  "
                     f"{r['uncertainty_ev']:>8.4f}  {rng:>15}  {r['r_squared']:>7.4f}  {r['n_points']:>4}")
    return "\n".join(lines)


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=RESULT_COLUMNS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", help="spectrum file: wavelength (nm) or energy (eV) in column 1, signal in the others")
    ap.add_argument("-o", "--output", default="tauc.png", help="plot file: .png, .svg or .pdf (default tauc.png)")
    ap.add_argument("--csv", help="also write the results to this CSV")
    ap.add_argument("--columns", nargs="+", help="signal columns to analyze (default: all)")
    ap.add_argument("--transition", choices=sorted(EXPONENTS), default="direct", help="transition type (default: direct)")
    ap.add_argument("--signal", choices=["absorbance", "reflectance", "reflectance-percent"], default="absorbance",
                    help="what the columns contain (reflectance uses Kubelka-Munk)")
    ap.add_argument("--path-length-cm", type=float, default=1.0, help="cuvette path length for absorbance (default 1 cm)")
    ap.add_argument("--x-unit", choices=["nm", "ev"], default="nm", help="unit of column 1 (default nm)")
    ap.add_argument("--zero-at", type=float, metavar="NM", help="subtract the absorbance at this wavelength first (e.g. 700)")
    ap.add_argument("--fit-range", type=float, nargs=2, metavar=("EMIN", "EMAX"), help="fit this energy range (eV) instead of searching")
    ap.add_argument("--emin", type=float, help="lowest energy (eV) the automatic search may use")
    ap.add_argument("--emax", type=float, help="highest energy (eV) the automatic search may use")
    ap.add_argument("--window", type=float, default=0.25, help="width in eV of the sliding fit window (default 0.25)")
    ap.add_argument("--r2-min", type=float, default=0.99, help="minimum R2 for the automatic fit (default 0.99)")
    ap.add_argument("--title", help="figure title")
    ap.add_argument("--figsize", type=float, nargs=2, metavar=("W", "H"), help="figure size in inches")
    ap.add_argument("--dpi", type=int, default=300, help="PNG resolution (default 300)")
    args = ap.parse_args()

    out = Path(args.output)
    if out.suffix.lower() not in (".png", ".svg", ".pdf"):
        print("error: output must end in .png, .svg or .pdf", file=sys.stderr)
        return 2
    try:
        spectra = read_spectra(Path(args.input))
        if args.columns:
            missing = [c for c in args.columns if c not in {s.name for s in spectra}]
            if missing:
                raise ValueError(f"column(s) not found: {', '.join(missing)}. Available: {', '.join(s.name for s in spectra)}")
            spectra = [s for s in spectra if s.name in args.columns]
        items = [analyze(s, transition=args.transition, signal=args.signal, path_length_cm=args.path_length_cm,
                         x_unit=args.x_unit, zero_at_nm=args.zero_at, fit_range=args.fit_range, emin=args.emin,
                         emax=args.emax, width_ev=args.window, r2_min=args.r2_min) for s in spectra]
        plot_tauc(items, out, title=args.title, dpi=args.dpi, figsize=tuple(args.figsize) if args.figsize else None)
    except ValueError as err:
        print(f"error: {err}", file=sys.stderr)
        return 2

    rows = [item[0] for item in items]
    print(f"[ok] {len(rows)} spectra -> {out}")
    print(result_table(rows))
    for r in rows:
        if r["warning"]:
            print(f"[warn] {r['sample']}: {r['warning']}", file=sys.stderr)
    if args.csv:
        write_csv(Path(args.csv), rows)
        print(f"Wrote {len(rows)} results to {args.csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
