"""Day 17: XRD Scherrer crystallite size calculator.

Fit the diffraction peaks of a powder pattern, subtract the instrumental broadening, and get the crystallite
size from the Scherrer equation, tau = K * lambda / (beta * cos(theta)). With three or more peaks, a
Williamson-Hall plot separates size broadening from strain broadening.

Usage:
    python tool.py example/synthetic_xrd.csv --inst-fwhm 0.08 --williamson-hall -o xrd.png
    python tool.py pattern.csv --source cu-ka --range 30 80 --csv sizes.csv
    python tool.py --peak 38.5 0.45 --peak 44.7 0.47 --inst-fwhm 0.08       # peaks you measured yourself
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
from scipy.optimize import curve_fit
from scipy.signal import find_peaks, savgol_filter

K_DEFAULT = 0.9
MIN_POINTS = 10
RESULT_COLUMNS = ["two_theta_deg", "d_angstrom", "fwhm_obs_deg", "fwhm_corrected_deg", "size_nm", "size_unc_nm",
                  "eta", "r_squared"]

# Wavelengths in angstrom, as listed by the UNSW XRD Lab (Tian, "Choice of X-ray sources") and the UNT lecture
# notes on XRD (Golden): the K-alpha1 line and the K-alpha weighted average of the alpha1/alpha2 doublet.
WAVELENGTHS = {
    "cu-ka1": 1.54056, "cu-ka": 1.54184, "co-ka1": 1.78897, "co-ka": 1.79026,
    "cr-ka1": 2.28970, "cr-ka": 2.29100, "mo-ka1": 0.70930, "mo-ka": 0.71073,
}


@dataclass
class Pattern:
    x: np.ndarray  # 2-theta, degrees, ascending
    y: np.ndarray


# ---------- reading ----------
def read_pattern(path: Path, column: str | None = None) -> Pattern:
    """Read 2-theta (column 1) and an intensity column; header and metadata lines are skipped."""
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
        raise ValueError(f"{path}: no numeric data found (expected 2-theta and intensity columns)")
    data = np.array(rows)
    data = data[np.argsort(data[:, 0], kind="stable")]
    _, keep = np.unique(data[:, 0], return_index=True)
    data = data[np.sort(keep)]
    if data.shape[0] < 50:
        raise ValueError(f"{path}: only {data.shape[0]} data points, need at least 50")
    names = header[1:] if header and len(header) == ncols else [f"column {i}" for i in range(1, ncols)]
    index = 0
    if column is not None:
        if column in names:
            index = names.index(column)
        elif column.isdigit() and 1 <= int(column) <= ncols - 1:
            index = int(column) - 1
        else:
            raise ValueError(f"{path}: column {column!r} not found. Available: {', '.join(names)}")
    return Pattern(data[:, 0].copy(), data[:, index + 1].copy())


# ---------- the physics ----------
def bragg_d(two_theta_deg: float, wavelength: float) -> float:
    """Lattice spacing in angstrom from Bragg's law, d = lambda / (2 sin theta)."""
    return wavelength / (2.0 * np.sin(np.radians(two_theta_deg) / 2.0))


def correct_fwhm(observed: float, instrument: float, method: str = "gaussian") -> float | None:
    """Remove the instrumental broadening. None if the peak is not wider than the instrument."""
    if method == "none" or instrument <= 0:
        return observed
    if observed <= instrument:
        return None
    return float(np.sqrt(observed**2 - instrument**2)) if method == "gaussian" else observed - instrument


def scherrer_size_nm(two_theta_deg: float, fwhm_deg: float, wavelength: float, k: float = K_DEFAULT) -> float:
    """Scherrer crystallite size in nm. `fwhm_deg` is the corrected width in degrees 2-theta."""
    beta = np.radians(fwhm_deg)
    theta = np.radians(two_theta_deg) / 2.0
    return float(k * wavelength / (beta * np.cos(theta)) / 10.0)


def analyze_peak(two_theta: float, fwhm_obs: float, wavelength: float, *, k: float, inst: float, method: str,
                 fwhm_unc: float | None = None) -> dict:
    corrected = correct_fwhm(fwhm_obs, inst, method)
    row = {"two_theta_deg": two_theta, "d_angstrom": bragg_d(two_theta, wavelength), "fwhm_obs_deg": fwhm_obs,
           "fwhm_corrected_deg": corrected, "size_nm": None, "size_unc_nm": None}
    if corrected is not None:
        row["size_nm"] = scherrer_size_nm(two_theta, corrected, wavelength, k)
        if fwhm_unc is not None:
            factor = fwhm_obs / corrected if method == "gaussian" and inst > 0 else 1.0
            row["size_unc_nm"] = row["size_nm"] * fwhm_unc * factor / corrected
    return row


def williamson_hall(rows: list[dict], wavelength: float, k: float = K_DEFAULT) -> dict:
    """Fit beta*cos(theta) = K*lambda/D + 4*strain*sin(theta) through the peaks with a corrected width."""
    use = [r for r in rows if r["fwhm_corrected_deg"] is not None]
    if len(use) < 3:
        raise ValueError("a Williamson-Hall fit needs at least 3 peaks wider than the instrumental broadening")
    theta = np.radians([r["two_theta_deg"] for r in use]) / 2.0
    x = 4.0 * np.sin(theta)
    y = np.radians([r["fwhm_corrected_deg"] for r in use]) * np.cos(theta)
    n = x.size
    sxx = float(np.sum((x - x.mean()) ** 2))
    slope = float(np.sum((x - x.mean()) * (y - y.mean())) / sxx)
    intercept = float(y.mean() - slope * x.mean())
    resid = y - (slope * x + intercept)
    s2 = float(np.sum(resid**2)) / (n - 2)
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    sigma_intercept = float(np.sqrt(s2 * (1.0 / n + x.mean() ** 2 / sxx)))
    size = k * wavelength / intercept / 10.0 if intercept > 0 else None
    return {"strain": slope, "strain_unc": float(np.sqrt(s2 / sxx)), "intercept": intercept, "size_nm": size,
            "size_unc_nm": None if size is None else size * sigma_intercept / intercept,
            "r2": 1.0 - float(np.sum(resid**2)) / ss_tot if ss_tot > 0 else None, "n": n, "x": x, "y": y}


# ---------- fitting the peaks ----------
def pseudo_voigt(x, amp, x0, w, eta, b0, b1):
    """Pseudo-Voigt peak (one FWHM `w` for the Lorentzian and Gaussian parts) on a linear background."""
    u = (x - x0) / w
    return amp * (eta / (1 + 4 * u**2) + (1 - eta) * np.exp(-4 * np.log(2) * u**2)) + b0 + b1 * (x - x0)


def auto_window(p: Pattern, width_deg: float = 0.25) -> int:
    """Odd smoothing window (at least 5 points) spanning about `width_deg` degrees."""
    step = float(np.median(np.diff(p.x)))
    window = max(5, round(width_deg / step))
    return window + 1 if window % 2 == 0 else window


def detect_peaks(p: Pattern, prominence: float = 0.05, smooth: int | None = None, min_separation: float = 0.5) -> list[float]:
    """Positions of the peaks. The pattern is smoothed for detection only (default: about 0.25 degrees), because
    on raw counting noise a random spike can look prominent next to the deepest noise dips."""
    window = auto_window(p) if smooth is None else smooth
    y = savgol_filter(p.y, window, 3) if window else p.y
    step = float(np.median(np.diff(p.x)))
    idx, _ = find_peaks(y, prominence=max(prominence * float(y.max() - y.min()), 1e-12),
                        distance=max(1, round(min_separation / step)))
    return [float(p.x[i]) for i in idx]


def fit_peak(p: Pattern, center: float, lo: float, hi: float) -> dict:
    """Pseudo-Voigt fit of the peak near `center` using only the points between lo and hi."""
    mask = (p.x >= lo) & (p.x <= hi)
    if mask.sum() < MIN_POINTS:
        raise ValueError(f"too few points around {center:.2f} degrees")
    x, y = p.x[mask], p.y[mask]
    base = float(np.median(np.concatenate([y[:3], y[-3:]])))
    amp0 = max(float(y.max()) - base, 1e-9)
    above = x[y >= base + amp0 / 2]
    w0 = max(float(above.max() - above.min()), 2 * float(np.median(np.diff(x))))
    guess = [amp0, center, w0, 0.5, base, 0.0]
    bounds = ([0, center - w0, float(np.median(np.diff(x))), 0.0, -np.inf, -np.inf],
              [np.inf, center + w0, hi - lo, 1.0, np.inf, np.inf])
    try:
        popt, pcov = curve_fit(pseudo_voigt, x, y, p0=guess, bounds=bounds, maxfev=20000)
    except (RuntimeError, ValueError) as err:
        raise ValueError(f"the fit near {center:.2f} degrees did not converge") from err
    err = np.sqrt(np.clip(np.diag(pcov), 0, None))
    # R2 only close to the peak (within 1.5 FWHM), so a wide, mostly flat window does not hide a good fit
    near = np.abs(x - popt[1]) <= 1.5 * popt[2]
    resid = y[near] - pseudo_voigt(x[near], *popt)
    ss_tot = float(np.sum((y[near] - y[near].mean()) ** 2)) if near.sum() > 2 else 0.0
    return {"amp": float(popt[0]), "x0": float(popt[1]), "fwhm": float(popt[2]), "eta": float(popt[3]),
            "fwhm_unc": float(err[2]), "x0_unc": float(err[1]), "popt": popt, "window": (float(x.min()), float(x.max())),
            "r2": 1.0 - float(np.sum(resid**2)) / ss_tot if ss_tot > 0 else None}


def fit_all(p: Pattern, centers: list[float], max_half_window: float = 3.0) -> list[dict]:
    """Fit every peak in a window that stops halfway to the neighbouring peaks."""
    centers = sorted(centers)
    fits = []
    for i, c in enumerate(centers):
        half = max_half_window
        lo = c - half if i == 0 else max(c - half, (centers[i - 1] + c) / 2)
        hi = c + half if i == len(centers) - 1 else min(c + half, (centers[i + 1] + c) / 2)
        try:
            fits.append(fit_peak(p, c, lo, hi))
        except ValueError as err:
            print(f"[warn] {err}; peak skipped", file=sys.stderr)
    return fits


def size_warnings(rows: list[dict], inst: float, fits: list[dict] | None) -> list[str]:
    notes = []
    for r in rows:
        tt = f"{r['two_theta_deg']:.2f}"
        if r["fwhm_corrected_deg"] is None:
            notes.append(f"peak at {tt} degrees is not wider than the instrumental broadening: no size can be given")
        elif r["size_nm"] is not None and r["size_nm"] > 100:
            notes.append(f"peak at {tt} degrees gives {r['size_nm']:.0f} nm: Scherrer is unreliable above about "
                         "100 nm, the broadening is close to the instrument's")
        if inst > 0 and r["fwhm_corrected_deg"] is not None and r["fwhm_obs_deg"] < 1.5 * inst:
            notes.append(f"peak at {tt} degrees is only {r['fwhm_obs_deg'] / inst:.1f} times the instrumental width: "
                         "the size is very sensitive to the instrument value")
    for f in fits or []:
        if f["r2"] is not None and f["r2"] < 0.8:
            notes.append(f"the fit of the peak at {f['x0']:.2f} degrees is poor (R2 = {f['r2']:.3f}): check for "
                         "overlapping peaks or a noisy pattern")
    return notes


# ---------- output ----------
def plot_xrd(p: Pattern, fits: list[dict], rows: list[dict], wh: dict | None, path: Path, *, title: str | None,
             dpi: int, figsize: tuple[float, float] | None) -> None:
    both = wh is not None
    fig, axes = plt.subplots(1, 2 if both else 1, figsize=figsize or ((10.5, 4.0) if both else (7.0, 4.0)),
                             squeeze=False, gridspec_kw={"width_ratios": [1.7, 1]} if both else None)
    ax = axes[0][0]
    ax.plot(p.x, p.y, color="#555555", lw=0.8, label="data")
    for k, f in enumerate(fits):
        lo, hi = f["window"]
        xs = np.linspace(lo, hi, 300)
        ax.plot(xs, pseudo_voigt(xs, *f["popt"]), color="#D55E00", lw=1.2, label="fit" if k == 0 else None)
    for r in rows:
        if r["size_nm"] is not None:
            top = float(np.interp(r["two_theta_deg"], p.x, p.y))
            ax.annotate(f"{r['size_nm']:.0f} nm", (r["two_theta_deg"], top), textcoords="offset points", xytext=(0, 6),
                        ha="center", fontsize=8, color="#D55E00")
    ax.set_xlabel("2$\\theta$ (degrees)")
    ax.set_ylabel("Intensity (counts)")
    ax.legend(frameon=False, fontsize=8)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    if both:
        ax2 = axes[0][1]
        ax2.plot(wh["x"], wh["y"], "o", color="#0072B2")
        xs = np.array([0.0, wh["x"].max() * 1.05])
        ax2.plot(xs, wh["strain"] * xs + wh["intercept"], "--", color="black", lw=1.1)
        size_text = f"D = {wh['size_nm']:.1f} nm, " if wh["size_nm"] else ""
        ax2.set_title(f"{size_text}strain = {wh['strain']:.4f}", fontsize=9)
        ax2.set_xlabel("4 sin $\\theta$")
        ax2.set_ylabel("$\\beta \\cos\\theta$ (rad)")
        ax2.set_xlim(left=0)
        ax2.set_ylim(bottom=0)
        for side in ("top", "right"):
            ax2.spines[side].set_visible(False)
    if title:
        fig.suptitle(title)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi)
    plt.close(fig)


def result_table(rows: list[dict]) -> str:
    lines = [f"{'2theta':>7}  {'d (A)':>7}  {'FWHM obs':>9}  {'FWHM corr':>9}  {'size (nm)':>9}  {'+/- (nm)':>8}"]
    for r in rows:
        corr = "n/a" if r["fwhm_corrected_deg"] is None else f"{r['fwhm_corrected_deg']:.4f}"
        size = "n/a" if r["size_nm"] is None else f"{r['size_nm']:.1f}"
        unc = "" if r["size_unc_nm"] is None else f"{r['size_unc_nm']:.1f}"
        lines.append(f"{r['two_theta_deg']:>7.2f}  {r['d_angstrom']:>7.4f}  {r['fwhm_obs_deg']:>9.4f}  {corr:>9}  "
                     f"{size:>9}  {unc:>8}")
    return "\n".join(lines)


def write_csv(path: Path, rows: list[dict], fits: list[dict] | None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=RESULT_COLUMNS, lineterminator="\n")
        writer.writeheader()
        for k, r in enumerate(rows):
            fit = fits[k] if fits else None
            writer.writerow({
                "two_theta_deg": round(r["two_theta_deg"], 4), "d_angstrom": round(r["d_angstrom"], 5),
                "fwhm_obs_deg": round(r["fwhm_obs_deg"], 5),
                "fwhm_corrected_deg": "" if r["fwhm_corrected_deg"] is None else round(r["fwhm_corrected_deg"], 5),
                "size_nm": "" if r["size_nm"] is None else round(r["size_nm"], 2),
                "size_unc_nm": "" if r["size_unc_nm"] is None else round(r["size_unc_nm"], 2),
                "eta": "" if fit is None else round(fit["eta"], 3),
                "r_squared": "" if fit is None or fit["r2"] is None else round(fit["r2"], 5)})


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", nargs="?", help="pattern file: 2-theta (degrees) in column 1, intensity in column 2")
    ap.add_argument("--peak", type=float, nargs=2, action="append", metavar=("2THETA", "FWHM"),
                    help="a peak you measured yourself: position and FWHM in degrees (repeat for several peaks)")
    ap.add_argument("--source", choices=sorted(WAVELENGTHS), default="cu-ka1", help="X-ray line (default cu-ka1)")
    ap.add_argument("--wavelength", type=float, help="wavelength in angstrom (overrides --source)")
    ap.add_argument("--k", type=float, default=K_DEFAULT, help="Scherrer shape factor (default 0.9)")
    ap.add_argument("--inst-fwhm", type=float, default=0.0, help="instrumental FWHM in degrees 2-theta (default 0)")
    ap.add_argument("--correction", choices=["gaussian", "lorentzian", "none"], default="gaussian",
                    help="how the instrumental width is removed (default gaussian)")
    ap.add_argument("--range", type=float, nargs=2, metavar=("LOW", "HIGH"), help="use only this 2-theta range")
    ap.add_argument("--peak-positions", type=float, nargs="+", help="fit these 2-theta positions instead of detecting peaks")
    ap.add_argument("--prominence", type=float, default=0.05, help="peak detection threshold, fraction of the intensity range")
    ap.add_argument("--smooth", type=int, default=None, metavar="N",
                    help="smoothing window used only to detect peaks (default: about 0.25 degrees, 0 turns it off)")
    ap.add_argument("--min-separation", type=float, default=0.5, help="minimum distance between peaks in degrees")
    ap.add_argument("--window", type=float, default=3.0,
                    help="half-width in degrees of the fit window around each peak (default 3; widen it for very broad peaks)")
    ap.add_argument("--column", help="intensity column (header name or number)")
    ap.add_argument("--williamson-hall", action="store_true", help="also fit a Williamson-Hall plot (3 or more peaks)")
    ap.add_argument("-o", "--output", help="save the pattern with fitted peaks (.png, .svg, .pdf)")
    ap.add_argument("--csv", help="save the per-peak results to a CSV")
    ap.add_argument("--title", help="plot title")
    ap.add_argument("--figsize", type=float, nargs=2, metavar=("W", "H"), help="figure size in inches")
    ap.add_argument("--dpi", type=int, default=300, help="PNG resolution (default 300)")
    args = ap.parse_args()

    try:
        if bool(args.input) == bool(args.peak):
            raise ValueError("give either a pattern file or one or more --peak 2THETA FWHM")
        if args.output and Path(args.output).suffix.lower() not in (".png", ".svg", ".pdf"):
            raise ValueError("output must end in .png, .svg or .pdf")
        if args.output and args.peak:
            raise ValueError("-o needs a pattern file")
        wavelength = args.wavelength if args.wavelength is not None else WAVELENGTHS[args.source]
        if wavelength <= 0 or args.k <= 0 or args.inst_fwhm < 0:
            raise ValueError("--wavelength and --k must be positive and --inst-fwhm must not be negative")
        pattern, fits = None, None
        if args.peak:
            if any(tt <= 0 or tt >= 180 or w <= 0 for tt, w in args.peak):
                raise ValueError("peak positions must be between 0 and 180 degrees and widths must be positive")
            rows = [analyze_peak(tt, w, wavelength, k=args.k, inst=args.inst_fwhm, method=args.correction)
                    for tt, w in sorted(args.peak)]
        else:
            pattern = read_pattern(Path(args.input), args.column)
            if args.range:
                lo, hi = args.range
                if lo >= hi:
                    raise ValueError("--range needs LOW below HIGH")
                keep = (pattern.x >= lo) & (pattern.x <= hi)
                if keep.sum() < 50:
                    raise ValueError("fewer than 50 points in the selected range")
                pattern = Pattern(pattern.x[keep], pattern.y[keep])
            if args.smooth and (args.smooth % 2 == 0 or args.smooth < 5 or args.smooth > pattern.y.size):
                raise ValueError("--smooth needs an odd window of at least 5 points, shorter than the data")
            centers = args.peak_positions or detect_peaks(pattern, args.prominence, args.smooth, args.min_separation)
            if not centers:
                raise ValueError("no peaks found; try a smaller --prominence or give --peak-positions")
            if args.window <= 0:
                raise ValueError("--window must be positive")
            fits = fit_all(pattern, centers, args.window)
            if not fits:
                raise ValueError("no peak could be fitted")
            rows = [analyze_peak(f["x0"], f["fwhm"], wavelength, k=args.k, inst=args.inst_fwhm, method=args.correction,
                                 fwhm_unc=f["fwhm_unc"]) for f in fits]
        wh = williamson_hall(rows, wavelength, args.k) if args.williamson_hall else None
        if args.output:
            plot_xrd(pattern, fits, rows, wh, Path(args.output), title=args.title, dpi=args.dpi,
                     figsize=tuple(args.figsize) if args.figsize else None)
    except ValueError as err:
        print(f"error: {err}", file=sys.stderr)
        return 2

    print(f"wavelength {wavelength:g} A, K = {args.k:g}, instrumental FWHM {args.inst_fwhm:g} deg ({args.correction})")
    print(result_table(rows))
    sizes = [r["size_nm"] for r in rows if r["size_nm"] is not None]
    if len(sizes) >= 2:
        print(f"Scherrer size over {len(sizes)} peaks: mean {np.mean(sizes):.1f} nm, SD {np.std(sizes, ddof=1):.1f} nm")
    if wh:
        size = f"{wh['size_nm']:.1f} +/- {wh['size_unc_nm']:.1f} nm" if wh["size_nm"] else "undefined (intercept not positive)"
        r2 = "" if wh["r2"] is None else f", R2 {wh['r2']:.3f}"
        print(f"Williamson-Hall ({wh['n']} peaks): size {size}, strain {wh['strain']:.5f} +/- {wh['strain_unc']:.5f}{r2}")
        if wh["strain"] < 0:
            print("[warn] negative strain is not physical: the peaks may be dominated by size broadening or the "
                  "widths are noisy", file=sys.stderr)
    for note in size_warnings(rows, args.inst_fwhm, fits):
        print(f"[warn] {note}", file=sys.stderr)
    if args.output:
        print(f"[ok] plot -> {args.output}")
    if args.csv:
        write_csv(Path(args.csv), rows, fits)
        print(f"Wrote {len(rows)} peaks to {args.csv}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except BrokenPipeError:  # e.g. piping into `head`
        sys.exit(0)
