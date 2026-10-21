"""Day 19: baseline correction and smoothing toolkit for spectra.

Smooth a spectrum and remove its baseline, with the usual methods, and see what they did: how much noise
went away, how much the peaks were flattened, and how well a known baseline was recovered.

    Smoothing:  savgol, moving-average, gaussian, median, whittaker
    Baselines:  poly (modified polynomial fit), asls (asymmetric least squares), snip (peak clipping), anchors

Usage:
    python tool.py example/synthetic_spectrum.csv --baseline asls --baseline-lambda 1e7 -o corrected.png
    python tool.py spectrum.csv --smooth savgol --window 11 --baseline poly --degree 4 --csv corrected.csv
    python tool.py spectrum.csv --baseline snip --width 80 --smooth gaussian --sigma 3
    python tool.py spectrum.csv --baseline anchors --anchors 420 700 1130 1800 --anchor-width 20
"""
import argparse
import csv
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import sparse
from scipy.interpolate import CubicSpline
from scipy.ndimage import gaussian_filter1d, median_filter, uniform_filter1d
from scipy.signal import find_peaks, peak_prominences, savgol_filter
from scipy.sparse.linalg import spsolve

MIN_POINTS = 10
SMOOTHERS = ["none", "savgol", "moving-average", "gaussian", "median", "whittaker"]
BASELINES = ["none", "poly", "asls", "snip", "anchors"]
CSV_COLUMNS = ["x", "raw", "smoothed", "baseline", "corrected"]


@dataclass
class Spectrum:
    name: str
    x: np.ndarray  # ascending
    y: np.ndarray


# ---------- reading ----------
def read_spectrum(path: Path, column: str | None = None) -> Spectrum:
    """Read x in column 1 and a signal column; header and metadata lines are skipped, x is sorted ascending."""
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
        raise ValueError(f"{path}: no numeric data found (expected x and signal columns)")
    data = np.array(rows)
    data = data[np.argsort(data[:, 0], kind="stable")]
    _, keep = np.unique(data[:, 0], return_index=True)
    data = data[np.sort(keep)]
    if data.shape[0] < MIN_POINTS:
        raise ValueError(f"{path}: only {data.shape[0]} data points, need at least {MIN_POINTS}")
    names = header[1:] if header and len(header) == ncols else [f"column {i}" for i in range(1, ncols)]
    index = 0
    if column is not None:
        if column in names:
            index = names.index(column)
        elif column.isdigit() and 1 <= int(column) <= ncols - 1:
            index = int(column) - 1
        else:
            raise ValueError(f"{path}: column {column!r} not found. Available: {', '.join(names)}")
    return Spectrum(names[index], data[:, 0].copy(), data[:, index + 1].copy())


def step_of(x: np.ndarray) -> float:
    return float(np.median(np.diff(x)))


def spacing_uniformity(x: np.ndarray) -> float:
    """Largest relative deviation of the spacing from its median (0 for evenly spaced data)."""
    d = np.diff(x)
    return float(np.max(np.abs(d - np.median(d))) / np.median(d))


def noise_sigma(y: np.ndarray) -> float:
    """Robust noise level from the point-to-point differences (not affected by smooth baselines or peaks)."""
    d = np.diff(y)
    return float(1.4826 * np.median(np.abs(d - np.median(d))) / np.sqrt(2.0))


# ---------- smoothing ----------
def penalty_matrix(n: int, order: int = 2) -> sparse.csc_matrix:
    """D'D for the difference matrix D of the given order (the roughness penalty of the Whittaker smoother)."""
    d = sparse.eye(n, format="csr")
    for _ in range(order):
        d = d[1:] - d[:-1]
    return (d.T @ d).tocsc()


def whittaker_smooth(y: np.ndarray, lam: float, weights: np.ndarray | None = None, order: int = 2) -> np.ndarray:
    """Penalized least squares smoother (Eilers 2003): minimizes sum(w (y - z)^2) + lam sum((diff^order z)^2)."""
    if lam < 0:
        raise ValueError("the smoothing parameter must not be negative")
    if y.size < order + 2:
        raise ValueError(f"need at least {order + 2} points to smooth")
    w = np.ones_like(y) if weights is None else weights
    return spsolve((sparse.diags(w, format="csc") + lam * penalty_matrix(y.size, order)).tocsc(), w * y)


def smooth(y: np.ndarray, method: str, *, window: int = 11, order: int = 3, sigma_points: float = 2.0,
           lam: float = 100.0) -> np.ndarray:
    if method == "none":
        return y.copy()
    if method in ("savgol", "median") and (window % 2 == 0 or window < 3):
        raise ValueError(f"--window must be odd and at least 3 for {method}")
    if method == "savgol":
        if window <= order or window > y.size:
            raise ValueError(f"--window ({window}) must be larger than --order ({order}) and not longer than the data")
        return savgol_filter(y, window, order)
    if method == "moving-average":
        if window < 1 or window > y.size:
            raise ValueError("--window must be between 1 and the number of points")
        return uniform_filter1d(y, size=window, mode="nearest")
    if method == "gaussian":
        if sigma_points <= 0:
            raise ValueError("--sigma must be positive")
        return gaussian_filter1d(y, sigma_points, mode="nearest")
    if method == "median":
        return median_filter(y, size=window, mode="nearest")
    if method == "whittaker":
        return whittaker_smooth(y, lam)
    raise ValueError(f"unknown smoothing method {method!r}")


# ---------- baselines ----------
def modpoly(x: np.ndarray, y: np.ndarray, degree: int = 4, niter: int = 100, tol: float = 1e-3) -> tuple[np.ndarray, int]:
    """Modified polynomial fit (Lieber and Mahadevan-Jansen 2003): fit a polynomial, keep the lower of the data and
    the fit at every point, and repeat until the fit stops changing. Peaks are shaved off, so the fit settles under them."""
    if degree < 0 or degree >= y.size:
        raise ValueError("--degree must be at least 0 and below the number of points")
    work, previous = y.copy(), None
    for iteration in range(1, niter + 1):
        fit = np.polynomial.Polynomial.fit(x, work, degree)(x)
        if previous is not None and np.linalg.norm(fit - previous) <= tol * max(np.linalg.norm(previous), 1e-300):
            return fit, iteration
        previous = fit
        work = np.minimum(work, fit)
    return previous, niter


def asls(y: np.ndarray, lam: float = 1e5, p: float = 0.01, niter: int = 20, tol: float = 1e-6) -> tuple[np.ndarray, int]:
    """Asymmetric least squares baseline (Eilers and Boelens 2005). Points above the baseline get weight p and points
    below get 1 - p, so the smooth curve settles under the peaks."""
    if not 0 < p < 1:
        raise ValueError("--asymmetry must be between 0 and 1")
    if lam <= 0:
        raise ValueError("--baseline-lambda must be positive")
    w = np.ones_like(y)
    z = y
    for iteration in range(1, niter + 1):
        z = whittaker_smooth(y, lam, w)
        w_new = np.where(y > z, p, 1.0 - p)
        if np.linalg.norm(w_new - w) <= tol * np.sqrt(y.size):
            return z, iteration
        w = w_new
    return z, niter


def snip(y: np.ndarray, half_window: int, lls: bool = True) -> np.ndarray:
    """SNIP peak clipping (Ryan et al. 1988, Morhac et al. 1997) with an increasing clipping window, optionally on
    the log-log-square transformed data so that peaks of very different heights are clipped alike."""
    if half_window < 1:
        raise ValueError("the SNIP half-window must be at least 1 point")
    shift = float(y.min())
    v = y - shift
    if lls:
        v = np.log(np.log(np.sqrt(v + 1.0) + 1.0) + 1.0)
    v = np.pad(v, half_window, mode="edge")
    for p in range(1, half_window + 1):
        previous = v.copy()
        v[p:-p] = np.minimum(previous[p:-p], 0.5 * (previous[:-2 * p] + previous[2 * p:]))
    v = v[half_window:-half_window]
    if lls:
        v = (np.exp(np.exp(v) - 1.0) - 1.0) ** 2 - 1.0
    return v + shift


def anchor_baseline(x: np.ndarray, y: np.ndarray, anchors: list[float], width: float, kind: str = "linear") -> np.ndarray:
    """Baseline through the medians of the data around chosen peak-free positions."""
    if len(anchors) < 2 or (kind == "cubic" and len(anchors) < 4):
        raise ValueError("--anchors needs at least 2 positions (4 for a cubic curve)")
    if width <= 0:
        raise ValueError("--anchor-width must be positive")
    ax, ay = [], []
    for a in sorted(anchors):
        inside = np.abs(x - a) <= width / 2
        if not inside.any():
            raise ValueError(f"no data within {width / 2:g} of the anchor at {a:g}")
        ax.append(a)
        ay.append(float(np.median(y[inside])))
    if kind == "cubic":
        curve = CubicSpline(ax, ay, extrapolate=False)(x)
        return np.where(np.isnan(curve), np.where(x < ax[0], ay[0], ay[-1]), curve)
    return np.interp(x, ax, ay)  # constant beyond the first and last anchor


def estimate_baseline(x: np.ndarray, y: np.ndarray, method: str, args) -> tuple[np.ndarray, str]:
    """Baseline of the spectrum and a short description of what was done."""
    if method == "none":
        return np.zeros_like(y), "none"
    if method == "poly":
        base, n = modpoly(x, y, args.degree)
        return base, f"modified polynomial, degree {args.degree}, {n} iterations"
    if method == "asls":
        base, n = asls(y, args.baseline_lambda, args.asymmetry)
        return base, f"asymmetric least squares, lambda {args.baseline_lambda:g}, p {args.asymmetry:g}, {n} iterations"
    if method == "snip":
        if args.width is None:
            raise ValueError("--baseline snip needs --width (the width of the widest peak, in x units)")
        half = max(1, round((args.width / step_of(x) - 1) / 2))
        return snip(y, half, not args.no_lls), f"SNIP, clipping half-window {half} points" + ("" if not args.no_lls else ", no LLS transform")
    if method == "anchors":
        if not args.anchors:
            raise ValueError("--baseline anchors needs --anchors (peak-free positions, in x units)")
        return anchor_baseline(x, y, args.anchors, args.anchor_width, args.anchor_kind), f"{args.anchor_kind} through {len(args.anchors)} anchors"
    raise ValueError(f"unknown baseline method {method!r}")


# ---------- the whole pipeline ----------
def process(s: Spectrum, args) -> dict:
    """Smooth, estimate the baseline and subtract it. The baseline is estimated on the smoothed spectrum."""
    step = step_of(s.x)
    sigma_points = args.sigma / step if args.sigma is not None else 2.0
    smoothed = smooth(s.y, args.smooth, window=args.window, order=args.order, sigma_points=sigma_points,
                      lam=args.smooth_lambda)
    source = s.y if args.baseline_from == "raw" else smoothed
    baseline, description = estimate_baseline(s.x, source, args.baseline, args)
    corrected = smoothed - baseline
    unsmoothed_peak = float(np.max(s.y - baseline))
    sigma = noise_sigma(s.y)
    return {"smoothed": smoothed, "baseline": baseline, "corrected": corrected, "description": description,
            "noise_raw": sigma, "noise_smoothed": noise_sigma(smoothed),
            "peak_retained": float(np.max(corrected) / unsmoothed_peak) if unsmoothed_peak > 0 else None,
            "prominence_kept": prominence_kept(smoothed, corrected, sigma) if args.baseline != "none" else None}


def prominence_kept(before: np.ndarray, after: np.ndarray, sigma: float) -> tuple[float, int] | None:
    """Median fraction of each clear peak's prominence that is still there after the baseline was removed.

    About 1 means the baseline stayed under the peaks. A small value means the baseline followed the peaks and
    took them out. Returns None when the spectrum has no peak standing 10 noise levels above its surroundings."""
    if sigma <= 0:
        return None
    idx, props = find_peaks(before, prominence=10 * sigma)
    if idx.size == 0:
        return None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # scipy warns when a peak is gone (prominence 0), which is the answer here
        remaining = peak_prominences(after, idx)[0]
    return float(np.median(remaining / props["prominences"])), int(idx.size)


def baseline_error(x: np.ndarray, baseline: np.ndarray, truth_x: np.ndarray, truth: np.ndarray) -> dict:
    """Difference between the estimated baseline and a known one (interpolated onto the spectrum's x)."""
    reference = np.interp(x, truth_x, truth)
    diff = baseline - reference
    return {"rms": float(np.sqrt(np.mean(diff**2))), "max": float(np.max(np.abs(diff))), "mean": float(np.mean(diff))}


def warnings_for(s: Spectrum, result: dict, args) -> list[str]:
    notes = []
    if spacing_uniformity(s.x) > 0.01:
        notes.append("the x spacing is not uniform: window sizes in points then mean different widths in x, and the "
                     "asymmetric least squares and Whittaker methods assume even spacing")
    retained = result["peak_retained"]
    if args.smooth != "none" and retained is not None and retained < 0.95:
        notes.append(f"smoothing lowered the highest peak to {100 * retained:.1f}% of its height: the window is too wide "
                     "for the peaks (and they will also be broader)")
    sigma = result["noise_raw"]
    below = float(np.mean(result["corrected"] < -3 * sigma)) if sigma > 0 else 0.0
    if args.baseline != "none" and below > 0.05:
        notes.append(f"{100 * below:.0f}% of the corrected points lie more than 3 noise levels below zero: the baseline "
                     "runs above the data in places (too flexible, or not asymmetric enough)")
    kept = result["prominence_kept"]
    if kept is not None and kept[0] < 0.7:
        notes.append(f"the baseline removed {100 * (1 - kept[0]):.0f}% of the height of the peaks (median over {kept[1]} "
                     "peaks): it is following the peaks. Make it stiffer (larger lambda or a lower degree) or "
                     "more asymmetric (smaller p), or widen the SNIP width")
    if args.baseline == "poly" and args.degree > 8:
        notes.append("a high polynomial degree starts to follow the peaks as well as the baseline")
    return notes


# ---------- output ----------
def plot_result(s: Spectrum, result: dict, path: Path, *, truth: tuple | None, smoothed: bool, reverse_x: bool,
                xlabel: str, ylabel: str, title: str | None, figsize: tuple[float, float] | None, dpi: int) -> None:
    fig, (top, bottom) = plt.subplots(2, 1, sharex=True, figsize=figsize or (7.0, 6.0))
    top.plot(s.x, s.y, color="#999999", lw=0.8, label="raw")
    if smoothed:
        top.plot(s.x, result["smoothed"], color="#0072B2", lw=1.1, label="smoothed")
    top.plot(s.x, result["baseline"], color="#D55E00", lw=1.8, label="estimated baseline")
    if truth is not None:
        top.plot(truth[0], truth[1], color="black", lw=1.0, ls="--", label="true baseline")
    top.set_ylabel(ylabel)
    top.legend(frameon=False, fontsize=8)
    bottom.plot(s.x, result["corrected"], color="#009E73", lw=1.1)
    bottom.axhline(0, color="black", lw=0.6)
    bottom.set_ylabel(f"{ylabel} (corrected)")
    bottom.set_xlabel(xlabel)
    for ax in (top, bottom):
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    if reverse_x:
        bottom.set_xlim(s.x.max(), s.x.min())
    if title:
        fig.suptitle(title)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi)
    plt.close(fig)


def write_csv(path: Path, s: Spectrum, result: dict, corrected_only: bool = False) -> None:
    """All columns, or just x and the corrected spectrum (a file the other tools in this repo can read directly)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        if corrected_only:
            writer.writerow(["x", "corrected"])
            rows = zip(s.x, result["corrected"])
        else:
            writer.writerow(CSV_COLUMNS)
            rows = zip(s.x, s.y, result["smoothed"], result["baseline"], result["corrected"])
        for row in rows:
            writer.writerow([f"{v:.8g}" for v in row])


def read_truth(path: Path) -> tuple[np.ndarray, np.ndarray]:
    spec = read_spectrum(path)
    return spec.x, spec.y


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", help="spectrum file: x in column 1, signal in column 2")
    ap.add_argument("--column", help="signal column (header name or number)")
    ap.add_argument("--smooth", choices=SMOOTHERS, default="none", help="smoothing method (default none)")
    ap.add_argument("--window", type=int, default=11, help="window in points for savgol, moving-average and median (default 11)")
    ap.add_argument("--order", type=int, default=3, help="polynomial order for savgol (default 3)")
    ap.add_argument("--sigma", type=float, help="width of the gaussian smoother in x units")
    ap.add_argument("--smooth-lambda", type=float, default=100.0, help="smoothness of the whittaker smoother (default 100)")
    ap.add_argument("--baseline", choices=BASELINES, default="none", help="baseline method (default none)")
    ap.add_argument("--baseline-from", choices=["smoothed", "raw"], default="smoothed",
                    help="estimate the baseline on the smoothed (default) or the raw spectrum")
    ap.add_argument("--degree", type=int, default=4, help="polynomial degree for --baseline poly (default 4)")
    ap.add_argument("--baseline-lambda", type=float, default=1e5, help="stiffness for --baseline asls (default 1e5)")
    ap.add_argument("--asymmetry", type=float, default=0.01, help="p for --baseline asls (default 0.01)")
    ap.add_argument("--width", type=float, help="width of the widest peak in x units, for --baseline snip")
    ap.add_argument("--no-lls", action="store_true", help="SNIP without the log-log-square transform")
    ap.add_argument("--anchors", type=float, nargs="+", help="peak-free x positions for --baseline anchors")
    ap.add_argument("--anchor-width", type=float, default=10.0, help="x range averaged around each anchor (default 10)")
    ap.add_argument("--anchor-kind", choices=["linear", "cubic"], default="linear", help="curve through the anchors")
    ap.add_argument("--truth", help="file with the true baseline (x, baseline) to measure the error against")
    ap.add_argument("--xmin", type=float, help="lowest x to use")
    ap.add_argument("--xmax", type=float, help="highest x to use")
    ap.add_argument("-o", "--output", help="save the plot (.png, .svg or .pdf)")
    ap.add_argument("--csv", help="save x, raw, smoothed, baseline and corrected columns")
    ap.add_argument("--corrected-only", action="store_true", help="make the --csv file just x and corrected")
    ap.add_argument("--reverse-x", action="store_true", help="plot with x decreasing (usual for FTIR)")
    ap.add_argument("--xlabel", default="x", help="x axis label")
    ap.add_argument("--ylabel", default="Intensity", help="y axis label")
    ap.add_argument("--title", help="plot title")
    ap.add_argument("--figsize", type=float, nargs=2, metavar=("W", "H"), help="figure size in inches")
    ap.add_argument("--dpi", type=int, default=300, help="PNG resolution (default 300)")
    args = ap.parse_args()

    try:
        if args.output and Path(args.output).suffix.lower() not in (".png", ".svg", ".pdf"):
            raise ValueError("output must end in .png, .svg or .pdf")
        s = read_spectrum(Path(args.input), args.column)
        keep = np.ones_like(s.x, dtype=bool)
        if args.xmin is not None:
            keep &= s.x >= args.xmin
        if args.xmax is not None:
            keep &= s.x <= args.xmax
        if keep.sum() < MIN_POINTS:
            raise ValueError(f"fewer than {MIN_POINTS} points in the selected range")
        s = Spectrum(s.name, s.x[keep], s.y[keep])
        result = process(s, args)
        truth = read_truth(Path(args.truth)) if args.truth else None
        if args.output:
            plot_result(s, result, Path(args.output), truth=truth, smoothed=args.smooth != "none", reverse_x=args.reverse_x,
                        xlabel=args.xlabel, ylabel=args.ylabel, title=args.title,
                        figsize=tuple(args.figsize) if args.figsize else None, dpi=args.dpi)
    except ValueError as err:
        print(f"error: {err}", file=sys.stderr)
        return 2

    print(f"{s.y.size} points, x {s.x.min():g} to {s.x.max():g} (step {step_of(s.x):g}), noise level {result['noise_raw']:.4g}")
    if args.smooth != "none":
        kept = "" if result["peak_retained"] is None else f", highest peak kept at {100 * result['peak_retained']:.1f}%"
        print(f"smoothing: {args.smooth} -> noise level {result['noise_smoothed']:.4g}{kept}")
    if args.baseline != "none":
        base = result["baseline"]
        print(f"baseline: {result['description']} (range {base.min():.4g} to {base.max():.4g})")
        print(f"corrected: range {result['corrected'].min():.4g} to {result['corrected'].max():.4g}")
        if result["prominence_kept"] is not None:
            frac, count = result["prominence_kept"]
            print(f"peaks: {100 * frac:.0f}% of the peak prominence kept (median over {count} clear peaks)")
    if truth is not None:
        err = baseline_error(s.x, result["baseline"], truth[0], truth[1])
        print(f"against the true baseline: RMS error {err['rms']:.4g} ({err['rms'] / result['noise_raw']:.1f} noise levels), "
              f"largest error {err['max']:.4g}")
    for note in warnings_for(s, result, args):
        print(f"[warn] {note}", file=sys.stderr)
    if args.output:
        print(f"[ok] plot -> {args.output}")
    if args.csv:
        write_csv(Path(args.csv), s, result, args.corrected_only)
        print(f"Wrote {s.y.size} rows to {args.csv}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except BrokenPipeError:  # e.g. piping into `head`
        sys.exit(0)
