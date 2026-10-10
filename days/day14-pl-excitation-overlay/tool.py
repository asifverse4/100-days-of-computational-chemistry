"""Day 14: PL spectrum plotter with excitation-dependent emission overlay.

Plot emission spectra recorded at several excitation wavelengths, colored by excitation, and find
how the emission maximum shifts with excitation (a hallmark of many carbon dots). Scattered
excitation light (Rayleigh and second-order) can be masked so it is not mistaken for emission.

Usage:
    python tool.py example/synthetic_pl.csv -o pl.png --mask-scatter 15
    python tool.py ex340.csv ex360.csv ex380.csv --normalize --csv maxima.csv
    python tool.py spectra.csv --ex 300 320 340 360 --stack 0.5 --overlay-only
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
from scipy.signal import savgol_filter

MIN_POINTS = 5
NUMBER_RE = re.compile(r"[-+]?\d+(?:\.\d+)?")
PEAK_COLUMNS = ["excitation_nm", "emission_max_nm", "intensity", "fwhm_nm"]


@dataclass
class Raw:
    name: str  # column header, or the file name when there is none
    stem: str
    x: np.ndarray
    y: np.ndarray
    has_header: bool


@dataclass
class Spectrum:
    ex: float
    x: np.ndarray
    y: np.ndarray  # may contain NaN where scatter is masked
    label: str


# ---------- reading ----------
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


def read_raw(path: Path) -> list[Raw]:
    """Read a text file with emission wavelength in column 1 and one or more intensity columns."""
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
        raise ValueError(f"{path}: no numeric data found (expected emission wavelength and intensity columns)")
    data = np.array(rows)
    data = data[np.argsort(data[:, 0], kind="stable")]
    _, keep = np.unique(data[:, 0], return_index=True)
    data = data[np.sort(keep)]
    if data.shape[0] < MIN_POINTS * 2:
        raise ValueError(f"{path}: only {data.shape[0]} data points, need at least {MIN_POINTS * 2}")
    names = header[1:] if header and len(header) == ncols else []
    return [Raw(names[i - 1] if i - 1 < len(names) else path.stem, path.stem, data[:, 0].copy(), data[:, i].copy(),
                i - 1 < len(names)) for i in range(1, ncols)]


def number_in(text: str) -> float | None:
    found = NUMBER_RE.search(text)
    return float(found[0]) if found else None


def assign_excitation(raws: list[Raw], override: list[float] | None) -> list[Spectrum]:
    """Excitation wavelength from --ex, else from the column header, else from the file name."""
    if override is not None:
        if len(override) != len(raws):
            raise ValueError(f"--ex gives {len(override)} wavelengths but there are {len(raws)} spectra")
        values = list(override)
    else:
        values = []
        for r in raws:
            candidates = ([r.name] if r.has_header else []) + [r.stem]
            ex = next((v for v in (number_in(c) for c in candidates) if v is not None), None)
            if ex is None:
                raise ValueError(f"cannot find the excitation wavelength for {r.name!r}: put it in the column "
                                 "header or file name (for example ex340), or use --ex")
            values.append(ex)
    if len(set(values)) != len(values):
        raise ValueError("two spectra have the same excitation wavelength")
    pairs = sorted(zip(values, raws), key=lambda p: p[0])
    return [Spectrum(ex, r.x, r.y, f"{ex:g} nm") for ex, r in pairs]


# ---------- processing ----------
def crop(s: Spectrum, xmin: float | None, xmax: float | None) -> Spectrum:
    mask = np.ones_like(s.x, dtype=bool)
    if xmin is not None:
        mask &= s.x >= xmin
    if xmax is not None:
        mask &= s.x <= xmax
    if mask.sum() < MIN_POINTS:
        raise ValueError(f"{s.label}: fewer than {MIN_POINTS} points between {xmin} and {xmax} nm")
    return Spectrum(s.ex, s.x[mask], s.y[mask], s.label)


def smooth(s: Spectrum, window: int) -> Spectrum:
    if window % 2 == 0 or window < MIN_POINTS:
        raise ValueError("--smooth needs an odd window of at least 5 points")
    if window > s.y.size:
        raise ValueError(f"{s.label}: smoothing window {window} is longer than the data ({s.y.size} points)")
    return Spectrum(s.ex, s.x, savgol_filter(s.y, window, 3), s.label)


def mask_scatter(s: Spectrum, margin: float) -> Spectrum:
    """Set the Rayleigh (at ex) and second-order (at 2 ex) scatter regions to NaN."""
    y = s.y.copy()
    y[(np.abs(s.x - s.ex) <= margin) | (np.abs(s.x - 2 * s.ex) <= margin)] = np.nan
    return Spectrum(s.ex, s.x, y, s.label)


def normalize(s: Spectrum) -> Spectrum:
    top = np.nanmax(s.y)
    if not top > 0:
        raise ValueError(f"{s.label}: cannot normalize, the maximum is not positive")
    return Spectrum(s.ex, s.x, s.y / top, s.label)


def prepare(s: Spectrum, args) -> Spectrum:
    s = crop(s, args.xmin, args.xmax)
    if args.smooth:
        s = smooth(s, args.smooth)
    if args.mask_scatter > 0:
        s = mask_scatter(s, args.mask_scatter)
    if args.normalize:
        s = normalize(s)
    return s


# ---------- emission maximum ----------
def half_crossing(x: np.ndarray, y: np.ndarray, start: int, half: float, step: int) -> float | None:
    """Wavelength where y crosses `half` walking from index `start` in direction `step` (-1 left, +1 right)."""
    j = start
    while 0 <= j + step < y.size:
        a, b = y[j], y[j + step]
        if np.isnan(b):
            return None
        if b < half <= a:
            return float(x[j] + (a - half) / (a - b) * (x[j + step] - x[j]))
        j += step
    return None


def refine_peak(x: np.ndarray, y: np.ndarray, i: int, top_fraction: float = 0.9) -> tuple[float, float]:
    """Peak position and height from a parabola fitted to the top of the peak (points above 90% of the maximum).

    A fit over many points is much less noisy than a three-point parabola on a broad band. Falls back to the
    sampled maximum when there are too few points."""
    threshold = top_fraction * y[i]
    lo = hi = i
    while lo > 0 and not np.isnan(y[lo - 1]) and y[lo - 1] >= threshold:
        lo -= 1
    while hi < y.size - 1 and not np.isnan(y[hi + 1]) and y[hi + 1] >= threshold:
        hi += 1
    if hi - lo + 1 >= 5:
        a, b, c = np.polyfit(x[lo:hi + 1] - x[i], y[lo:hi + 1], 2)
        if a < 0:
            vertex = -b / (2 * a)
            if x[lo] - x[i] <= vertex <= x[hi] - x[i]:
                return float(x[i] + vertex), float(c - b * b / (4 * a))
    return float(x[i]), float(y[i])


def find_maximum(s: Spectrum, search_offset: float = 10.0) -> dict | None:
    """Emission maximum and FWHM, searching only above ex + search_offset."""
    region = (s.x >= s.ex + search_offset) & ~np.isnan(s.y)
    if region.sum() < 3:
        return None
    i = int(np.argmax(np.where(region, s.y, -np.inf)))
    peak_x, peak_y = refine_peak(s.x, s.y, i)
    half = peak_y / 2
    left, right = half_crossing(s.x, s.y, i, half, -1), half_crossing(s.x, s.y, i, half, +1)
    return {"excitation_nm": s.ex, "emission_max_nm": peak_x, "intensity": peak_y,
            "fwhm_nm": None if left is None or right is None else right - left}


def fit_shift(peaks: list[dict]) -> dict | None:
    """Straight-line fit of emission maximum against excitation wavelength (needs 3 or more spectra)."""
    if len(peaks) < 3:
        return None
    ex = np.array([p["excitation_nm"] for p in peaks])
    em = np.array([p["emission_max_nm"] for p in peaks])
    slope, intercept = np.polyfit(ex, em, 1)
    ss_tot = float(np.sum((em - em.mean()) ** 2))
    r2 = 1.0 - float(np.sum((em - (slope * ex + intercept)) ** 2)) / ss_tot if ss_tot > 0 else None
    return {"slope": float(slope), "intercept": float(intercept), "r2": r2,
            "shift_nm": float(em[-1] - em[0]), "ex_range_nm": float(ex[-1] - ex[0])}


def scatter_warnings(peaks: list[dict], margin: float = 15.0) -> list[str]:
    notes = []
    for p in peaks:
        ex, em = p["excitation_nm"], p["emission_max_nm"]
        if abs(em - ex) < margin or abs(em - 2 * ex) < margin:
            order = "second-order" if abs(em - 2 * ex) < margin else "Rayleigh"
            notes.append(f"ex {ex:g} nm: the maximum ({em:.1f} nm) is within {margin:g} nm of {order} scatter; "
                         "use --mask-scatter")
    return notes


# ---------- output ----------
def plot_pl(spectra: list[Spectrum], peaks: list[dict], fit: dict | None, path: Path, *, cmap: str, stack: float,
            overlay_only: bool, normalized: bool, title: str | None, figsize: tuple[float, float] | None,
            dpi: int) -> None:
    both = not overlay_only and len(peaks) >= 2
    fig, axes = plt.subplots(1, 2 if both else 1, figsize=figsize or ((10.0, 4.0) if both else (5.5, 4.0)),
                             squeeze=False, gridspec_kw={"width_ratios": [1.5, 1]} if both else None)
    ax = axes[0][0]
    ex_all = np.array([s.ex for s in spectra])
    norm = plt.Normalize(ex_all.min(), ex_all.max() if ex_all.max() > ex_all.min() else ex_all.min() + 1)
    colormap = matplotlib.colormaps[cmap]
    by_ex = {p["excitation_nm"]: p for p in peaks}
    for k, s in enumerate(spectra):
        color = colormap(norm(s.ex))
        offset = k * stack
        ax.plot(s.x, s.y + offset, color=color, lw=1.3)
        p = by_ex.get(s.ex)
        if p:
            ax.plot(p["emission_max_nm"], p["intensity"] + offset, "v", color=color, ms=4, mec="black", mew=0.4)
    ax.set_xlabel("Emission wavelength (nm)")
    ax.set_ylabel("Normalized PL intensity" if normalized else "PL intensity (a.u.)")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    bar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=colormap), ax=ax, pad=0.02)
    bar.set_label("Excitation wavelength (nm)")
    if both:
        ax2 = axes[0][1]
        ex = np.array([p["excitation_nm"] for p in peaks])
        em = np.array([p["emission_max_nm"] for p in peaks])
        ax2.scatter(ex, em, c=colormap(norm(ex)), edgecolors="black", linewidths=0.5, s=40, zorder=3)
        if fit:
            xs = np.array([ex.min(), ex.max()])
            ax2.plot(xs, fit["slope"] * xs + fit["intercept"], "--", color="black", lw=1.1)
            ax2.set_title(f"slope {fit['slope']:.2f} nm/nm, shift {fit['shift_nm']:+.0f} nm", fontsize=9)
        ax2.set_xlabel("Excitation wavelength (nm)")
        ax2.set_ylabel("Emission maximum (nm)")
        for side in ("top", "right"):
            ax2.spines[side].set_visible(False)
    if title:
        fig.suptitle(title)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi)
    plt.close(fig)


def peak_table(peaks: list[dict]) -> str:
    lines = [f"{'ex (nm)':>8}  {'em max (nm)':>11}  {'intensity':>9}  {'FWHM (nm)':>9}"]
    for p in peaks:
        fwhm = "n/a" if p["fwhm_nm"] is None else f"{p['fwhm_nm']:.1f}"
        lines.append(f"{p['excitation_nm']:>8g}  {p['emission_max_nm']:>11.1f}  {p['intensity']:>9.3f}  {fwhm:>9}")
    return "\n".join(lines)


def write_csv(path: Path, peaks: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=PEAK_COLUMNS, lineterminator="\n")
        writer.writeheader()
        for p in peaks:
            writer.writerow({"excitation_nm": p["excitation_nm"], "emission_max_nm": round(p["emission_max_nm"], 2),
                             "intensity": round(p["intensity"], 4),
                             "fwhm_nm": "" if p["fwhm_nm"] is None else round(p["fwhm_nm"], 2)})


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+", help="emission spectra: one file with a column per excitation, or one file each")
    ap.add_argument("-o", "--output", default="pl.png", help="plot file: .png, .svg or .pdf (default pl.png)")
    ap.add_argument("--csv", help="also write the emission maxima to this CSV")
    ap.add_argument("--ex", type=float, nargs="+", help="excitation wavelengths (nm), in the order of the spectra")
    ap.add_argument("--mask-scatter", type=float, default=0.0, metavar="NM",
                    help="blank out Rayleigh (ex) and second-order (2 ex) scatter within +/- NM (try 15)")
    ap.add_argument("--search-offset", type=float, default=10.0,
                    help="look for the maximum only above ex + this many nm (default 10)")
    ap.add_argument("--normalize", action="store_true", help="scale each spectrum to a maximum of 1")
    ap.add_argument("--smooth", type=int, default=0, metavar="N", help="Savitzky-Golay window (odd, >= 5 points)")
    ap.add_argument("--xmin", type=float, help="shortest emission wavelength to use (nm)")
    ap.add_argument("--xmax", type=float, help="longest emission wavelength to use (nm)")
    ap.add_argument("--stack", type=float, default=0.0, help="vertical offset between spectra")
    ap.add_argument("--cmap", default="viridis", help="matplotlib colormap for the excitation scale (default viridis)")
    ap.add_argument("--overlay-only", action="store_true", help="skip the emission-vs-excitation panel")
    ap.add_argument("--title", help="figure title")
    ap.add_argument("--figsize", type=float, nargs=2, metavar=("W", "H"), help="figure size in inches")
    ap.add_argument("--dpi", type=int, default=300, help="PNG resolution (default 300)")
    args = ap.parse_args()

    out = Path(args.output)
    if out.suffix.lower() not in (".png", ".svg", ".pdf"):
        print("error: output must end in .png, .svg or .pdf", file=sys.stderr)
        return 2
    if args.cmap not in matplotlib.colormaps:
        print(f"error: unknown colormap {args.cmap!r}", file=sys.stderr)
        return 2
    try:
        raws = [r for p in args.inputs for r in read_raw(Path(p))]
        spectra = [prepare(s, args) for s in assign_excitation(raws, args.ex)]
        peaks = [p for p in (find_maximum(s, args.search_offset) for s in spectra) if p]
        fit = fit_shift(peaks)
        plot_pl(spectra, peaks, fit, out, cmap=args.cmap, stack=args.stack, overlay_only=args.overlay_only,
                normalized=args.normalize, title=args.title, figsize=tuple(args.figsize) if args.figsize else None,
                dpi=args.dpi)
    except ValueError as err:
        print(f"error: {err}", file=sys.stderr)
        return 2

    print(f"[ok] {len(spectra)} spectra -> {out}")
    if peaks:
        print(peak_table(peaks))
    else:
        print("no emission maximum found (check --search-offset and the emission range)")
    if fit:
        r2 = "" if fit["r2"] is None else f", R2 {fit['r2']:.3f}"
        print(f"emission maximum vs excitation: slope {fit['slope']:.3f}{r2}; "
              f"{fit['shift_nm']:+.1f} nm over {fit['ex_range_nm']:g} nm of excitation")
    elif len(peaks) < 3:
        print("fewer than 3 maxima: no shift fit")
    for note in scatter_warnings(peaks):
        print(f"[warn] {note}", file=sys.stderr)
    if args.csv:
        write_csv(Path(args.csv), peaks)
        print(f"Wrote {len(peaks)} maxima to {args.csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
