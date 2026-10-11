"""Day 16: FTIR peak picker with functional-group suggestions.

Find the peaks in an FTIR spectrum (transmittance or absorbance), and list the functional groups whose
characteristic band ranges they fall in. The suggestions are a starting point, not an assignment: many groups
overlap, and positions shift with hydrogen bonding, conjugation and ring strain.

Usage:
    python tool.py example/synthetic_ftir.csv -o ftir.png
    python tool.py spectrum.csv --prominence 0.02 --smooth 11 --csv peaks.csv
    python tool.py spectrum.csv --signal absorbance --xmin 600 --xmax 3800 --top 2
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
from scipy.signal import find_peaks, peak_widths, savgol_filter

MIN_POINTS = 10
FINGERPRINT_BELOW = 1500.0
PEAK_COLUMNS = ["wavenumber_cm1", "height_absorbance", "relative_intensity", "fwhm_cm1", "suggestions"]
INTENSITY_WORDS = {"s": {"strong"}, "m": {"medium"}, "w": {"weak"}, "m-s": {"medium", "strong"},
                   "w-m": {"weak", "medium"}, "v": {"weak", "medium", "strong"}}

SOURCES = {
    "UCLA": "UCLA WebSpectra IR table (https://webspectra.chem.ucla.edu/irtable.html)",
    "PMC": "Acta Biomater. 2020, 105, 280, Table 2 (PMC7134197)",
    "text": "standard textbook range; not checked against a source while writing this tool",
}


@dataclass(frozen=True)
class Band:
    group: str
    lo: float
    hi: float
    intensity: str  # s, m, w, m-s, w-m or v (variable)
    shape: str  # "broad", "sharp" or ""
    source: str


# Ranges in cm-1. "UCLA" entries use the ranges of that table as printed (single "~" values are widened by about 40).
TABLE = [
    Band("Alcohol/phenol O-H stretch (hydrogen-bonded)", 3200, 3550, "s", "broad", "UCLA"),
    Band("Amide N-H stretch", 3500, 3700, "m", "", "UCLA"),
    Band("Amine N-H stretch", 3300, 3500, "m", "", "UCLA"),
    Band("Alkynyl C-H stretch (about 3300)", 3260, 3340, "s", "sharp", "UCLA"),
    Band("Alkenyl C-H stretch", 3010, 3100, "m", "", "UCLA"),
    Band("Aromatic C-H stretch (about 3030)", 3000, 3100, "v", "", "UCLA"),
    Band("Alkyl C-H stretch", 2850, 2950, "m-s", "", "UCLA"),
    Band("Carboxylic acid O-H stretch", 2500, 3000, "v", "broad", "UCLA"),
    Band("Nitrile C\u2261N stretch", 2220, 2260, "m", "sharp", "UCLA"),
    Band("Alkynyl C\u2261C stretch", 2100, 2260, "v", "", "UCLA"),
    Band("Carboxylic acid C=O stretch", 1710, 1780, "s", "", "UCLA"),
    Band("Ester C=O stretch", 1735, 1750, "s", "", "UCLA"),
    Band("Ketone C=O stretch", 1680, 1750, "s", "", "UCLA"),
    Band("Aldehyde C=O stretch", 1690, 1740, "s", "", "UCLA"),
    Band("Amide C=O stretch (amide I)", 1630, 1690, "s", "", "UCLA"),
    Band("Alkenyl C=C stretch", 1620, 1680, "v", "", "UCLA"),
    Band("Aromatic C=C stretch", 1500, 1700, "m", "", "UCLA"),
    Band("Aromatic C-H bend (out of plane)", 680, 860, "s", "", "UCLA"),
    Band("Amine N-H bend", 1550, 1650, "m", "", "PMC"),
    Band("Free O-H stretch", 3584, 3650, "m", "sharp", "text"),
    Band("Aldehyde C-H stretch", 2700, 2850, "w", "", "text"),
    Band("Thiol S-H stretch", 2550, 2600, "w", "", "text"),
    Band("Imine C=N stretch", 1640, 1690, "m", "", "text"),
    Band("Amide N-H bend (amide II)", 1515, 1570, "m-s", "", "text"),
    Band("Nitro N=O stretch (asymmetric)", 1475, 1550, "s", "", "text"),
    Band("Nitro N=O stretch (symmetric)", 1290, 1360, "s", "", "text"),
    Band("Carboxylate COO- stretch (asymmetric)", 1540, 1650, "s", "", "text"),
    Band("Carboxylate COO- stretch (symmetric)", 1360, 1450, "s", "", "text"),
    Band("Carboxylic acid O-H bend", 1395, 1440, "m", "", "text"),
    Band("CH2/CH3 C-H bend", 1430, 1480, "m", "", "text"),
    Band("CH3 C-H bend (symmetric)", 1365, 1385, "m", "", "text"),
    Band("C-N stretch (amines)", 1020, 1250, "m", "", "text"),
    Band("C-O stretch (alcohols, ethers, esters, acids)", 1000, 1300, "s", "", "text"),
    Band("C-F stretch", 1000, 1400, "s", "", "text"),
    Band("C-Cl stretch", 550, 850, "s", "", "text"),
]

REGIONS = [(3700, 3200, "O-H, N-H"), (3100, 2800, "C-H"), (2300, 2100, "C\u2261C, C\u2261N"),
           (1850, 1650, "C=O"), (1650, 1450, "C=C, N-H"), (1500, 400, "fingerprint")]


@dataclass
class Spectrum:
    name: str
    x: np.ndarray  # wavenumber, cm-1, ascending
    y: np.ndarray  # raw signal as read


# ---------- reading and signal handling ----------
def read_spectrum(path: Path, column: str | None = None) -> Spectrum:
    """Read wavenumber (column 1) and one signal column; metadata lines and headers are skipped."""
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
        raise ValueError(f"{path}: no numeric data found (expected wavenumber and signal columns)")
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


def detect_signal(y: np.ndarray, mode: str = "auto") -> str:
    """Return 'absorbance', 'transmittance' (percent) or 'transmittance-fraction'."""
    if mode != "auto":
        return mode
    return "transmittance" if float(np.nanmax(y)) > 5 else "absorbance"


def to_absorbance(y: np.ndarray, signal: str) -> np.ndarray:
    if signal == "absorbance":
        return y.copy()
    scale = 100.0 if signal == "transmittance" else 1.0
    if np.any(y <= 0):
        raise ValueError("transmittance values must be above 0 (a spectrum with zero or negative values may be "
                         "absorbance: use --signal absorbance)")
    return -np.log10(y / scale)


# ---------- peaks ----------
def refine_position(x: np.ndarray, y: np.ndarray, i: int, top_fraction: float = 0.9) -> float:
    """Peak position from a parabola fitted to the top of the peak (at least 5 points)."""
    threshold = y[i] - (1 - top_fraction) * abs(y[i])
    lo = hi = i
    while lo > 0 and y[lo - 1] >= threshold:
        lo -= 1
    while hi < y.size - 1 and y[hi + 1] >= threshold:
        hi += 1
    lo, hi = min(lo, max(i - 2, 0)), max(hi, min(i + 2, y.size - 1))
    a, b, _ = np.polyfit(x[lo:hi + 1] - x[i], y[lo:hi + 1], 2)
    if a < 0 and x[lo] - x[i] <= -b / (2 * a) <= x[hi] - x[i]:
        return float(x[i] - b / (2 * a))
    return float(x[i])


def classify(relative: float) -> str:
    return "strong" if relative >= 0.6 else "medium" if relative >= 0.25 else "weak"


def find_peaks_in(x: np.ndarray, a: np.ndarray, prominence: float | None = None, min_separation: float = 10.0) -> list[dict]:
    """Peaks of an absorbance-like curve. Returns dicts sorted from high to low wavenumber."""
    prom = prominence if prominence is not None else 0.05 * float(a.max() - a.min())
    step = float(np.median(np.diff(x)))
    idx, props = find_peaks(a, prominence=max(prom, 1e-12), distance=max(1, round(min_separation / step)))
    if idx.size == 0:
        return []
    _, _, left, right = peak_widths(a, idx, rel_height=0.5, prominence_data=(
        props["prominences"], props["left_bases"], props["right_bases"]))
    grid = np.arange(x.size)
    fwhm = np.interp(right, grid, x) - np.interp(left, grid, x)
    top = float(props["prominences"].max())
    peaks = [{"wavenumber_cm1": refine_position(x, a, int(i)), "height": float(a[i]), "prominence": float(p),
              "relative": float(p) / top, "intensity_class": classify(float(p) / top), "fwhm_cm1": float(w)}
             for i, p, w in zip(idx, props["prominences"], fwhm)]
    return sorted(peaks, key=lambda p: -p["wavenumber_cm1"])


# ---------- functional-group suggestions ----------
def range_text(b: Band) -> str:
    shape = f", {b.shape}" if b.shape else ""
    return f"{b.lo:g}-{b.hi:g}, {b.intensity}{shape}"


def suggest(position: float, intensity_class: str, fwhm: float | None, top: int = 3, tolerance: float = 10.0) -> list[dict]:
    """Rank the table entries that contain the position (or lie within `tolerance` cm-1 of it)."""
    found = []
    for b in TABLE:
        if not (b.lo - tolerance <= position <= b.hi + tolerance):
            continue
        outside = max(b.lo - position, position - b.hi, 0.0)
        score = 1.0 - 0.6 * min(outside / tolerance, 1.0) if tolerance > 0 else 1.0  # fades smoothly past the edge
        # Position dominates. Intensity is only relative to the strongest peak here, so it counts for little.
        score += (0.1 if b.intensity == "v" else 0.2) if intensity_class in INTENSITY_WORDS[b.intensity] else -0.1
        if fwhm is not None and b.shape == "broad":
            score += 0.3 if fwhm >= 100 else -0.3
        if fwhm is not None and b.shape == "sharp":
            score += 0.2 if fwhm < 60 else -0.2
        score += 0.3 * (1 - min((b.hi - b.lo) / 1500.0, 1.0))  # narrower ranges are more specific
        found.append({"group": b.group, "range": range_text(b), "source": b.source, "score": round(score, 3)})
    return sorted(found, key=lambda s: -s["score"])[:top]


# ---------- output ----------
def plot_ftir(spec: Spectrum, signal: str, peaks: list[dict], path: Path, *, regions: bool, max_labels: int,
              title: str | None, figsize: tuple[float, float], dpi: int) -> None:
    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(spec.x, spec.y, color="#0072B2", lw=1.1)
    transmittance = signal != "absorbance"
    if regions:
        for k, (hi, lo, label) in enumerate(REGIONS):
            ax.axvspan(lo, hi, color="gray", alpha=0.07 if k % 2 == 0 else 0.0, lw=0)
            ax.text((hi + lo) / 2, 1.01, label, transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=7)
    scale = 100.0 if signal == "transmittance" else 1.0
    shown = sorted(peaks, key=lambda p: -p["prominence"])[:max_labels]
    for p in shown:
        pos = p["wavenumber_cm1"]
        y = float(np.interp(pos, spec.x, spec.y))
        ax.plot(pos, y, "v" if not transmittance else "^", color="#D55E00", ms=4)
        ax.annotate(f"{pos:.0f}", (pos, y), textcoords="offset points", xytext=(0, -9 if transmittance else 6),
                    ha="center", va="top" if transmittance else "bottom", rotation=90, fontsize=7, color="#D55E00")
    ax.set_xlim(spec.x.max(), spec.x.min())
    ax.set_xlabel("Wavenumber (cm$^{-1}$)")
    ax.set_ylabel("Transmittance (%)" if signal == "transmittance" else
                  "Transmittance" if transmittance else "Absorbance")
    lo_y, hi_y = float(spec.y.min()), float(spec.y.max())
    span = hi_y - lo_y or scale
    ax.set_ylim(lo_y - (0.22 * span if transmittance else 0.0), hi_y + 0.1 * span + (0.0 if transmittance else 0.12 * span))
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    if title:
        ax.set_title(title, pad=16 if regions else 6)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi)
    plt.close(fig)


def describe(peaks: list[dict], top: int, tolerance: float) -> list[str]:
    lines = []
    for p in peaks:
        shape = "broad" if p["fwhm_cm1"] >= 100 else "sharp" if p["fwhm_cm1"] < 60 else "medium width"
        lines.append(f"{p['wavenumber_cm1']:7.0f} cm-1  {p['intensity_class']}, {shape} (FWHM {p['fwhm_cm1']:.0f}), "
                     f"A = {p['height']:.3f}")
        matches = suggest(p["wavenumber_cm1"], p["intensity_class"], p["fwhm_cm1"], top, tolerance)
        if p["wavenumber_cm1"] < FINGERPRINT_BELOW:
            lines.append("            (fingerprint region: assignments are tentative)")
        if not matches:
            lines.append("            no band in the table")
        for k, m in enumerate(matches, 1):
            lines.append(f"            {k}. {m['group']} ({m['range']})")
    return lines


def write_csv(path: Path, peaks: list[dict], top: int, tolerance: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=PEAK_COLUMNS, lineterminator="\n")
        writer.writeheader()
        for p in peaks:
            names = " | ".join(m["group"] for m in suggest(p["wavenumber_cm1"], p["intensity_class"], p["fwhm_cm1"],
                                                            top, tolerance))
            writer.writerow({"wavenumber_cm1": round(p["wavenumber_cm1"], 1), "height_absorbance": round(p["height"], 4),
                             "relative_intensity": p["intensity_class"], "fwhm_cm1": round(p["fwhm_cm1"], 1),
                             "suggestions": names})


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", help="spectrum file: wavenumber (cm-1) in column 1, signal in column 2")
    ap.add_argument("-o", "--output", help="save the annotated spectrum (.png, .svg or .pdf)")
    ap.add_argument("--csv", help="save the peak table to a CSV")
    ap.add_argument("--column", help="signal column to use (header name or number; default the first)")
    ap.add_argument("--signal", choices=["auto", "absorbance", "transmittance", "transmittance-fraction"], default="auto",
                    help="auto treats values above 5 as percent transmittance (default)")
    ap.add_argument("--prominence", type=float, help="minimum peak prominence in absorbance units (default 5%% of the range)")
    ap.add_argument("--min-separation", type=float, default=10.0, help="minimum distance between peaks in cm-1 (default 10)")
    ap.add_argument("--smooth", type=int, default=0, metavar="N", help="Savitzky-Golay window (odd, at least 5 points)")
    ap.add_argument("--xmin", type=float, help="lowest wavenumber to use (cm-1)")
    ap.add_argument("--xmax", type=float, help="highest wavenumber to use (cm-1)")
    ap.add_argument("--top", type=int, default=3, help="suggestions per peak (default 3)")
    ap.add_argument("--tolerance", type=float, default=10.0, help="accept peaks this far outside a band range (cm-1, default 10)")
    ap.add_argument("--max-labels", type=int, default=25, help="most peaks to label on the plot (default 25)")
    ap.add_argument("--no-regions", action="store_true", help="do not draw the region labels on the plot")
    ap.add_argument("--title", help="plot title")
    ap.add_argument("--figsize", type=float, nargs=2, default=(8.0, 4.0), metavar=("W", "H"), help="figure size in inches")
    ap.add_argument("--dpi", type=int, default=300, help="PNG resolution (default 300)")
    args = ap.parse_args()

    try:
        if args.output and Path(args.output).suffix.lower() not in (".png", ".svg", ".pdf"):
            raise ValueError("output must end in .png, .svg or .pdf")
        if args.top < 1 or args.tolerance < 0:
            raise ValueError("--top must be at least 1 and --tolerance must not be negative")
        spec = read_spectrum(Path(args.input), args.column)
        mask = np.ones_like(spec.x, dtype=bool)
        if args.xmin is not None:
            mask &= spec.x >= args.xmin
        if args.xmax is not None:
            mask &= spec.x <= args.xmax
        if mask.sum() < MIN_POINTS:
            raise ValueError(f"fewer than {MIN_POINTS} points in the selected range")
        spec = Spectrum(spec.name, spec.x[mask], spec.y[mask])
        signal = detect_signal(spec.y, args.signal)
        absorbance = to_absorbance(spec.y, signal)
        if args.smooth:
            if args.smooth % 2 == 0 or args.smooth < 5 or args.smooth > absorbance.size:
                raise ValueError("--smooth needs an odd window of at least 5 points, shorter than the data")
            absorbance = savgol_filter(absorbance, args.smooth, 3)
        peaks = find_peaks_in(spec.x, absorbance, args.prominence, args.min_separation)
        if args.output:
            plot_ftir(spec, signal, peaks, Path(args.output), regions=not args.no_regions, max_labels=args.max_labels,
                      title=args.title, figsize=tuple(args.figsize), dpi=args.dpi)
    except ValueError as err:
        print(f"error: {err}", file=sys.stderr)
        return 2

    shown = {"absorbance": "absorbance", "transmittance": "percent transmittance",
             "transmittance-fraction": "fractional transmittance"}[signal]
    print(f"{len(peaks)} peaks in {spec.x.min():.0f}-{spec.x.max():.0f} cm-1 (signal read as {shown}"
          + ("; treated as absorbance, use --signal transmittance-fraction for 0-1 transmittance"
             if args.signal == "auto" and signal == "absorbance" and float(np.nanmax(spec.y)) <= 1.2 else "") + ")")
    if signal == "absorbance" and float(np.nanmax(spec.y)) > 5:
        print("[warn] values above 5 look like percent transmittance, not absorbance: check --signal", file=sys.stderr)
    print("\n".join(describe(peaks, args.top, args.tolerance)))
    if args.output:
        print(f"[ok] plot -> {args.output}")
    if args.csv:
        write_csv(Path(args.csv), peaks, args.top, args.tolerance)
        print(f"Wrote {len(peaks)} peaks to {args.csv}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except BrokenPipeError:  # e.g. piping into `head`
        sys.exit(0)
