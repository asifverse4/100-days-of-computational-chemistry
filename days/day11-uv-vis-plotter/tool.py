"""Day 11: UV-vis plotter with peak detection.

Plot one or more UV-vis spectra (CSV/TSV/TXT, wavelength in nm) and find the absorption
peaks: position, height, width and prominence. Handles instrument exports with header lines.

Usage:
    python tool.py example/synthetic_cqd.csv -o uvvis.png
    python tool.py a.csv b.csv --normalize --stack 0.3 -o compare.png
    python tool.py sample.csv --smooth 11 --prominence 0.02 --peaks-csv peaks.csv
    python tool.py sample.csv --zero-at 800 --xmin 250 --xmax 600 -o plot.svg
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

# Okabe-Ito palette: distinguishable with the common types of color blindness
PALETTE = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#000000"]
PEAK_COLUMNS = ["sample", "wavelength_nm", "absorbance", "fwhm_nm", "prominence"]
MIN_POINTS = 5


@dataclass
class Spectrum:
    name: str
    x: np.ndarray  # wavelength, nm, ascending
    y: np.ndarray  # absorbance (or %T before conversion)


def _delimiter(line: str) -> str | None:
    for d in ("\t", ";", ","):
        if d in line:
            return d
    return None  # whitespace


def _numbers(tokens: list[str], decimal_comma: bool) -> list[float] | None:
    try:
        return [float(t.replace(",", ".") if decimal_comma else t) for t in tokens]
    except ValueError:
        return None


def read_spectra(path: Path) -> list[Spectrum]:
    """Read a text file with the wavelength in column 1 and one or more intensity columns."""
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
            header = tokens  # keep the last text line before the data starts
    if not rows:
        raise ValueError(f"{path}: no numeric data found (expected wavelength and absorbance columns)")
    data = np.array(rows)
    order = np.argsort(data[:, 0], kind="stable")
    data = data[order]
    _, keep = np.unique(data[:, 0], return_index=True)
    data = data[np.sort(keep)]
    if data.shape[0] < MIN_POINTS:
        raise ValueError(f"{path}: only {data.shape[0]} data points, need at least {MIN_POINTS}")
    names = header[1:] if header and len(header) == ncols else []
    spectra = []
    for i in range(1, ncols):
        name = names[i - 1] if i - 1 < len(names) else (path.stem if ncols == 2 else f"{path.stem}_{i}")
        spectra.append(Spectrum(name, data[:, 0].copy(), data[:, i].copy()))
    return spectra


def percent_t_to_absorbance(t: np.ndarray) -> np.ndarray:
    if np.any(t <= 0):
        raise ValueError("--percent-t needs transmittance values above 0")
    return -np.log10(t / 100.0)


def crop(spec: Spectrum, xmin: float | None, xmax: float | None) -> Spectrum:
    mask = np.ones_like(spec.x, dtype=bool)
    if xmin is not None:
        mask &= spec.x >= xmin
    if xmax is not None:
        mask &= spec.x <= xmax
    if mask.sum() < MIN_POINTS:
        raise ValueError(f"{spec.name}: fewer than {MIN_POINTS} points between {xmin} and {xmax} nm")
    return Spectrum(spec.name, spec.x[mask], spec.y[mask])


def zero_at(spec: Spectrum, wavelength: float) -> Spectrum:
    """Subtract the absorbance at the wavelength closest to `wavelength` (scattering baseline)."""
    idx = int(np.argmin(np.abs(spec.x - wavelength)))
    return Spectrum(spec.name, spec.x, spec.y - spec.y[idx])


def smooth(spec: Spectrum, window: int) -> Spectrum:
    if window % 2 == 0 or window < MIN_POINTS:
        raise ValueError("--smooth needs an odd window of at least 5 points")
    if window > spec.y.size:
        raise ValueError(f"{spec.name}: smoothing window {window} is longer than the data ({spec.y.size} points)")
    return Spectrum(spec.name, spec.x, savgol_filter(spec.y, window, 3))


def normalize(spec: Spectrum) -> Spectrum:
    top = spec.y.max()
    if top <= 0:
        raise ValueError(f"{spec.name}: cannot normalize, the maximum is not positive")
    return Spectrum(spec.name, spec.x, spec.y / top)


def detect_peaks(spec: Spectrum, prominence: float | None = None, min_separation_nm: float = 10.0) -> list[dict]:
    """Find peaks. FWHM is the width at half of each peak's prominence, converted to nm."""
    y, x = spec.y, spec.x
    prom = prominence if prominence is not None else 0.05 * float(y.max() - y.min())
    step = float(np.median(np.diff(x)))
    distance = max(1, round(min_separation_nm / step))
    idx, props = find_peaks(y, prominence=max(prom, 1e-12), distance=distance)
    if idx.size == 0:
        return []
    _, _, left, right = peak_widths(y, idx, rel_height=0.5, prominence_data=(props["prominences"], props["left_bases"], props["right_bases"]))
    grid = np.arange(x.size)
    fwhm = np.interp(right, grid, x) - np.interp(left, grid, x)
    return [{"sample": spec.name, "wavelength_nm": float(x[i]), "absorbance": float(y[i]),
             "fwhm_nm": float(w), "prominence": float(p)}
            for i, w, p in zip(idx, fwhm, props["prominences"])]


def plot(spectra: list[Spectrum], peaks: list[dict], path: Path, *, stack: float = 0.0,
         figsize=(6.0, 4.0), dpi: int = 300, title: str | None = None, ylabel: str = "Absorbance (a.u.)",
         annotate: bool = True) -> None:
    fig, ax = plt.subplots(figsize=figsize)
    offsets = {s.name: i * stack for i, s in enumerate(spectra)}
    for i, s in enumerate(spectra):
        color = PALETTE[i % len(PALETTE)]
        ax.plot(s.x, s.y + offsets[s.name], color=color, lw=1.4, label=s.name)
        mine = [p for p in peaks if p["sample"] == s.name]
        if mine:
            ax.plot([p["wavelength_nm"] for p in mine], [p["absorbance"] + offsets[s.name] for p in mine],
                    "v", color=color, ms=5)
            if annotate:
                top = max(p["absorbance"] for p in mine) + offsets[s.name]
                for p in mine:
                    ax.annotate(f"{p['wavelength_nm']:.0f}", (p["wavelength_nm"], p["absorbance"] + offsets[s.name]),
                                textcoords="offset points", xytext=(0, 7), ha="center", fontsize=8, color=color)
                ax.set_ylim(top=max(ax.get_ylim()[1], top * 1.08))
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
    if len(spectra) > 1:
        ax.legend(frameon=False, fontsize=8)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi)
    plt.close(fig)


def write_peaks_csv(path: Path, peaks: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=PEAK_COLUMNS, lineterminator="\n")
        writer.writeheader()
        for p in peaks:
            writer.writerow({"sample": p["sample"], "wavelength_nm": round(p["wavelength_nm"], 2),
                             "absorbance": round(p["absorbance"], 4), "fwhm_nm": round(p["fwhm_nm"], 2),
                             "prominence": round(p["prominence"], 4)})


def peak_table(peaks: list[dict]) -> str:
    width = max([len("sample")] + [len(p["sample"]) for p in peaks])
    lines = [f"{'sample':<{width}}  {'peak (nm)':>9}  {'abs':>7}  {'FWHM (nm)':>9}  {'prominence':>10}"]
    for p in peaks:
        lines.append(f"{p['sample']:<{width}}  {p['wavelength_nm']:>9.1f}  {p['absorbance']:>7.3f}  "
                     f"{p['fwhm_nm']:>9.1f}  {p['prominence']:>10.3f}")
    return "\n".join(lines)


def process(spec: Spectrum, args) -> Spectrum:
    """Apply the processing chain: %T conversion, crop, baseline, smoothing, normalization."""
    if args.percent_t:
        spec = Spectrum(spec.name, spec.x, percent_t_to_absorbance(spec.y))
    spec = crop(spec, args.xmin, args.xmax)
    if args.zero_at is not None:
        spec = zero_at(spec, args.zero_at)
    if args.smooth:
        spec = smooth(spec, args.smooth)
    if args.normalize:
        spec = normalize(spec)
    return spec


def unique_names(groups: list[tuple[Path, list[Spectrum]]]) -> list[Spectrum]:
    """Flatten spectra from several files; prefix with the file name when files share column names."""
    flat = [(path, s) for path, spectra in groups for s in spectra]
    names = [s.name for _, s in flat]
    for path, s in flat:
        if len(groups) > 1 and (names.count(s.name) > 1 or s.name == path.stem):
            s.name = s.name if s.name == path.stem else f"{path.stem} {s.name}"
    return [s for _, s in flat]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+", help="spectrum files: wavelength (nm) in column 1, absorbance in the others")
    ap.add_argument("-o", "--output", default="uvvis.png", help="plot file: .png, .svg or .pdf (default uvvis.png)")
    ap.add_argument("--peaks-csv", help="also write the peak table to this CSV")
    ap.add_argument("--prominence", type=float, help="minimum peak prominence in plotted units (default: 5%% of the range)")
    ap.add_argument("--min-separation", type=float, default=10.0, help="minimum distance between peaks in nm (default 10)")
    ap.add_argument("--smooth", type=int, default=0, metavar="N", help="Savitzky-Golay smoothing window (odd, >= 5 points)")
    ap.add_argument("--zero-at", type=float, metavar="NM", help="subtract the absorbance at this wavelength (e.g. 800)")
    ap.add_argument("--normalize", action="store_true", help="scale each spectrum to a maximum of 1")
    ap.add_argument("--percent-t", action="store_true", help="input is %% transmittance, convert to absorbance")
    ap.add_argument("--xmin", type=float, help="shortest wavelength to use (nm)")
    ap.add_argument("--xmax", type=float, help="longest wavelength to use (nm)")
    ap.add_argument("--stack", type=float, default=0.0, help="vertical offset between spectra")
    ap.add_argument("--no-labels", action="store_true", help="do not write peak wavelengths on the plot")
    ap.add_argument("--title", help="plot title")
    ap.add_argument("--figsize", type=float, nargs=2, default=(6.0, 4.0), metavar=("W", "H"), help="figure size in inches")
    ap.add_argument("--dpi", type=int, default=300, help="resolution of PNG output (default 300)")
    args = ap.parse_args()

    out = Path(args.output)
    if out.suffix.lower() not in (".png", ".svg", ".pdf"):
        print("error: output must end in .png, .svg or .pdf", file=sys.stderr)
        return 2
    try:
        groups = [(Path(p), read_spectra(Path(p))) for p in args.inputs]
        spectra = [process(s, args) for s in unique_names(groups)]
        peaks = [p for s in spectra for p in detect_peaks(s, args.prominence, args.min_separation)]
        ylabel = "Normalized absorbance" if args.normalize else "Absorbance (a.u.)"
        plot(spectra, peaks, out, stack=args.stack, figsize=tuple(args.figsize), dpi=args.dpi,
             title=args.title, ylabel=ylabel, annotate=not args.no_labels)
    except ValueError as err:
        print(f"error: {err}", file=sys.stderr)
        return 2

    print(f"[ok] {len(spectra)} spectra -> {out}")
    if peaks:
        print(peak_table(peaks))
    else:
        print("no peaks found (try a smaller --prominence)")
    if args.peaks_csv:
        write_peaks_csv(Path(args.peaks_csv), peaks)
        print(f"Wrote {len(peaks)} peaks to {args.peaks_csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
