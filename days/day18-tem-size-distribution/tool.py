"""Day 18: TEM particle size distribution (histogram plus lognormal fit).

Turn a list of particle sizes measured on TEM images into a histogram with a fitted lognormal curve and
the numbers a paper needs: mean and SD, median, geometric mean and geometric SD, with confidence
intervals, plus checks of how well a lognormal (and a normal) distribution describes the data.

Usage:
    python tool.py example/synthetic_sizes.csv --column diameter_nm -o sizes.png
    python tool.py example/synthetic_imagej.csv --column Area --from-area --unit-scale 0.25 --csv summary.csv
    python tool.py sizes.csv --bin-width 0.5 --density --also-normal --hist-csv bins.csv
"""
import argparse
import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

MIN_PARTICLES = 5
BIN_METHODS = {"auto", "fd", "sturges", "scott", "doane", "rice", "sqrt"}
SUMMARY_COLUMNS = [
    "n", "mean", "sd", "sem", "mean_ci95_low", "mean_ci95_high", "cv_percent", "pdi", "median", "min", "max",
    "geometric_mean", "geometric_mean_ci95_low", "geometric_mean_ci95_high", "geometric_sd", "mu_ln", "sigma_ln",
    "mu_ln_se", "sigma_ln_se", "lognormal_mean", "lognormal_mode", "shapiro_p_normal", "shapiro_p_lognormal",
    "delta_aic",
]
COLORS = {"hist": "#9ecae1", "edge": "#3182bd", "lognormal": "#D55E00", "normal": "#009E73"}


# ---------- reading ----------
def _split(line: str) -> tuple[list[str], str | None]:
    for delim in ("\t", ";", ","):
        if delim in line:
            return [t.strip() for t in line.split(delim)], delim
    return line.split(), None


def read_columns(path: Path) -> tuple[list[str], np.ndarray]:
    """Read a text table of numbers. Returns (column names, array of shape rows x columns)."""
    try:
        lines = [ln for ln in path.read_text(errors="replace").splitlines() if ln.strip()]
    except OSError as err:
        raise ValueError(f"cannot read {path}: {err.strerror}") from None
    rows, header, ncols = [], None, 0
    for line in lines:
        tokens, delim = _split(line)
        try:
            values = [float(t.replace(",", ".") if delim == ";" else t) if t != "" else np.nan for t in tokens]
        except ValueError:
            values = None
        if values is not None and not all(np.isnan(v) for v in values):
            if not rows:
                ncols = len(values)
            if len(values) == ncols:
                rows.append(values)
        elif not rows:
            header = tokens
    if not rows:
        raise ValueError(f"{path}: no numbers found")
    names = header if header and len(header) == ncols else []
    names = [names[i] if i < len(names) and names[i] else f"column {i + 1}" for i in range(ncols)]
    return names, np.array(rows)


def pick_column(names: list[str], data: np.ndarray, column: str | None, path: Path) -> tuple[str, np.ndarray]:
    """The column to use: the only one, or the one named (or numbered from 1) with --column."""
    if data.shape[1] == 1:
        return names[0], data[:, 0]
    if column is None:
        raise ValueError(f"{path}: the file has {data.shape[1]} columns ({', '.join(names)}); choose one with --column")
    lower = [n.lower() for n in names]
    if column.lower() in lower:
        index = lower.index(column.lower())
    elif column.isdigit() and 1 <= int(column) <= len(names):
        index = int(column) - 1
    else:
        raise ValueError(f"{path}: column {column!r} not found. Available: {', '.join(names)}")
    return names[index], data[:, index]


def to_diameters(values: np.ndarray, from_area: bool = False, scale: float = 1.0) -> np.ndarray:
    """Particle diameters in length units. `scale` converts the input length unit (for example pixels) to the output."""
    if scale <= 0:
        raise ValueError("--unit-scale must be positive")
    values = values[np.isfinite(values)]
    if values.size < MIN_PARTICLES:
        raise ValueError(f"only {values.size} values, need at least {MIN_PARTICLES}")
    if np.any(values <= 0):
        raise ValueError(f"sizes must be above 0 ({int(np.sum(values <= 0))} values are 0 or negative)")
    return (2.0 * np.sqrt(values / np.pi) if from_area else values) * scale


# ---------- statistics ----------
def describe(d: np.ndarray) -> dict:
    """Number-weighted statistics, the lognormal fit with standard errors, and goodness-of-fit checks."""
    n = d.size
    if float(np.ptp(d)) == 0.0:
        raise ValueError("all the values are identical: there is no distribution to describe")
    t = stats.t.ppf(0.975, n - 1)
    mean, sd = float(d.mean()), float(d.std(ddof=1))
    sem = sd / np.sqrt(n)
    ln = np.log(d)
    mu, sigma = float(ln.mean()), float(ln.std(ddof=1))  # lognormal parameters from the sample statistics of ln d
    mu_se, sigma_se = sigma / np.sqrt(n), sigma / np.sqrt(2.0 * (n - 1))
    # log-likelihoods (maximum-likelihood widths) for a normal and a lognormal model, for AIC
    sd_mle, sigma_mle = float(d.std()), float(ln.std())
    loglik_normal = float(np.sum(stats.norm.logpdf(d, mean, sd_mle)))
    loglik_lognormal = float(np.sum(stats.norm.logpdf(ln, mu, sigma_mle) - ln))
    shapiro = (lambda x: float(stats.shapiro(x).pvalue) if 3 <= n <= 5000 else None)
    return {
        "n": n, "mean": mean, "sd": sd, "sem": sem, "mean_ci95_low": mean - t * sem, "mean_ci95_high": mean + t * sem,
        "cv_percent": 100.0 * sd / mean, "pdi": (sd / mean) ** 2, "median": float(np.median(d)),
        "min": float(d.min()), "max": float(d.max()),
        "geometric_mean": float(np.exp(mu)), "geometric_mean_ci95_low": float(np.exp(mu - t * mu_se)),
        "geometric_mean_ci95_high": float(np.exp(mu + t * mu_se)), "geometric_sd": float(np.exp(sigma)),
        "mu_ln": mu, "sigma_ln": sigma, "mu_ln_se": mu_se, "sigma_ln_se": sigma_se,
        "lognormal_mean": float(np.exp(mu + sigma**2 / 2)), "lognormal_mode": float(np.exp(mu - sigma**2)),
        "shapiro_p_normal": shapiro(d), "shapiro_p_lognormal": shapiro(ln),
        "delta_aic": float((2 * 2 - 2 * loglik_normal) - (2 * 2 - 2 * loglik_lognormal)),
    }


def bin_edges(d: np.ndarray, bins: str = "auto", width: float | None = None) -> np.ndarray:
    """Histogram edges: a fixed bin width, a number of bins, or a numpy rule (auto, fd, sturges, ...)."""
    if width is not None:
        if width <= 0:
            raise ValueError("--bin-width must be positive")
        start = np.floor(d.min() / width) * width
        return np.arange(start, d.max() + width, width)
    if bins.isdigit():
        if int(bins) < 1:
            raise ValueError("--bins must be at least 1")
        return np.histogram_bin_edges(d, bins=int(bins))
    if bins not in BIN_METHODS:
        raise ValueError(f"unknown --bins {bins!r}: use a number or one of {', '.join(sorted(BIN_METHODS))}")
    return np.histogram_bin_edges(d, bins=bins)


def warnings_for(d: np.ndarray, s: dict, edges: np.ndarray, resolution: float | None) -> list[str]:
    notes = []
    if s["n"] < 30:
        notes.append(f"only {s['n']} particles: the mean, and especially the width of the distribution, are very uncertain")
    elif s["n"] < 100:
        notes.append(f"{s['n']} particles is a small sample: the fitted width and the histogram shape are uncertain")
    if resolution and float(np.median(np.diff(edges))) < resolution:
        notes.append(f"the bins ({np.median(np.diff(edges)):.3g}) are narrower than the measurement resolution "
                     f"({resolution:g}): the histogram will show empty bins and false structure")
    if s["max"] > 5 * s["median"]:
        notes.append(f"the largest particle ({s['max']:.3g}) is more than 5 times the median: check for aggregates "
                     "or mis-measured particles")
    if s["delta_aic"] < 0:
        notes.append(f"a normal distribution fits better than a lognormal one (AIC difference {s['delta_aic']:.1f}): "
                     "report the plain mean and SD")
    if s["shapiro_p_lognormal"] is not None and s["shapiro_p_lognormal"] < 0.05:
        notes.append(f"the log of the sizes is not normal at the 5% level (Shapiro-Wilk p = {s['shapiro_p_lognormal']:.3g}): "
                     "a single lognormal may not describe the data (look for two populations or aggregates)")
    return notes


# ---------- output ----------
def lognormal_pdf(x: np.ndarray, mu: float, sigma: float) -> np.ndarray:
    return np.exp(-((np.log(x) - mu) ** 2) / (2 * sigma**2)) / (x * sigma * np.sqrt(2 * np.pi))


def plot_distribution(d: np.ndarray, s: dict, edges: np.ndarray, path: Path, *, density: bool, also_normal: bool,
                      qq: bool, unit: str, title: str | None, figsize: tuple[float, float] | None, dpi: int) -> None:
    fig, axes = plt.subplots(1, 2 if qq else 1, figsize=figsize or ((9.5, 4.0) if qq else (5.5, 4.0)),
                             squeeze=False, gridspec_kw={"width_ratios": [1.6, 1]} if qq else None)
    ax = axes[0][0]
    ax.hist(d, bins=edges, density=density, color=COLORS["hist"], edgecolor=COLORS["edge"], lw=0.8)
    factor = 1.0 if density else s["n"] * float(np.mean(np.diff(edges)))
    x = np.linspace(max(d.min() * 0.5, 1e-9), d.max() * 1.15, 400)
    ax.plot(x, factor * lognormal_pdf(x, s["mu_ln"], s["sigma_ln"]), color=COLORS["lognormal"], lw=1.8, label="lognormal fit")
    if also_normal:
        ax.plot(x, factor * stats.norm.pdf(x, s["mean"], s["sd"]), color=COLORS["normal"], lw=1.4, ls="--", label="normal fit")
    ax.axvline(s["median"], color="black", lw=1.0, ls=":")
    ax.set_xlabel(f"Particle diameter ({unit})")
    ax.set_ylabel("Probability density" if density else "Number of particles")
    ax.legend(frameon=False, fontsize=8, loc="center right")
    ax.text(0.97, 0.97, f"n = {s['n']}\nmean = {s['mean']:.2f} $\\pm$ {s['sd']:.2f} {unit}\n"
                        f"median = {s['median']:.2f} {unit}\ngeometric SD = {s['geometric_sd']:.2f}",
            transform=ax.transAxes, ha="right", va="top", fontsize=8,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.8})
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    if qq:
        ax2 = axes[0][1]
        (theory, ordered), (slope, intercept, _) = stats.probplot(np.log(d), dist="norm")
        ax2.plot(theory, ordered, "o", color=COLORS["edge"], ms=3)
        ax2.plot(theory, slope * theory + intercept, color=COLORS["lognormal"], lw=1.4)
        p = s["shapiro_p_lognormal"]
        ax2.set_title("Normal Q-Q plot of ln(d)" + ("" if p is None else f"  (Shapiro p = {p:.2g})"), fontsize=9)
        ax2.set_xlabel("Normal quantiles")
        ax2.set_ylabel("ln(d)")
        for side in ("top", "right"):
            ax2.spines[side].set_visible(False)
    if title:
        fig.suptitle(title)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi)
    plt.close(fig)


def fmt_p(p: float | None) -> str:
    return "n/a" if p is None else f"{p:.3g}"


def describe_text(s: dict, source: str, unit: str) -> str:
    lines = [
        f"n = {s['n']} particles ({source}, {unit})",
        "number-weighted statistics",
        f"  mean +/- SD        {s['mean']:.3g} +/- {s['sd']:.3g} {unit}   (CV {s['cv_percent']:.0f} %, PDI = (SD/mean)^2 = {s['pdi']:.3f})",
        f"  mean, 95% CI       {s['mean_ci95_low']:.3g} to {s['mean_ci95_high']:.3g} {unit}",
        f"  median             {s['median']:.3g} {unit}   (range {s['min']:.3g} to {s['max']:.3g})",
        "lognormal fit (mean and SD of ln d)",
        (f"  geometric mean     {s['geometric_mean']:.3g} {unit}   (95% CI {s['geometric_mean_ci95_low']:.3g} to "
         f"{s['geometric_mean_ci95_high']:.3g})"),
        f"  geometric SD       {s['geometric_sd']:.3g}   (ln sigma = {s['sigma_ln']:.3f} +/- {s['sigma_ln_se']:.3f})",
        f"  mean / mode        {s['lognormal_mean']:.3g} / {s['lognormal_mode']:.3g} {unit}",
        "checks",
        f"  Shapiro-Wilk p     {fmt_p(s['shapiro_p_normal'])} for d (normal), {fmt_p(s['shapiro_p_lognormal'])} for ln d (lognormal)",
        f"  AIC normal - lognormal = {s['delta_aic']:.1f}   (positive favours the lognormal)",
    ]
    return "\n".join(lines)


def write_summary_csv(path: Path, s: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=SUMMARY_COLUMNS, lineterminator="\n")
        writer.writeheader()
        writer.writerow({k: "" if s[k] is None else (s[k] if k == "n" else float(f"{s[k]:.6g}")) for k in SUMMARY_COLUMNS})


def write_hist_csv(path: Path, d: np.ndarray, edges: np.ndarray) -> None:
    counts, _ = np.histogram(d, bins=edges)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(["bin_low", "bin_high", "bin_center", "count", "fraction"])
        for lo, hi, c in zip(edges[:-1], edges[1:], counts):
            writer.writerow([f"{lo:.6g}", f"{hi:.6g}", f"{(lo + hi) / 2:.6g}", int(c), f"{c / d.size:.6g}"])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", help="text file with one measured value per row (a column of a table such as an ImageJ export)")
    ap.add_argument("--column", help="column to use (header name or number from 1); needed if the file has several")
    ap.add_argument("--from-area", action="store_true", help="the column holds particle areas: use the equivalent circle diameter")
    ap.add_argument("--unit-scale", type=float, default=1.0,
                    help="multiply lengths by this (for example nm per pixel); with --from-area it scales the diameter")
    ap.add_argument("--unit-label", default="nm", help="length unit shown in labels (default nm)")
    ap.add_argument("--bins", default="auto", help="number of bins or a rule: auto, fd, sturges, scott, doane, rice, sqrt")
    ap.add_argument("--bin-width", type=float, help="fixed bin width (overrides --bins)")
    ap.add_argument("--density", action="store_true", help="plot probability density instead of counts")
    ap.add_argument("--also-normal", action="store_true", help="also draw the normal fit")
    ap.add_argument("--no-qq", action="store_true", help="leave out the Q-Q plot of ln(d)")
    ap.add_argument("-o", "--output", help="save the plot (.png, .svg or .pdf)")
    ap.add_argument("--csv", help="save the summary statistics to a one-row CSV")
    ap.add_argument("--hist-csv", help="save the histogram bins and counts to a CSV")
    ap.add_argument("--title", help="plot title")
    ap.add_argument("--figsize", type=float, nargs=2, metavar=("W", "H"), help="figure size in inches")
    ap.add_argument("--dpi", type=int, default=300, help="PNG resolution (default 300)")
    args = ap.parse_args()

    try:
        if args.output and Path(args.output).suffix.lower() not in (".png", ".svg", ".pdf"):
            raise ValueError("output must end in .png, .svg or .pdf")
        path = Path(args.input)
        names, data = read_columns(path)
        column, values = pick_column(names, data, args.column, path)
        d = to_diameters(values, args.from_area, args.unit_scale)
        s = describe(d)
        edges = bin_edges(d, args.bins, args.bin_width)
        if args.output:
            plot_distribution(d, s, edges, Path(args.output), density=args.density, also_normal=args.also_normal,
                              qq=not args.no_qq, unit=args.unit_label, title=args.title,
                              figsize=tuple(args.figsize) if args.figsize else None, dpi=args.dpi)
    except ValueError as err:
        print(f"error: {err}", file=sys.stderr)
        return 2

    source = f"column {column!r}" + (" as equivalent circle diameter" if args.from_area else "")
    if args.unit_scale != 1.0:
        source += f", scaled by {args.unit_scale:g}"
    print(describe_text(s, source, args.unit_label))
    resolution = args.unit_scale if args.unit_scale != 1.0 else None
    for note in warnings_for(d, s, edges, resolution):
        print(f"[warn] {note}", file=sys.stderr)
    if args.output:
        print(f"[ok] plot -> {args.output}")
    if args.csv:
        write_summary_csv(Path(args.csv), s)
        print(f"Wrote the summary to {args.csv}")
    if args.hist_csv:
        write_hist_csv(Path(args.hist_csv), d, edges)
        print(f"Wrote {len(edges) - 1} bins to {args.hist_csv}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except BrokenPipeError:  # e.g. piping into `head`
        sys.exit(0)
