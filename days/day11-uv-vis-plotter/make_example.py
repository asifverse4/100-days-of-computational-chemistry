"""Generate the SYNTHETIC demo spectra in example/synthetic_cqd.csv (fixed seed, reproducible).

The bands (near 240 nm and 340 nm, on a decaying baseline) only imitate the general shape often
reported for carbon dots. They are not measured data.
"""
import csv
from pathlib import Path

import numpy as np


def band(x, center, sigma, height):
    return height * np.exp(-0.5 * ((x - center) / sigma) ** 2)


def main() -> None:
    rng = np.random.default_rng(11)
    wl = np.arange(220.0, 601.0, 1.0)
    tail = 0.25 * np.exp(-(wl - 220.0) / 120.0)  # decaying baseline from scattering
    a = band(wl, 240, 18, 1.00) + band(wl, 340, 30, 0.35) + tail + rng.normal(0, 0.004, wl.size)
    b = band(wl, 244, 18, 0.80) + band(wl, 348, 32, 0.45) + tail + rng.normal(0, 0.004, wl.size)
    out = Path(__file__).parent / "example" / "synthetic_cqd.csv"
    with out.open("w", newline="") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(["wavelength_nm", "sample_A", "sample_B"])
        for row in zip(wl, a, b):
            writer.writerow([f"{row[0]:.0f}", f"{row[1]:.4f}", f"{row[2]:.4f}"])
    print(f"Wrote {wl.size} rows to {out}")


if __name__ == "__main__":
    main()
