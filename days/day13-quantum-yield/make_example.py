"""Generate the SYNTHETIC demo data in example/ (fixed seed, reproducible).

A fluorescein-like reference in water (QY 0.925, n = 1.333) and a sample in ethanol (n = 1.361) whose
true quantum yield is set to 0.35. Intensities follow F = K * QY * f / n^2 (f = 1 - 10^-A) with 1 % noise.
The numbers are invented for the demo. They are not measured data.
"""
import csv
from pathlib import Path

import numpy as np

K = 2.0e7  # arbitrary instrument constant
TRUE_SAMPLE_QY = 0.35


def intensity(absorbance, qy, n, rng):
    f = 1 - 10 ** (-np.asarray(absorbance))
    return K * qy * f / n**2 * (1 + rng.normal(0, 0.01, f.size))


def write(path, absorbance, counts):
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(["absorbance", "integrated_intensity"])
        for a, i in zip(absorbance, counts):
            writer.writerow([f"{a:.3f}", f"{i:.0f}"])


def main() -> None:
    rng = np.random.default_rng(13)
    out = Path(__file__).parent / "example"
    ref_a = [0.012, 0.025, 0.041, 0.058, 0.077, 0.095]
    sample_a = [0.010, 0.022, 0.039, 0.055, 0.071, 0.090]
    write(out / "synthetic_reference.csv", ref_a, intensity(ref_a, 0.925, 1.333, rng))
    write(out / "synthetic_sample.csv", sample_a, intensity(sample_a, TRUE_SAMPLE_QY, 1.361, rng))
    print(f"Wrote demo files to {out}")


if __name__ == "__main__":
    main()
