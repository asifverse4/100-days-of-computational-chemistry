"""Generate the SYNTHETIC demo files in example/ (fixed seed, reproducible).

Absorption: a strong band at 240 nm and a weaker, longest-wavelength band at 345 nm on a small tail.
Emission: one broad band at 450 nm. The invented numbers show why the choice of absorption band matters.
They are not measured data. samples.csv is a made-up table of band positions for the batch mode.
"""
import csv
from pathlib import Path

import numpy as np


def gaussian(x, center, sigma, height):
    return height * np.exp(-0.5 * ((x - center) / sigma) ** 2)


def write_spectrum(path, header, x, y):
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(["wavelength_nm", header])
        for a, b in zip(x, y):
            writer.writerow([f"{a:.0f}", f"{b:.4f}"])


def main() -> None:
    rng = np.random.default_rng(15)
    out = Path(__file__).parent / "example"
    x_abs = np.arange(220.0, 521.0, 1.0)
    absorbance = (gaussian(x_abs, 240, 18, 1.0) + gaussian(x_abs, 345, 26, 0.42)
                  + 0.05 * np.exp(-(x_abs - 220) / 150) + rng.normal(0, 0.003, x_abs.size))
    write_spectrum(out / "synthetic_absorption.csv", "absorbance", x_abs, absorbance)
    x_em = np.arange(380.0, 701.0, 1.0)
    emission = gaussian(x_em, 450, 38, 1.0) + rng.normal(0, 0.004, x_em.size)
    write_spectrum(out / "synthetic_emission.csv", "intensity", x_em, emission)
    rows = [("sample_A", 345, 450), ("sample_B", 352, 468), ("sample_C", 330, 415), ("sample_D", 360, 500)]
    with (out / "samples.csv").open("w", newline="") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(["name", "abs_max_nm", "em_max_nm"])
        writer.writerows(rows)
    print(f"Wrote demo files to {out}")


if __name__ == "__main__":
    main()
