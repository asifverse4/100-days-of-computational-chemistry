"""Generate the SYNTHETIC demo FTIR spectrum in example/synthetic_ftir.csv (fixed seed, reproducible).

Absorbance bands (Gaussians) at positions often reported for oxygen- and nitrogen-containing carbon dots:
broad O-H/N-H near 3400, C-H near 2925 and 2855, C=O near 1710, a band near 1620, and bands near 1400,
1260 and 1100, on a sloping baseline with noise, saved as % transmittance in descending wavenumber order
like many instrument exports. The numbers are invented for the demo. They are not measured data.
"""
import csv
from pathlib import Path

import numpy as np

BANDS = [  # (center cm-1, FWHM cm-1, absorbance)
    (3400, 380, 0.55), (2925, 40, 0.18), (2855, 30, 0.12), (1710, 45, 0.28), (1620, 70, 0.38),
    (1400, 80, 0.30), (1260, 70, 0.20), (1100, 120, 0.45),
]


def main() -> None:
    rng = np.random.default_rng(16)
    x = np.arange(4000.0, 399.0, -2.0)
    absorbance = 0.02 + 0.00002 * (4000 - x)  # slowly rising baseline
    for center, fwhm, height in BANDS:
        absorbance += height * np.exp(-0.5 * ((x - center) / (fwhm / 2.3548)) ** 2)
    absorbance += rng.normal(0, 0.002, x.size)
    transmittance = 100 * 10 ** (-absorbance)
    out = Path(__file__).parent / "example" / "synthetic_ftir.csv"
    with out.open("w", newline="") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(["wavenumber_cm-1", "transmittance_percent"])
        for a, t in zip(x, transmittance):
            writer.writerow([f"{a:.0f}", f"{t:.3f}"])
    print(f"Wrote {x.size} points to {out}")


if __name__ == "__main__":
    main()
