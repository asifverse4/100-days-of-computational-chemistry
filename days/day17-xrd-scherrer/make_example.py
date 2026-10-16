"""Generate the SYNTHETIC demo XRD pattern in example/synthetic_xrd.csv (fixed seed, reproducible).

A face-centered cubic crystal (a = 4.05 A) measured with Cu K-alpha1. The peak widths follow the
Williamson-Hall relation  beta*cos(theta) = K*lambda/D + 4*strain*sin(theta)  with D = 20 nm and strain = 0.002,
combined in quadrature with an instrumental width of 0.08 degrees, and pseudo-Voigt peak shapes (eta = 0.5),
a sloping background and counting noise. The numbers are invented for the demo. They are not measured data.
"""
import csv
from pathlib import Path

import numpy as np

WAVELENGTH, K, SIZE_NM, STRAIN, INSTRUMENT_DEG, A_LATTICE = 1.54056, 0.9, 20.0, 0.002, 0.08, 4.05
REFLECTIONS = [((1, 1, 1), 1000), ((2, 0, 0), 450), ((2, 2, 0), 250), ((3, 1, 1), 280), ((2, 2, 2), 110), ((4, 0, 0), 90)]


def pseudo_voigt(x, center, fwhm, eta=0.5):
    u = (x - center) / fwhm
    return eta / (1 + 4 * u**2) + (1 - eta) * np.exp(-4 * np.log(2) * u**2)


def main() -> None:
    rng = np.random.default_rng(17)
    x = np.arange(20.0, 100.0, 0.02)
    signal = 60 + 0.4 * (100 - x)  # sloping background
    for (h, k, l), height in REFLECTIONS:
        d = A_LATTICE / np.sqrt(h * h + k * k + l * l)
        theta = np.arcsin(WAVELENGTH / (2 * d))
        beta = (K * WAVELENGTH / (SIZE_NM * 10) + 4 * STRAIN * np.sin(theta)) / np.cos(theta)  # radians
        fwhm = np.hypot(np.degrees(beta), INSTRUMENT_DEG)
        signal = signal + height * pseudo_voigt(x, np.degrees(2 * theta), fwhm)
    counts = rng.poisson(signal)
    out = Path(__file__).parent / "example" / "synthetic_xrd.csv"
    with out.open("w", newline="") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(["two_theta_deg", "intensity_counts"])
        for a, c in zip(x, counts):
            writer.writerow([f"{a:.2f}", int(c)])
    print(f"Wrote {x.size} points to {out}")


if __name__ == "__main__":
    main()
