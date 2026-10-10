"""Generate the SYNTHETIC demo emission spectra in example/synthetic_pl.csv (fixed seed, reproducible).

Emission bands are Gaussians whose maximum moves with the excitation wavelength:
    emission max = 430 nm + 0.45 * (excitation - 300 nm).
Narrow Rayleigh scatter (at the excitation wavelength) and second-order scatter (at twice that) are added on
purpose, as they appear in real fluorimeter scans. The numbers are invented for the demo, not measured data.
"""
import csv
from pathlib import Path

import numpy as np

SLOPE, EM0, EX0 = 0.45, 430.0, 300.0


def spectrum(wavelength, excitation, rng):
    center = EM0 + SLOPE * (excitation - EX0)
    amplitude = np.exp(-0.5 * ((excitation - 370.0) / 45.0) ** 2)  # brightest near 370 nm excitation
    band = amplitude * np.exp(-0.5 * ((wavelength - center) / 38.0) ** 2)
    rayleigh = 0.9 * np.exp(-0.5 * ((wavelength - excitation) / 4.0) ** 2)
    second = 0.35 * np.exp(-0.5 * ((wavelength - 2 * excitation) / 6.0) ** 2)
    return band + rayleigh + second + rng.normal(0, 0.004, wavelength.size)


def main() -> None:
    rng = np.random.default_rng(14)
    wl = np.arange(320.0, 701.0, 1.0)
    excitations = [300, 320, 340, 360, 380, 400, 420, 440]
    columns = [spectrum(wl, ex, rng) for ex in excitations]
    out = Path(__file__).parent / "example" / "synthetic_pl.csv"
    with out.open("w", newline="") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(["emission_nm"] + [f"ex{ex}" for ex in excitations])
        for k, w in enumerate(wl):
            writer.writerow([f"{w:.0f}"] + [f"{col[k]:.4f}" for col in columns])
    print(f"Wrote {wl.size} rows x {len(excitations)} excitation columns to {out}")


if __name__ == "__main__":
    main()
