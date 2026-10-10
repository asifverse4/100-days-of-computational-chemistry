"""Generate the SYNTHETIC demo spectra in example/synthetic_tauc.csv (fixed seed, reproducible).

Both samples follow the direct-gap Tauc law, (alpha*h*nu)^2 = B*(h*nu - Eg), above their gap, with a
small exponential tail below it, a constant baseline offset and noise. The gaps (3.10 and 2.85 eV)
are chosen for the demo. They are not measured values.
"""
import csv
from pathlib import Path

import numpy as np

HC_EV_NM = 1239.841984  # h*c in eV nm


def absorbance(wavelength_nm, gap_ev, b=170.0, offset=0.012, noise=0.002, seed=0):
    e = HC_EV_NM / wavelength_nm
    alpha = np.sqrt(b * np.clip(e - gap_ev, 0, None)) / e  # cm^-1, path length 1 cm
    tail = 0.004 * np.exp((e - gap_ev) / 0.06) * (e <= gap_ev)  # Urbach-like tail below the gap
    rng = np.random.default_rng(seed)
    return (alpha + tail) / 2.303 + offset + rng.normal(0, noise, e.size)


def main() -> None:
    wl = np.arange(250.0, 701.0, 1.0)
    a, b = absorbance(wl, 3.10, seed=1), absorbance(wl, 2.85, seed=2)
    out = Path(__file__).parent / "example" / "synthetic_tauc.csv"
    with out.open("w", newline="") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(["wavelength_nm", "sample_A", "sample_B"])
        for row in zip(wl, a, b):
            writer.writerow([f"{row[0]:.0f}", f"{row[1]:.4f}", f"{row[2]:.4f}"])
    print(f"Wrote {wl.size} rows to {out}")


if __name__ == "__main__":
    main()
