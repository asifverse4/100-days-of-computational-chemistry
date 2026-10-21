"""Generate the SYNTHETIC demo spectrum in example/ (fixed seed, reproducible).

synthetic_spectrum.csv: six Gaussian peaks (positions 520 to 1620, widths 20 to 40) on a smooth curved baseline
that looks like a fluorescence background (a slow rise plus a broad hump), plus noise (sigma 0.006).
synthetic_baseline_truth.csv: the true baseline, so methods can be compared with --truth.
The numbers are invented for the demo. They are not measured data.
"""
import csv
from pathlib import Path

import numpy as np


def main() -> None:
    rng = np.random.default_rng(19)
    x = np.arange(400.0, 1801.0, 1.0)
    peaks = [(520, 22, 0.55), (760, 18, 0.30), (1004, 14, 1.00), (1250, 28, 0.45), (1450, 35, 0.70), (1620, 24, 0.50)]
    signal = sum(h * np.exp(-0.5 * ((x - c) / w) ** 2) for c, w, h in peaks)
    u = (x - 400.0) / 1400.0
    baseline = 0.15 + 0.55 * u - 0.35 * u**2 + 0.30 * np.exp(-0.5 * ((x - 1200.0) / 380.0) ** 2)
    y = signal + baseline + rng.normal(0, 0.006, x.size)
    out = Path(__file__).parent / "example"
    for name, header, values in (("synthetic_spectrum.csv", "intensity", y), ("synthetic_baseline_truth.csv", "baseline", baseline)):
        with (out / name).open("w", newline="") as fh:
            writer = csv.writer(fh, lineterminator="\n")
            writer.writerow(["wavenumber_cm-1", header])
            for a, b in zip(x, values):
                writer.writerow([f"{a:.0f}", f"{b:.5f}"])
    print(f"Wrote {x.size} points to {out}")


if __name__ == "__main__":
    main()
