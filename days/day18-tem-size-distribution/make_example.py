"""Generate the SYNTHETIC demo data in example/ (fixed seed, reproducible).

synthetic_sizes.csv: 250 particle diameters drawn from a lognormal distribution (median 4.2 nm,
geometric standard deviation 1.32), rounded to 0.1 nm like a manual measurement.
synthetic_imagej.csv: 220 particles in the layout of an ImageJ "Results" table (Area in square pixels,
0.25 nm per pixel), drawn from a different lognormal distribution (median 6.0 nm, geometric SD 1.20).
The numbers are invented for the demo. They are not measured data.
"""
import csv
from pathlib import Path

import numpy as np


def main() -> None:
    rng = np.random.default_rng(18)
    out = Path(__file__).parent / "example"

    sizes = np.round(rng.lognormal(np.log(4.2), np.log(1.32), 250), 1)
    with (out / "synthetic_sizes.csv").open("w", newline="") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(["particle", "diameter_nm"])
        for i, d in enumerate(sizes, 1):
            writer.writerow([i, f"{d:.1f}"])

    nm_per_pixel = 0.25
    diameters_nm = rng.lognormal(np.log(6.0), np.log(1.20), 220)
    area_px2 = np.pi * (diameters_nm / nm_per_pixel / 2) ** 2
    with (out / "synthetic_imagej.csv").open("w", newline="") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow([" ", "Area", "Mean", "Perim."])  # ImageJ leaves the first header cell blank
        for i, a in enumerate(area_px2, 1):
            writer.writerow([i, f"{a:.0f}", f"{rng.uniform(80, 200):.3f}", f"{np.pi * np.sqrt(4 * a / np.pi):.2f}"])
    print(f"Wrote {sizes.size} sizes and {area_px2.size} ImageJ rows to {out}")


if __name__ == "__main__":
    main()
