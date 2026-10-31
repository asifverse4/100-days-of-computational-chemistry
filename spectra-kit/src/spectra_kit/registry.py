"""Which daily tools make up spectra-kit, and where each one comes from.

The daily folders under days/ stay the single source of truth. spectra-kit/sync.py copies each tool.py
into this package, and CI fails if a copy is out of date.
"""
from dataclasses import dataclass


@dataclass(frozen=True)
class Tool:
    key: str  # the name after "spectra " or "spectra-"
    module: str  # module name inside spectra_kit
    day: int
    source: str  # folder under the repository root
    title: str
    summary: str
    examples: tuple[str, ...]  # files from the day's example/ folder that ship with the package (all synthetic)


TOOLS = (
    Tool("uvvis", "uvvis", 11, "days/day11-uv-vis-plotter", "UV-vis plotter",
         "plot UV-vis spectra and find the absorption peaks", ("synthetic_cqd.csv",)),
    Tool("tauc", "tauc", 12, "days/day12-tauc-band-gap", "Tauc band gap",
         "Tauc plot and optical band gap from absorbance or reflectance", ("synthetic_tauc.csv",)),
    Tool("qy", "quantum_yield", 13, "days/day13-quantum-yield", "Quantum yield",
         "relative fluorescence quantum yield against a standard", ("synthetic_reference.csv", "synthetic_sample.csv")),
    Tool("pl", "pl_overlay", 14, "days/day14-pl-excitation-overlay", "PL excitation overlay",
         "excitation-dependent emission overlay and emission maximum shift", ("synthetic_pl.csv",)),
    Tool("stokes", "stokes", 15, "days/day15-stokes-shift", "Stokes shift",
         "Stokes shift in nm, cm-1 and eV", ("synthetic_absorption.csv", "synthetic_emission.csv", "samples.csv")),
    Tool("ftir", "ftir", 16, "days/day16-ftir-peak-picker", "FTIR peak picker",
         "find FTIR peaks and suggest functional groups", ("synthetic_ftir.csv",)),
    Tool("xrd", "xrd", 17, "days/day17-xrd-scherrer", "XRD Scherrer",
         "XRD peak fits, Scherrer crystallite size and Williamson-Hall strain", ("synthetic_xrd.csv",)),
    Tool("tem", "tem_sizes", 18, "days/day18-tem-size-distribution", "TEM size distribution",
         "particle size histogram with a lognormal fit", ("synthetic_sizes.csv", "synthetic_imagej.csv")),
    Tool("baseline", "baseline", 19, "days/day19-baseline-smoothing", "Baseline and smoothing",
         "smooth spectra and remove baselines", ("synthetic_spectrum.csv", "synthetic_baseline_truth.csv")),
)

BY_KEY = {t.key: t for t in TOOLS}


def command(tool: Tool) -> str:
    """The stand-alone command name, for example spectra-uvvis."""
    return f"spectra-{tool.key}"
