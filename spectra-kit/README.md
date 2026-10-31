# spectra-kit

Spectroscopy and materials characterization tools for the command line, packaged from Days 11 to 19 of [100 Days of Computational Chemistry](https://github.com/asifverse4/100-days-of-computational-chemistry). This is mini release **v0.1.0** (Day 20).

Nine tools, each also available as `spectra <tool>`:

```
  uvvis     Day 11  plot UV-vis spectra and find the absorption peaks  (spectra-uvvis)
  tauc      Day 12  Tauc plot and optical band gap from absorbance or reflectance  (spectra-tauc)
  qy        Day 13  relative fluorescence quantum yield against a standard  (spectra-qy)
  pl        Day 14  excitation-dependent emission overlay and emission maximum shift  (spectra-pl)
  stokes    Day 15  Stokes shift in nm, cm-1 and eV  (spectra-stokes)
  ftir      Day 16  find FTIR peaks and suggest functional groups  (spectra-ftir)
  xrd       Day 17  XRD peak fits, Scherrer crystallite size and Williamson-Hall strain  (spectra-xrd)
  tem       Day 18  particle size histogram with a lognormal fit  (spectra-tem)
  baseline  Day 19  smooth spectra and remove baselines  (spectra-baseline)
```

The tools work on the files you already have: a text or CSV file with x in the first column and the signal in the next. Each one writes a plot, a table and, if you ask, a CSV.

## Install
Python 3.10 or newer. The dependencies are numpy, scipy and matplotlib. From a clone of the repository:
```bash
pip install ./spectra-kit
```
Straight from GitHub (this works once the folder is on GitHub; I could not test it before it was pushed):
```bash
pip install "git+https://github.com/asifverse4/100-days-of-computational-chemistry#subdirectory=spectra-kit"
```
It is not published on PyPI.

## Quick start
Copy the example files and try a tool on them:
```bash
spectra examples
spectra-tauc spectra-examples/tauc/synthetic_tauc.csv -o tauc.png
```
```
sample    transition   Eg (eV)  +/- (eV)   fit range (eV)       R2     n
sample_A  direct         3.091    0.0019        3.07-3.32   0.9930    31
sample_B  direct         2.841    0.0018        2.82-3.06   0.9923    36
```
```bash
spectra xrd spectra-examples/xrd/synthetic_xrd.csv --inst-fwhm 0.08 --williamson-hall
```
```
Scherrer size over 6 peaks: mean 12.3 nm, SD 1.5 nm
Williamson-Hall (6 peaks): size 19.6 +/- 1.3 nm, strain 0.00197 +/- 0.00020, R2 0.960
```
```bash
spectra stokes --abs spectra-examples/stokes/synthetic_absorption.csv --em spectra-examples/stokes/synthetic_emission.csv --abs-peak longest
```
```
sample  abs (nm)   em (nm)  shift (nm)  shift (cm-1)  shift (eV)
sample     344.6     450.0       105.4          6794       0.842
```
**All the example data are synthetic.** They were generated from known models with fixed random seeds, so you can check that a tool recovers what was put in. They are not measured data. Use your own files for real work.

The output of one tool can feed another. Remove a baseline, then pick the peaks:
```bash
spectra baseline spectrum.csv --baseline asls --baseline-lambda 1e7 --csv corrected.csv --corrected-only
spectra ftir corrected.csv --signal absorbance
```

## The commands
Each command is the same program as the daily tool, with the same options. Use `--help` for the options, and see the day's README for the method, the checks it makes and its limitations.

| Command | Same as | Day | What it does |
|---|---|---|---|
| `spectra-uvvis` | `spectra uvvis` | [11](../days/day11-uv-vis-plotter/) | plot UV-vis spectra and find the absorption peaks |
| `spectra-tauc` | `spectra tauc` | [12](../days/day12-tauc-band-gap/) | Tauc plot and optical band gap, with a fit uncertainty |
| `spectra-qy` | `spectra qy` | [13](../days/day13-quantum-yield/) | relative fluorescence quantum yield (IUPAC equation) |
| `spectra-pl` | `spectra pl` | [14](../days/day14-pl-excitation-overlay/) | excitation-dependent emission and the shift of the maximum |
| `spectra-stokes` | `spectra stokes` | [15](../days/day15-stokes-shift/) | Stokes shift in nm, cm-1 and eV |
| `spectra-ftir` | `spectra ftir` | [16](../days/day16-ftir-peak-picker/) | FTIR peaks and candidate functional groups |
| `spectra-xrd` | `spectra xrd` | [17](../days/day17-xrd-scherrer/) | peak fits, Scherrer crystallite size, Williamson-Hall strain |
| `spectra-tem` | `spectra tem` | [18](../days/day18-tem-size-distribution/) | TEM size histogram with a lognormal fit |
| `spectra-baseline` | `spectra baseline` | [19](../days/day19-baseline-smoothing/) | smoothing and baseline removal |

`spectra list` shows them, `spectra examples [DIR]` copies the example files (add `--force` to overwrite), and `spectra --version` prints the version.

## How it is built
The daily folders under `days/` are the single source of truth. `spectra-kit/sync.py` copies each `tool.py` into the package, with only its module docstring changed so that `--help` shows `spectra-uvvis` and the packaged example paths instead of `python tool.py` and `example/`. The code itself is not touched. The example files and the license are copied byte for byte.

If a daily tool changes, run `python spectra-kit/sync.py` and commit the result. CI runs `python spectra-kit/sync.py --check` and fails when a copy is out of date, so the package cannot quietly drift from the tools.

The release is checked with `python days/day20-spectra-kit-release/tool.py`, which builds the wheel and the source distribution, looks inside them, installs the wheel into a fresh virtual environment, runs every command and runs the package tests.

## Tests
```bash
pip install ./spectra-kit
pytest spectra-kit
```
The tests check the sync, the registry, the dispatcher, each tool on its packaged example files, that the installed commands and entry points exist and run, and that tools can be chained.

## Limitations
- I ran it on Linux with Python 3.12. The repository's CI uses Python 3.11. It has not been tried on Windows or macOS, and the lower bounds of numpy and scipy are not tested (only matplotlib 3.5 or newer is required, for a feature the tools use).
- The tools are command-line programs that write files. They have no graphical window.
- The method limits described in each day's README apply as they are. In particular, the example data are synthetic, so the recovered values only show that the tools work, not that your measurement is right.
- Version 0.1.0 is a first release. The command options may change in later versions.

## License
MIT, see `LICENSE`.
