---
name: new-day
description: Builds one day of the 100 Days of Computational Chemistry repo end to end. Scaffolds the day folder, writes tool.py, tests, README and sourced or clearly synthetic example data, updates the README tables, CHANGELOG and INDEX.md, then runs ruff, pytest and the repo check. Use when asked to add, start, continue or fix a day, for example "do Day 18".
license: MIT
compatibility: Needs git and Python 3.11 or newer with numpy, scipy, pandas, matplotlib, rdkit, pytest and ruff (see requirements.txt).
metadata:
  version: "1.0"
---

# Build one day of 100 Days of Computational Chemistry

Read `AGENTS.md` first. It holds the rules. This skill is the order of work.

## 1. Understand the day
- Find the day in `ROADMAP.md`. Note its phase and the tool it names.
- Look at the previous two days in `days/` to match the style (README sections, CLI options, test layout).
- Decide inputs, outputs and the method. If a published equation or a reference table is involved, find the source before writing code.

## 2. Scaffold
```bash
python tools/new_day.py N slug      # for example: python tools/new_day.py 18 tem-size-distribution
```
This copies `template/dayNN-name/` and updates the README Progress table and the badge. Delete the template's placeholder text.

## 3. Data first
- Reference tables: cite the source, check every value against it, tag what you could not check.
- Example data: write `make_example.py` with a fixed seed. Label the data **synthetic** in its docstring and in the README. Use real, sourced data only if you have it.

## 4. Write `tool.py`
- Docstring with usage lines, `argparse`, `main() -> int`, exit codes 0, 1 and 2, errors on stderr.
- Standalone: no imports from other days.
- Units in names. State what each uncertainty includes. Give real measurement uncertainties to fits (`absolute_sigma=True`).
- Handle bad input with a clear message, not a traceback.

## 5. Write `test_dayNN.py`
- Hand-check each formula, recover known values from synthetic data, and compare reported uncertainties with a simulation.
- Cover every error message and exit code, the plot formats and the CSV columns.
- If a test fails, work out whether the tool or the test is wrong. Report which it was.

## 6. Write the README
Title `# Day N: Title`, `**What:**`, `**Why:**`, then Install, Usage (real output pasted from a run, with an image if the tool makes one), Method, Test, Limitations and Next. Look at any plot you commit before keeping it.

## 7. Repo files
- Add rows to the tables in `README.md`, and a line to `CHANGELOG.md`.
- `python days/day10-index-page/tool.py`, then `python days/day10-index-page/tool.py --check`.

## 8. Verify
```bash
ruff check .
for d in days/*/; do pytest -q "$d"; done
python days/day10-index-page/tool.py --check
```
Everything must pass. CI also audits `tool.py` files for unit and uncertainty defects.

## 9. Report
Commit as `Day N: <tool>`. In the summary, list what you ran, which bugs the tests found, and what you could not verify (real instrument file formats, reference values, behavior on other platforms). Do not claim more than you checked.
