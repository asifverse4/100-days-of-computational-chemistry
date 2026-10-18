# Changelog

## Repo tooling
- `AGENTS.md`, the `new-day` project skill, a pull request template and a "New day or tool" issue form, so AI coding agents follow the repo's rules.
- `docs/agent-skills.md`: which K-Dense Scientific Agent Skills fit each phase, and how to install them safely.
- CI audit of every `tool.py` with the pinned `uncertainty-and-units` auditor. It found a real problem in Day 17 (fit uncertainties about 43% too small in a simulation), which is fixed: raw-count patterns are now fitted with counting-statistics weights.
- The Day 10 check now validates `AGENTS.md` and the format of any project skill.

## Phase 2: Nanomaterial and spectroscopy data (Days 11-20, in progress)
- Day 11: UV-vis plotter with peak detection.
- Day 12: Tauc plot and optical band gap calculator.
- Day 13: relative fluorescence quantum yield calculator (IUPAC equation, gradient and single-point methods).
- Day 14: PL spectrum plotter with excitation-dependent emission overlay and scatter masking.
- Day 15: Stokes shift calculator (nm, cm-1, eV) with band selection.
- Day 16: FTIR peak picker with functional-group suggestions.
- Day 17: XRD Scherrer crystallite size calculator with Williamson-Hall strain separation.

## Phase 1: Foundations (Days 1-10)
- Day 1: repo setup with license, README, requirements and the day template.
- Day 2: SMILES to optimized 3D `.xyz` geometries (RDKit).
- Day 3: molecular descriptors table with Lipinski rule-of-five checks.
- Day 4: batch xTB optimization wrapper.
- Day 5: xTB output parser for energies, gaps and charges, with conformer ranking.
- Day 6: XYZ toolkit to inspect, split, merge, center, align and convert structure files.
- Day 7: unit converter (hartree, eV, kcal/mol, kJ/mol, cm-1, nm, THz, K).
- Day 8: offline periodic-table lookup and formula molar mass.
- Day 9: molecule grid image generator for figures.
- Day 10: index page builder and repo structure check in CI.