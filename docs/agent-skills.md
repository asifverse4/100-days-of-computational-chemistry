# Agent skills for this repo

An **agent skill** is a folder of instructions (a `SKILL.md` file, sometimes with scripts) that an AI coding agent loads when a task matches. The [Scientific Agent Skills](https://github.com/K-Dense-AI/scientific-agent-skills) collection from K-Dense (MIT license, 177 skills when this page was written) has several that fit the work in this repo. This page says which ones, and how to use them safely.

## Install safely
Skills can run code and change how an agent behaves. The collection's own security disclaimer says it plainly: "Skills can execute code and influence your coding agent's behavior. Review what you install." Before installing one, read its `SKILL.md` and any scripts, and **pin a version** so it cannot change under you.

```bash
npx skills add K-Dense-AI/scientific-agent-skills
gh skill install K-Dense-AI/scientific-agent-skills rdkit --agent claude-code      # one skill, for one agent
gh skill install K-Dense-AI/scientific-agent-skills --pin <tag-or-commit>          # pin a version
gh skill update
```
`gh skill` needs GitHub CLI 2.90 or newer. Supported agents include `claude-code`, `cursor`, `codex` and `gemini`. See the collection's README for the other options and install paths.

## Which skills fit which part of the roadmap
These matches come from each skill's description. They are suggestions, and I have not tried every skill with this repo.

| Phase (days) | Skills | Useful for |
|---|---|---|
| 1 Foundations (1-10) | `rdkit`, `datamol`, `medchem`, `uncertainty-and-units` | SMILES, descriptors and 3D (Days 2, 3, 9). `medchem` adds drug-likeness rules and structural alerts beyond Lipinski. `uncertainty-and-units` for unit conversion and constants (Day 7) |
| 2 Spectroscopy (11-20) | `scientific-visualization`, `matplotlib`, `uncertainty-and-units`, `statistical-analysis` | figure design and checks for color and contrast, uncertainty budgets for band gap, quantum yield and crystallite size, fit and test selection |
| 3 DFT and xTB (21-30) | `pymatgen`, `rowan` | `pymatgen` for periodic structures and file conversion. `rowan` is a cloud platform that needs an account and sends structures to an outside service, so do not use it with unpublished structures unless that is acceptable |
| 4 Docking and MD (31-40) | `molecular-dynamics`, `diffdock` | OpenMM and MDAnalysis analysis (Day 37), and deep-learning docking as an alternative to the AutoDock Vina route |
| 5 Synthesis (41-50) | `experimental-design`, `statistical-analysis`, `datamol` | screening designs (Day 44), yield statistics (Day 41), reaction and structure handling (Day 48) |
| 6 Machine learning (51-60) | `scikit-learn`, `deepchem`, `molfeat`, `shap`, `torch-geometric` | fingerprints and featurization, scaffold splits, property models, feature importance, graph networks |
| 8 Apps and figures (71-80) | `scientific-visualization`, `matplotlib` | publication figures and palettes |
| 9-10 Quality and capstone (81-100) | `uncertainty-and-units`, `scientific-writing`, `citation-management`, `scientific-slides` | benchmarks against literature values, the blog post and thesis text, citations, talks |
| Later days that read raw data | `nmrglue`, `matchms` | NMR free-induction decays (Day 47) and tandem mass spectra |

## How this repo already uses the collection
- **CI audit.** `.github/workflows/ci.yml` runs the `uncertainty-and-units` auditor (`audit_units.py`) over every `days/*/tool.py`. It reads the code without running it, and flags stripped units, broken uncertainty propagation and hard-coded constants. CI fetches the script from one **pinned commit** (`a1043ef544e3a100c28e70a736b6cf32554d0677`, which I read through before pinning: it has no network or subprocess calls). To update it, review the new version, then change the commit in the workflow.
- **What the audit found.** On first use it flagged one high-severity item in the XRD tool (Day 17): the fit uncertainties were not based on the counting noise of the data. A simulation showed the old uncertainties were about 43% too small, and the tool now weights raw counts by their Poisson noise. It also listed five hard-coded physical constants in the unit converter (Day 7). They are exact SI-defining values, and the file now says so, with a `# audit-units: ignore-file CONST001 -- reason` note that the report counts as a suppression.
- **Silencing a finding.** Only after a real review, with a specific reason on the same line: `# audit-units: ignore RULE -- reason`.

## Using skills with this repo's own conventions
`AGENTS.md` holds the rules for every agent. The project skill `.claude/skills/new-day/` lists the steps for building a day. Claude Code reads `.claude/skills/`. For other agents, copy the folder to `.agents/skills/` or install it with the tool you use. The repo check (`days/day10-index-page/tool.py --check`) validates the format of any skill in either place.

A reasonable request to an agent that has `uncertainty-and-units` installed: "Review `days/day13-quantum-yield/tool.py` for unit and uncertainty problems and show me a simulation that checks the reported uncertainty."
