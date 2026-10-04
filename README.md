# 100 Days of Computational Chemistry

One small, working, documented tool per day for chemistry students and researchers: spectroscopy, DFT/xTB, docking/MD, synthesis analysis, and ML.

![streak](https://img.shields.io/badge/streak-day%204%2F100-blue)

## Quick start
```bash
git clone https://github.com/asifverse4/100-days-of-computational-chemistry
cd 100-days-of-computational-chemistry
pip install -r requirements.txt
python tools/new_day.py 2 smiles-to-3d      # scaffolds days/day02-smiles-to-3d/
```

## How it works
- Each day lives in `days/dayNN-name/` with a README, a script, and a tiny example.
- `ROADMAP.md` has the full 100-day plan.
- Every tenth day ends in a mini release.

## Progress
| Day | Topic | Folder |
|-----|-------|--------|
| 1 | Repo setup | [`template/`](template/) |
| 2 | Smiles To 3D | [`days/day02-smiles-to-3d/`](days/day02-smiles-to-3d/) |
| 3 | Descriptors Table | [`days/day03-descriptors-table/`](days/day03-descriptors-table/) |
| 4 | xTB Batch | [`days/day04-xtb-batch/`](days/day04-xtb-batch/) |

## Contributing
Issues and PRs are welcome. See `CONTRIBUTING.md`.

## Citation
If a tool helps your work, please cite via `CITATION.cff`.

## License
MIT
