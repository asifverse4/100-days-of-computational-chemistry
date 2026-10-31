# Day 20: Release check for spectra-kit v0.1 (mini release)

**What:** Today closes Phase 2 with the first mini release, **spectra-kit v0.1.0**: the nine spectroscopy tools from Days 11 to 19 as an installable package ([`spectra-kit/`](../../spectra-kit/)). The tool for the day is a release checker. It builds the package, looks inside what it built, installs it into a fresh virtual environment, runs every command and the package tests, and prints the git tag command only if everything passed.
**Why:** A release is easy to get subtly wrong: a version that disagrees between files, a copy that no longer matches the tool it came from, a missing data file in the wheel, a command that fails on import, a dependency that was never declared. A person checking by hand misses things. This turns the checks into one repeatable command that CI also runs.

## Install
```bash
pip install build
```
The check needs Python 3.11 or newer. Building downloads the build backend (setuptools) unless you pass `--no-isolation`.

## Usage
From the repository root:
```bash
python days/day20-spectra-kit-release/tool.py
```
```
[PASS] version format: spectra-kit 0.1.0
[PASS] __version__ matches pyproject: 0.1.0
[PASS] README is present
[PASS] license file is present
[PASS] console scripts are declared: 10 commands
[PASS] CHANGELOG mentions this release
[PASS] git tag spectra-kit-v0.1.0 is free
[PASS] copies are in step with the daily tools
[PASS] build the wheel and the source distribution: spectra_kit-0.1.0-py3-none-any.whl, spectra_kit-0.1.0.tar.gz
[PASS] wheel metadata has the right name and version
[PASS] wheel carries the license
[PASS] wheel entry points match pyproject: 10 commands
[PASS] wheel has no tests or build scripts
[PASS] wheel includes its data files: 14 CSV files
[PASS] source distribution has pyproject, license and README
[PASS] source distribution has no tests that cannot run outside the repository
[PASS] create a fresh virtual environment: dependencies taken from this Python
[PASS] install the wheel: spectra_kit-0.1.0-py3-none-any.whl
[PASS] every installed command runs (--help), 10 commands
[PASS] package tests pass against the installed wheel: 36 passed in 13.81s

20 checks, 0 failed
spectra-kit 0.1.0 is ready to tag:
  git tag -a spectra-kit-v0.1.0 -m "spectra-kit 0.1.0"
  git push origin spectra-kit-v0.1.0
  then on GitHub: Releases > Draft a new release > choose spectra-kit-v0.1.0, and paste the spectra-kit v0.1.0 notes from CHANGELOG.md
```
(`example/release_check_output.txt` holds this output from a real run.) A failed check is printed as `[FAIL]` with the reason, the exit code is 1, and the tag command is not shown.

| Option | Meaning |
|---|---|
| `PACKAGE` | the package folder under the repository root (default `spectra-kit`) |
| `--clean-venv` | build the test environment without this Python's packages, so numpy, scipy and matplotlib come from PyPI and pytest is installed. This is the honest test that the package declares all its own dependencies |
| `--no-isolation` | build with the tools already installed instead of downloading the build backend |
| `--skip-tests` | do not run the package tests |
| `--keep DIR` | keep the wheel and the source distribution in `DIR` |
| `--allow-existing-tag` | do not fail when the tag already exists (CI uses this, because after you tag the release the tag always exists) |
| `--tag-prefix` | the tag prefix (default `<package name>-v`) |
| `--root` | the repository root (default: this repository) |

## What it checks
1. **Metadata:** the version is `major.minor.patch`, `__version__` in the package equals the version in `pyproject.toml`, the README and license files exist, console scripts are declared, `CHANGELOG.md` has release notes that name this version (`spectra-kit v0.1.0`), and the git tag does not exist yet.
2. **Copies in step:** if the package has a `sync.py`, run its `--check`. For spectra-kit this proves the package copies match the daily tools.
3. **Build:** the wheel and the source distribution build.
4. **Inside the wheel:** the metadata has the right name and version, the license is included, the entry points are exactly the ones declared in `pyproject.toml`, no tests or build scripts leaked in, and the data files are there when the package declares them. The source distribution has `pyproject.toml`, the license and the README.
5. **Install:** a fresh virtual environment, then `pip install` of the wheel.
6. **Run:** every installed command answers `--help`, and the package tests pass against the installed wheel.

## Method
Every step that runs a program goes through one function, `run_command`, which returns the exit code and the output and never raises. This is what lets the tests replace it with a fake runner: the tests build a believable wheel and source distribution, record the commands, and check the logic of the whole run (which steps run, which are skipped after a failure, what changes with each option) without a real build. A separate test, skipped unless you set `RELEASE_CHECK_FULL=1`, runs the real thing.

## Test
```bash
pytest days/day20-spectra-kit-release
RELEASE_CHECK_FULL=1 pytest days/day20-spectra-kit-release     # also builds and installs the real package
```
The tests cover every metadata problem, the git tag check with a real git repository, the sync check, good and bad wheels and source distributions, every option, a failing step at every stage, the usage errors, and the static checks of the real spectra-kit.

## Limitations
- **By default the environment reuses this Python's packages**, so it can miss a dependency the package forgot to declare. Run with `--clean-venv` before you publish. I ran it that way once for v0.1.0 and it passed.
- The commands are only run with `--help` here. Their real behavior is covered by the package tests (each tool on its example files), which is why running those matters.
- The tag check only sees tags in your local clone. Run `git fetch --tags` first if someone else may have tagged.
- It does not upload anything, create a tag or make a GitHub release. It checks and prints the commands for you to run.
- It does not test other operating systems or Python versions. It checks the one Python it runs on.
- The tool needs Python 3.11 (it reads TOML with `tomllib`) and network access for the build, unless you use `--no-isolation`.
- The first CI run on GitHub is the real test of the CI step. I could only run the same commands here.

## Next
Day 21 starts Phase 3, DFT and xTB workflows, with an ORCA input generator.
