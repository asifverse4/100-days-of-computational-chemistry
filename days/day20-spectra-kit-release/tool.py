"""Day 20: release check for a mini-release package (spectra-kit v0.1).

Before you tag a release: check the package metadata and version, check that its copies of the daily tools are
in step, build the wheel and the source distribution, look inside them, install the wheel into a fresh virtual
environment, run every installed command and the package tests, and print the tag command.

Needs Python 3.11 or newer and the `build` package (pip install build). Building downloads the build backend
unless you pass --no-isolation.

Usage:
    python days/day20-spectra-kit-release/tool.py                      # check spectra-kit
    python days/day20-spectra-kit-release/tool.py --clean-venv         # also install the dependencies from PyPI
    python days/day20-spectra-kit-release/tool.py --skip-tests --keep dist
"""
import argparse
import os
import re
import subprocess
import sys
import tarfile
import tempfile
import tomllib
import zipfile
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
VERSION_RE = re.compile(r"^\d+\.\d+\.\d+$")
BUILD_TIMEOUT, VENV_TIMEOUT, SMOKE_TIMEOUT, TEST_TIMEOUT = 900, 1800, 180, 900


@dataclass
class Check:
    name: str
    ok: bool
    detail: str = ""
    warning: bool = False  # a problem worth showing that does not fail the release

    def line(self) -> str:
        status = "PASS" if self.ok else ("warn" if self.warning else "FAIL")
        return f"[{status}] {self.name}" + (f": {self.detail}" if self.detail else "")


def run_command(cmd: list[str], cwd: Path | None = None, env: dict | None = None, timeout: int = 600) -> tuple[int, str]:
    """Run a command and return (exit code, combined output). Tests replace this function."""
    try:
        done = subprocess.run(cmd, cwd=cwd, env=env, capture_output=True, text=True, timeout=timeout, check=False)
    except FileNotFoundError as err:
        return 127, str(err)
    except subprocess.TimeoutExpired:
        return 124, f"timed out after {timeout} s"
    return done.returncode, (done.stdout + done.stderr).strip()


def tail(text: str, lines: int = 4) -> str:
    return " | ".join(text.strip().splitlines()[-lines:])


# ---------- metadata and version ----------
def read_pyproject(package: Path) -> dict:
    """The [project] table of the package's pyproject.toml."""
    path = package / "pyproject.toml"
    if not path.is_file():
        raise ValueError(f"{path} not found")
    try:
        document = tomllib.loads(path.read_text())
    except tomllib.TOMLDecodeError as err:
        raise ValueError(f"{path} is not valid TOML: {err}") from None
    project = dict(document.get("project", {}))
    for field in ("name", "version"):
        if field not in project:
            raise ValueError(f"{path} has no project.{field}")
    # data files are declared under [tool.setuptools], not [project]: keep the answer with the project table
    project["_declares_package_data"] = bool(document.get("tool", {}).get("setuptools", {}).get("package-data"))
    return project


def module_version(package: Path, name: str) -> str | None:
    """The __version__ string in src/<module>/__init__.py, if there is one."""
    init = package / "src" / name.replace("-", "_") / "__init__.py"
    if not init.is_file():
        return None
    found = re.search(r'^__version__\s*=\s*["\']([^"\']+)["\']', init.read_text(), re.MULTILINE)
    return found[1] if found else None


def tag_exists(root: Path, tag: str) -> bool | None:
    """True or False from git, or None when git or the repository is not available."""
    code, out = run_command(["git", "tag", "-l", tag], cwd=root, timeout=30)
    return None if code != 0 else tag in out.split()


def check_metadata(package: Path, root: Path, tag_prefix: str | None = None, allow_existing_tag: bool = False) -> list[Check]:
    project = read_pyproject(package)
    name, version = project["name"], project["version"]
    tag = f"{tag_prefix if tag_prefix is not None else name + '-v'}{version}"
    checks = [Check("version format", bool(VERSION_RE.match(version)), f"{name} {version}" if VERSION_RE.match(version)
                    else f"{version!r} is not major.minor.patch")]
    in_module = module_version(package, name)
    checks.append(Check("__version__ matches pyproject", in_module == version,
                        f"{in_module!r} in the package, {version!r} in pyproject" if in_module != version else version))
    readme = package / str(project.get("readme", "README.md"))
    checks.append(Check("README is present", readme.is_file(), "" if readme.is_file() else f"{readme.name} is missing"))
    licences = [package / f for f in project.get("license-files", [])]
    missing = [p.name for p in licences if not p.is_file()]
    checks.append(Check("license file is present", bool(licences) and not missing,
                        f"missing {', '.join(missing)}" if missing else ("no license-files in pyproject" if not licences else "")))
    checks.append(Check("console scripts are declared", bool(project.get("scripts")),
                        f"{len(project.get('scripts', {}))} commands" if project.get("scripts") else "none"))
    changelog = root / "CHANGELOG.md"
    noted = changelog.is_file() and f"{name} v{version}" in changelog.read_text()
    checks.append(Check("CHANGELOG mentions this release", noted, f"expected the text '{name} v{version}'" if not noted else ""))
    exists = tag_exists(root, tag)
    if exists is None:
        checks.append(Check(f"git tag {tag} is free", True, "could not ask git, not checked", warning=True))
    elif exists and allow_existing_tag:
        checks.append(Check(f"git tag {tag} is free", False, "the tag already exists (allowed)", warning=True))
    else:
        checks.append(Check(f"git tag {tag} is free", not exists, "the tag already exists" if exists else ""))
    return checks


def check_sync(package: Path) -> Check:
    """Run the package's own sync check, if it has one."""
    script = package / "sync.py"
    if not script.is_file():
        return Check("copies are in step with the daily tools", True, "no sync.py in the package, nothing to check", warning=True)
    code, out = run_command([sys.executable, str(script), "--check"], timeout=120)
    return Check("copies are in step with the daily tools", code == 0, tail(out) if code else "")


# ---------- build and inspect ----------
def build_package(package: Path, outdir: Path, isolation: bool = True) -> tuple[list[Check], Path | None, Path | None]:
    cmd = [sys.executable, "-m", "build", "--outdir", str(outdir), str(package)] + ([] if isolation else ["--no-isolation"])
    code, out = run_command(cmd, timeout=BUILD_TIMEOUT)
    wheels, sdists = sorted(outdir.glob("*.whl")), sorted(outdir.glob("*.tar.gz"))
    if code != 0 or not wheels or not sdists:
        return [Check("build the wheel and the source distribution", False, tail(out))], None, None
    return [Check("build the wheel and the source distribution", True, f"{wheels[0].name}, {sdists[0].name}")], wheels[0], sdists[0]


def parse_entry_points(text: str) -> dict[str, str]:
    """The [console_scripts] section of an entry_points.txt file."""
    scripts, inside = {}, False
    for line in text.splitlines():
        line = line.strip()
        if line.startswith("["):
            inside = line == "[console_scripts]"
        elif inside and "=" in line:
            key, value = line.split("=", 1)
            scripts[key.strip()] = value.strip()
    return scripts


def inspect_wheel(wheel: Path, project: dict) -> list[Check]:
    name, version = project["name"], project["version"]
    with zipfile.ZipFile(wheel) as zf:
        names = zf.namelist()
        dist = f"{name.replace('-', '_')}-{version}.dist-info"
        metadata = zf.read(f"{dist}/METADATA").decode() if f"{dist}/METADATA" in names else ""
        entry = zf.read(f"{dist}/entry_points.txt").decode() if f"{dist}/entry_points.txt" in names else ""
    checks = [Check("wheel metadata has the right name and version",
                    f"Name: {name}" in metadata and f"Version: {version}" in metadata,
                    "" if metadata else "no METADATA file")]
    checks.append(Check("wheel carries the license", any(n.startswith(f"{dist}/licenses/") or n == f"{dist}/LICENSE" for n in names)))
    expected = dict(project.get("scripts", {}))
    found = parse_entry_points(entry)
    checks.append(Check("wheel entry points match pyproject", found == expected,
                        f"missing {sorted(set(expected) - set(found))}, unexpected {sorted(set(found) - set(expected))}"
                        if found != expected else f"{len(found)} commands"))
    leaked = [n for n in names if n.startswith("tests/") or "/tests/" in n or n.endswith("sync.py")]
    checks.append(Check("wheel has no tests or build scripts", not leaked, ", ".join(leaked[:3])))
    data = [n for n in names if n.endswith(".csv")]
    declared = project.get("_declares_package_data", False)
    checks.append(Check("wheel includes its data files", bool(data) or not declared,
                        f"{len(data)} CSV files" if data else ("package-data is declared but the wheel has no data files" if declared else "")))
    return checks


def inspect_sdist(sdist: Path, project: dict) -> list[Check]:
    with tarfile.open(sdist) as tf:
        names = tf.getnames()
    base = f"{project['name'].replace('-', '_')}-{project['version']}"
    needed = [f"{base}/pyproject.toml", f"{base}/LICENSE", f"{base}/README.md"]
    missing = [Path(n).name for n in needed if n not in names]
    checks = [Check("source distribution has pyproject, license and README", not missing,
                    f"missing {', '.join(missing)}" if missing else "")]
    tests = [n for n in names if "/tests/" in n]
    checks.append(Check("source distribution has no tests that cannot run outside the repository", not tests,
                        f"{len(tests)} test files shipped" if tests else "", warning=True))
    return checks


# ---------- install and run ----------
def venv_paths(venv: Path, windows: bool | None = None) -> tuple[Path, Path]:
    """(scripts folder, python executable) of a virtual environment."""
    windows = (os.name == "nt") if windows is None else windows
    folder = venv / ("Scripts" if windows else "bin")
    return folder, folder / ("python.exe" if windows else "python")


def install_wheel(venv: Path, wheel: Path, clean: bool) -> list[Check]:
    code, out = run_command([sys.executable, "-m", "venv", *([] if clean else ["--system-site-packages"]), str(venv)],
                            timeout=VENV_TIMEOUT)
    if code != 0:
        return [Check("create a fresh virtual environment", False, tail(out))]
    _, python = venv_paths(venv)
    code, out = run_command([str(python), "-m", "pip", "install", "--quiet", str(wheel)], timeout=VENV_TIMEOUT)
    mode = "dependencies installed from PyPI" if clean else "dependencies taken from this Python"
    return [Check("create a fresh virtual environment", True, mode),
            Check("install the wheel", code == 0, tail(out) if code else wheel.name)]


def install_pytest(venv: Path) -> Check:
    """A clean environment has no pytest: install it so the package tests can run."""
    _, python = venv_paths(venv)
    code, out = run_command([str(python), "-m", "pip", "install", "--quiet", "pytest"], timeout=VENV_TIMEOUT)
    return Check("install pytest into the clean environment", code == 0, tail(out) if code else "")


def smoke_scripts(venv: Path, scripts: list[str]) -> list[Check]:
    folder, _ = venv_paths(venv)
    bad = []
    for script in scripts:
        exe = folder / (script + (".exe" if os.name == "nt" else ""))
        code, out = run_command([str(exe), "--help"], timeout=SMOKE_TIMEOUT)
        if code != 0 or "usage" not in out.lower():
            bad.append(f"{script}: {tail(out, 2) or 'exit ' + str(code)}")
    return [Check(f"every installed command runs (--help), {len(scripts)} commands", not bad, "; ".join(bad))]


def run_package_tests(venv: Path, package: Path) -> Check:
    folder, python = venv_paths(venv)
    env = dict(os.environ, PATH=str(folder) + os.pathsep + os.environ.get("PATH", ""))
    code, out = run_command([str(python), "-m", "pytest", "-q", "--no-header", str(package)], env=env, timeout=TEST_TIMEOUT)
    summary = out.strip().splitlines()[-1] if out.strip() else ""
    return Check("package tests pass against the installed wheel", code == 0, summary)


def next_steps(project: dict, tag_prefix: str | None) -> list[str]:
    name, version = project["name"], project["version"]
    tag = f"{tag_prefix if tag_prefix is not None else name + '-v'}{version}"
    return [f'git tag -a {tag} -m "{name} {version}"', f"git push origin {tag}",
            f"then on GitHub: Releases > Draft a new release > choose {tag}, and paste the {name} v{version} notes from CHANGELOG.md"]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("package", nargs="?", default="spectra-kit", help="package folder under the repository root (default spectra-kit)")
    ap.add_argument("--root", default=str(ROOT), help="repository root (default: this repository)")
    ap.add_argument("--clean-venv", action="store_true", help="install the dependencies from PyPI instead of reusing this Python's")
    ap.add_argument("--allow-existing-tag", action="store_true",
                    help="do not fail when the release tag already exists (for CI runs after the release)")
    ap.add_argument("--no-isolation", action="store_true", help="build without downloading the build backend")
    ap.add_argument("--skip-tests", action="store_true", help="do not run the package tests")
    ap.add_argument("--keep", metavar="DIR", help="keep the built files in DIR instead of a temporary folder")
    ap.add_argument("--tag-prefix", help="git tag prefix (default: <package name>-v)")
    args = ap.parse_args()

    root = Path(args.root)
    package = root / args.package
    try:
        project = read_pyproject(package)
    except ValueError as err:
        print(f"error: {err}", file=sys.stderr)
        return 2

    checks = check_metadata(package, root, args.tag_prefix, args.allow_existing_tag)
    checks.append(check_sync(package))
    with tempfile.TemporaryDirectory(prefix="release-check-") as tmp:
        outdir = Path(args.keep) if args.keep else Path(tmp) / "dist"
        outdir.mkdir(parents=True, exist_ok=True)
        built, wheel, sdist = build_package(package, outdir, not args.no_isolation)
        checks += built
        if wheel and sdist:
            checks += inspect_wheel(wheel, project) + inspect_sdist(sdist, project)
            installed = install_wheel(Path(tmp) / "venv", wheel, args.clean_venv)
            checks += installed
            if all(c.ok for c in installed):
                checks += smoke_scripts(Path(tmp) / "venv", sorted(project.get("scripts", {})))
                if not args.skip_tests:
                    ready = [install_pytest(Path(tmp) / "venv")] if args.clean_venv else []
                    checks += ready
                    if all(c.ok for c in ready):
                        checks.append(run_package_tests(Path(tmp) / "venv", package))

    for check in checks:
        print(check.line())
    failed = [c for c in checks if not c.ok and not c.warning]
    print(f"\n{len(checks)} checks, {len(failed)} failed")
    if failed:
        return 1
    print(f"{project['name']} {project['version']} is ready to tag:")
    for line in next_steps(project, args.tag_prefix):
        print(f"  {line}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
