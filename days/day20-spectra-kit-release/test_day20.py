import io
import os
import subprocess
import tarfile
import zipfile
from pathlib import Path

import pytest
import tool

REAL_ROOT = tool.ROOT


# ---------- a fake package and a fake command runner ----------
PYPROJECT = """[project]
name = "demo-kit"
version = "1.2.3"
readme = "README.md"
license-files = ["LICENSE"]

[project.scripts]
demo = "demo_kit.cli:main"
demo-two = "demo_kit.two:main"

[tool.setuptools.package-data]
demo_kit = ["data/*.csv"]
"""


def make_repo(tmp_path, *, version="1.2.3", module_version="1.2.3", changelog=True, readme=True, licence=True,
              scripts=True, package_data=True):
    root = tmp_path / "repo"
    pkg = root / "demo-kit"
    (pkg / "src" / "demo_kit").mkdir(parents=True)
    text = PYPROJECT.replace('version = "1.2.3"', f'version = "{version}"')
    if not scripts:
        text = text.replace("[project.scripts]\ndemo = \"demo_kit.cli:main\"\ndemo-two = \"demo_kit.two:main\"\n", "")
    if not package_data:
        text = text.split("[tool.setuptools")[0]
    (pkg / "pyproject.toml").write_text(text)
    if module_version is not None:
        (pkg / "src" / "demo_kit" / "__init__.py").write_text(f'__version__ = "{module_version}"\n')
    if readme:
        (pkg / "README.md").write_text("# demo\n")
    if licence:
        (pkg / "LICENSE").write_text("MIT\n")
    if changelog:
        (root / "CHANGELOG.md").write_text("### demo-kit v1.2.3 release notes\n- first release\n")
    return root, pkg


SCRIPTS = {"demo": "demo_kit.cli:main", "demo-two": "demo_kit.two:main"}


def make_wheel(path, *, name="demo-kit", version="1.2.3", licence=True, scripts=("demo", "demo-two"), tests=False,
               data=True, metadata=True):
    dist = f"{name.replace('-', '_')}-{version}.dist-info"
    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr("demo_kit/__init__.py", "")
        if data:
            zf.writestr("demo_kit/data/example.csv", "a,b\n1,2\n")
        if tests:
            zf.writestr("tests/test_x.py", "")
        if metadata:
            zf.writestr(f"{dist}/METADATA", f"Name: {name}\nVersion: {version}\n")
        if licence:
            zf.writestr(f"{dist}/licenses/LICENSE", "MIT")
        lines = "\n".join(f"{s} = {SCRIPTS.get(s, 'x:y')}" for s in scripts)
        zf.writestr(f"{dist}/entry_points.txt", f"[console_scripts]\n{lines}\n\n[other]\nx = y:z\n")
    return path


def make_sdist(path, *, name="demo-kit", version="1.2.3", missing=(), tests=False):
    base = f"{name.replace('-', '_')}-{version}"
    with tarfile.open(path, "w:gz") as tf:
        for member in ("pyproject.toml", "LICENSE", "README.md", "tests/test_x.py" if tests else None):
            if member and member not in missing:
                info = tarfile.TarInfo(f"{base}/{member}")
                info.size = 1
                tf.addfile(info, io.BytesIO(b"x"))
    return path


class FakeRunner:
    """Stands in for tool.run_command: records calls and builds a believable wheel and sdist."""

    def __init__(self, project, fail=None):
        self.calls, self.fail, self.project = [], fail or {}, project

    def __call__(self, cmd, cwd=None, env=None, timeout=600):
        cmd = [str(c) for c in cmd]
        self.calls.append({"cmd": cmd, "env": env})
        joined = " ".join(cmd)
        for needle, result in self.fail.items():
            if needle in joined:
                return result
        if "-m build" in joined:
            outdir = Path(cmd[cmd.index("--outdir") + 1])
            make_wheel(outdir / "demo_kit-1.2.3-py3-none-any.whl")
            make_sdist(outdir / "demo_kit-1.2.3.tar.gz")
        if cmd[-1] == "--help":
            return 0, f"usage: {Path(cmd[0]).name} [-h]"
        if " -m pytest " in f" {joined} ":
            return 0, "36 passed in 12.3s"
        return 0, ""

    def commands(self):
        return [" ".join(c["cmd"]) for c in self.calls]


# ---------- the pyproject ----------
def test_read_pyproject(tmp_path):
    _root, pkg = make_repo(tmp_path)
    project = tool.read_pyproject(pkg)
    assert project["name"] == "demo-kit" and project["scripts"]["demo"] == "demo_kit.cli:main"
    assert project["_declares_package_data"] is True
    assert tool.read_pyproject(make_repo(tmp_path / "b", package_data=False)[1])["_declares_package_data"] is False


def test_read_pyproject_errors(tmp_path):
    with pytest.raises(ValueError, match="not found"):
        tool.read_pyproject(tmp_path)
    (tmp_path / "pyproject.toml").write_text("this is [not toml")
    with pytest.raises(ValueError, match="not valid TOML"):
        tool.read_pyproject(tmp_path)
    (tmp_path / "pyproject.toml").write_text('[project]\nname = "x"\n')
    with pytest.raises(ValueError, match="no project.version"):
        tool.read_pyproject(tmp_path)


def test_module_version(tmp_path):
    _root, pkg = make_repo(tmp_path)
    assert tool.module_version(pkg, "demo-kit") == "1.2.3"
    (pkg / "src" / "demo_kit" / "__init__.py").write_text("__version__ = '9.9.9'\n")
    assert tool.module_version(pkg, "demo-kit") == "9.9.9"
    (pkg / "src" / "demo_kit" / "__init__.py").write_text("# no version\n")
    assert tool.module_version(pkg, "demo-kit") is None
    assert tool.module_version(pkg, "other-name") is None


# ---------- metadata checks ----------
@pytest.fixture
def no_git(monkeypatch):
    monkeypatch.setattr(tool, "tag_exists", lambda root, tag: False)


def failed(checks):
    return [c.name for c in checks if not c.ok and not c.warning]


def test_a_correct_package_passes(tmp_path, no_git):
    root, pkg = make_repo(tmp_path)
    checks = tool.check_metadata(pkg, root)
    assert failed(checks) == [] and len(checks) == 7
    assert any("git tag demo-kit-v1.2.3 is free" in c.name for c in checks)


@pytest.mark.parametrize("kwargs, expected", [
    ({"version": "1.0", "module_version": "1.0"}, "version format"),
    ({"module_version": "1.2.4"}, "__version__ matches pyproject"),
    ({"module_version": None}, "__version__ matches pyproject"),
    ({"readme": False}, "README is present"),
    ({"licence": False}, "license file is present"),
    ({"scripts": False}, "console scripts are declared"),
    ({"changelog": False}, "CHANGELOG mentions this release"),
])
def test_each_metadata_problem_is_reported(tmp_path, no_git, kwargs, expected):
    root, pkg = make_repo(tmp_path, **kwargs)
    assert expected in failed(tool.check_metadata(pkg, root))


def test_changelog_must_name_the_version(tmp_path, no_git):
    root, pkg = make_repo(tmp_path)
    (root / "CHANGELOG.md").write_text("### demo-kit v1.2.0 release notes\n")
    assert "CHANGELOG mentions this release" in failed(tool.check_metadata(pkg, root))


def test_existing_tag_fails_unless_allowed(tmp_path, monkeypatch):
    root, pkg = make_repo(tmp_path)
    monkeypatch.setattr(tool, "tag_exists", lambda r, t: True)
    assert any("is free" in n for n in failed(tool.check_metadata(pkg, root)))
    allowed = tool.check_metadata(pkg, root, allow_existing_tag=True)
    assert failed(allowed) == []
    tag_check = next(c for c in allowed if "is free" in c.name)
    assert tag_check.warning and not tag_check.ok and "allowed" in tag_check.detail and tag_check.line().startswith("[warn]")


def test_unknown_git_state_is_only_a_warning(tmp_path, monkeypatch):
    root, pkg = make_repo(tmp_path)
    monkeypatch.setattr(tool, "tag_exists", lambda r, t: None)
    checks = tool.check_metadata(pkg, root)
    assert failed(checks) == [] and any(c.warning and "not checked" in c.detail for c in checks)


def test_custom_tag_prefix(tmp_path, no_git):
    root, pkg = make_repo(tmp_path)
    assert any(c.name == "git tag v1.2.3 is free" for c in tool.check_metadata(pkg, root, tag_prefix="v"))


def test_tag_exists_with_real_git(tmp_path):
    git = ["git", "-c", "user.name=t", "-c", "user.email=t@t"]
    subprocess.run([*git, "init", "-q"], cwd=tmp_path, check=True)
    subprocess.run([*git, "commit", "-q", "--allow-empty", "-m", "x"], cwd=tmp_path, check=True)
    subprocess.run([*git, "tag", "demo-kit-v1.0.0"], cwd=tmp_path, check=True)
    assert tool.tag_exists(tmp_path, "demo-kit-v1.0.0") is True
    assert tool.tag_exists(tmp_path, "demo-kit-v2.0.0") is False
    assert tool.tag_exists(tmp_path / "missing-folder", "x") is None  # git cannot run there


# ---------- the sync check ----------
def test_check_sync(tmp_path):
    _root, pkg = make_repo(tmp_path)
    none = tool.check_sync(pkg)
    assert none.ok and none.warning and "nothing to check" in none.detail
    (pkg / "sync.py").write_text("import sys\nprint('all in step')\nsys.exit(0)\n")
    assert tool.check_sync(pkg).ok
    (pkg / "sync.py").write_text("import sys\nprint('[error] out of date: uvvis.py')\nsys.exit(1)\n")
    bad = tool.check_sync(pkg)
    assert not bad.ok and "out of date: uvvis.py" in bad.detail


# ---------- inspecting the built files ----------
PROJECT = {"name": "demo-kit", "version": "1.2.3", "scripts": dict(SCRIPTS), "_declares_package_data": True}


def test_parse_entry_points():
    text = "[console_scripts]\na = x:y\nb = z:w\n\n[gui_scripts]\nc = q:r\n"
    assert tool.parse_entry_points(text) == {"a": "x:y", "b": "z:w"}


def test_a_good_wheel(tmp_path):
    assert failed(tool.inspect_wheel(make_wheel(tmp_path / "w.whl"), PROJECT)) == []


@pytest.mark.parametrize("kwargs, expected", [
    ({"metadata": False}, "wheel metadata has the right name and version"),
    ({"version": "9.9.9"}, "wheel metadata has the right name and version"),  # a wheel of a different version
    ({"licence": False}, "wheel carries the license"),
    ({"scripts": ("demo",)}, "wheel entry points match pyproject"),
    ({"scripts": ("demo", "demo-two", "extra")}, "wheel entry points match pyproject"),
    ({"tests": True}, "wheel has no tests or build scripts"),
    ({"data": False}, "wheel includes its data files"),
])
def test_wheel_problems(tmp_path, kwargs, expected):
    wheel = make_wheel(tmp_path / "w.whl", **kwargs)
    assert expected in failed(tool.inspect_wheel(wheel, PROJECT))


def test_wheel_without_declared_data_does_not_need_csv_files(tmp_path):
    project = dict(PROJECT, _declares_package_data=False)
    assert "wheel includes its data files" not in failed(tool.inspect_wheel(make_wheel(tmp_path / "w.whl", data=False), project))


def test_sdist_checks(tmp_path):
    assert failed(tool.inspect_sdist(make_sdist(tmp_path / "a.tar.gz"), PROJECT)) == []
    missing = tool.inspect_sdist(make_sdist(tmp_path / "b.tar.gz", missing=("LICENSE",)), PROJECT)
    assert "source distribution has pyproject, license and README" in failed(missing)
    with_tests = tool.inspect_sdist(make_sdist(tmp_path / "c.tar.gz", tests=True), PROJECT)
    assert failed(with_tests) == [] and any(c.warning and "1 test files" in c.detail for c in with_tests)


# ---------- building, installing and running ----------
def test_build_package(tmp_path, monkeypatch):
    runner = FakeRunner(PROJECT)
    monkeypatch.setattr(tool, "run_command", runner)
    out = tmp_path / "dist"
    out.mkdir()
    checks, wheel, sdist = tool.build_package(tmp_path / "pkg", out)
    assert checks[0].ok and wheel.name.endswith(".whl") and sdist.name.endswith(".tar.gz")
    assert "--no-isolation" not in runner.commands()[0]
    out2 = tmp_path / "dist2"
    out2.mkdir()
    tool.build_package(tmp_path / "pkg", out2, isolation=False)
    assert "--no-isolation" in runner.commands()[1]


def test_build_failure_and_missing_outputs(tmp_path, monkeypatch):
    out = tmp_path / "dist"
    out.mkdir()
    monkeypatch.setattr(tool, "run_command", FakeRunner(PROJECT, fail={"-m build": (1, "line1\nerror: backend failed")}))
    checks, wheel, sdist = tool.build_package(tmp_path, out)
    assert not checks[0].ok and "backend failed" in checks[0].detail and wheel is None and sdist is None
    monkeypatch.setattr(tool, "run_command", lambda *a, **k: (0, ""))  # "succeeds" but writes nothing
    assert not tool.build_package(tmp_path, out)[0][0].ok


def test_venv_paths():
    folder, python = tool.venv_paths(Path("v"), windows=False)
    assert folder == Path("v/bin") and python == Path("v/bin/python")
    folder, python = tool.venv_paths(Path("v"), windows=True)
    assert folder == Path("v/Scripts") and python == Path("v/Scripts/python.exe")


def test_install_wheel_modes_and_failures(tmp_path, monkeypatch):
    runner = FakeRunner(PROJECT)
    monkeypatch.setattr(tool, "run_command", runner)
    assert all(c.ok for c in tool.install_wheel(tmp_path / "v", tmp_path / "w.whl", clean=False))
    assert "--system-site-packages" in runner.commands()[0] and "pip install --quiet" in runner.commands()[1]
    runner2 = FakeRunner(PROJECT)
    monkeypatch.setattr(tool, "run_command", runner2)
    checks = tool.install_wheel(tmp_path / "v", tmp_path / "w.whl", clean=True)
    assert "--system-site-packages" not in runner2.commands()[0] and "PyPI" in checks[0].detail
    monkeypatch.setattr(tool, "run_command", FakeRunner(PROJECT, fail={"-m venv": (1, "no venv module")}))
    only = tool.install_wheel(tmp_path / "v", tmp_path / "w.whl", clean=False)
    assert len(only) == 1 and not only[0].ok
    monkeypatch.setattr(tool, "run_command", FakeRunner(PROJECT, fail={"pip install": (1, "no matching distribution")}))
    assert not tool.install_wheel(tmp_path / "v", tmp_path / "w.whl", clean=False)[1].ok


def test_smoke_scripts(tmp_path, monkeypatch):
    monkeypatch.setattr(tool, "run_command", FakeRunner(PROJECT))
    (check,) = tool.smoke_scripts(tmp_path, ["demo", "demo-two"])
    assert check.ok and "2 commands" in check.name
    monkeypatch.setattr(tool, "run_command", FakeRunner(PROJECT, fail={"demo-two": (1, "ImportError: no module named x")}))
    (check,) = tool.smoke_scripts(tmp_path, ["demo", "demo-two"])
    assert not check.ok and "demo-two" in check.detail and "ImportError" in check.detail
    monkeypatch.setattr(tool, "run_command", lambda *a, **k: (0, "no help text here"))
    assert not tool.smoke_scripts(tmp_path, ["demo"])[0].ok  # exit 0 but not a usage message


def test_package_tests_run_with_the_venv_on_the_path(tmp_path, monkeypatch):
    runner = FakeRunner(PROJECT)
    monkeypatch.setattr(tool, "run_command", runner)
    check = tool.run_package_tests(tmp_path / "venv", tmp_path / "pkg")
    assert check.ok and "36 passed" in check.detail
    assert str(tool.venv_paths(tmp_path / "venv")[0]) in runner.calls[0]["env"]["PATH"].split(os.pathsep)[0]
    monkeypatch.setattr(tool, "run_command", FakeRunner(PROJECT, fail={" -m pytest ": (1, "1 failed, 35 passed")}))
    assert not tool.run_package_tests(tmp_path / "venv", tmp_path / "pkg").ok


def test_next_steps():
    steps = tool.next_steps({"name": "demo-kit", "version": "1.2.3"}, None)
    assert steps[0] == 'git tag -a demo-kit-v1.2.3 -m "demo-kit 1.2.3"' and steps[1] == "git push origin demo-kit-v1.2.3"
    assert tool.next_steps({"name": "demo-kit", "version": "1.2.3"}, "v")[1] == "git push origin v1.2.3"


# ---------- the whole command, with a fake runner ----------
def run_main(monkeypatch, root, *argv, runner=None):
    runner = runner or FakeRunner(PROJECT)
    monkeypatch.setattr(tool, "run_command", runner)
    monkeypatch.setattr(tool, "tag_exists", lambda r, t: False)
    monkeypatch.setattr("sys.argv", ["tool.py", "demo-kit", "--root", str(root), *argv])
    return tool.main(), runner


def test_a_clean_release_passes_end_to_end(tmp_path, monkeypatch, capsys):
    root, _ = make_repo(tmp_path)
    code, runner = run_main(monkeypatch, root)
    out = capsys.readouterr().out
    assert code == 0 and "0 failed" in out and "demo-kit 1.2.3 is ready to tag" in out
    assert 'git tag -a demo-kit-v1.2.3 -m "demo-kit 1.2.3"' in out
    assert any("pytest" in c for c in runner.commands()) and not any("install --quiet pytest" in c for c in runner.commands())


def test_options_change_what_is_run(tmp_path, monkeypatch):
    root, _ = make_repo(tmp_path)
    code, runner = run_main(monkeypatch, root, "--skip-tests", "--no-isolation")
    assert code == 0 and not any("-m pytest" in c for c in runner.commands())
    assert any("--no-isolation" in c for c in runner.commands())
    code, runner = run_main(monkeypatch, root, "--clean-venv")
    commands = runner.commands()
    assert code == 0 and any("install --quiet pytest" in c for c in commands)
    assert not any("--system-site-packages" in c for c in commands)


def test_keep_puts_the_built_files_in_a_folder(tmp_path, monkeypatch):
    root, _ = make_repo(tmp_path)
    keep = tmp_path / "dist-kept"
    code, _ = run_main(monkeypatch, root, "--keep", str(keep))
    assert code == 0 and len(list(keep.glob("*.whl"))) == 1 and len(list(keep.glob("*.tar.gz"))) == 1


@pytest.mark.parametrize("needle, expected_line", [
    ("-m build", "build the wheel"),
    ("pip install --quiet", "install the wheel"),
    ("--help", "every installed command runs"),
    (" -m pytest ", "package tests pass"),
])
def test_a_failing_step_fails_the_release(tmp_path, monkeypatch, capsys, needle, expected_line):
    root, _ = make_repo(tmp_path)
    runner = FakeRunner(PROJECT, fail={needle: (1, "something broke")})
    code, _ = run_main(monkeypatch, root, runner=runner)
    out = capsys.readouterr().out
    assert code == 1 and "[FAIL] " + expected_line in out and "ready to tag" not in out


def test_metadata_problem_fails_but_the_rest_still_runs(tmp_path, monkeypatch, capsys):
    root, _ = make_repo(tmp_path, changelog=False)
    code, _ = run_main(monkeypatch, root)
    out = capsys.readouterr().out
    assert code == 1 and "[FAIL] CHANGELOG mentions this release" in out and "[PASS] install the wheel" in out


def test_a_failed_build_skips_the_install_steps(tmp_path, monkeypatch, capsys):
    root, _ = make_repo(tmp_path)
    code, runner = run_main(monkeypatch, root, runner=FakeRunner(PROJECT, fail={"-m build": (1, "no backend")}))
    assert code == 1 and not any("-m venv" in c for c in runner.commands())


def test_usage_errors(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr("sys.argv", ["tool.py", "no-such-package", "--root", str(tmp_path)])
    assert tool.main() == 2 and "not found" in capsys.readouterr().err
    (tmp_path / "bad").mkdir()
    (tmp_path / "bad" / "pyproject.toml").write_text("[broken")
    monkeypatch.setattr("sys.argv", ["tool.py", "bad", "--root", str(tmp_path)])
    assert tool.main() == 2 and "not valid TOML" in capsys.readouterr().err


def test_check_line_formats():
    assert tool.Check("a", True).line() == "[PASS] a"
    assert tool.Check("a", False, "why").line() == "[FAIL] a: why"
    assert tool.Check("a", False, "careful", warning=True).line() == "[warn] a: careful"


def test_run_command_reports_failures_instead_of_raising():
    assert tool.run_command(["definitely-not-a-program-xyz"])[0] == 127
    code, out = tool.run_command([__import__("sys").executable, "-c", "print('hi')"])
    assert code == 0 and out == "hi"
    code, out = tool.run_command([__import__("sys").executable, "-c", "import time; time.sleep(5)"], timeout=1)
    assert code == 124 and "timed out" in out


# ---------- the real package (no build) ----------
def test_the_real_spectra_kit_passes_the_static_checks():
    package = REAL_ROOT / "spectra-kit"
    checks = tool.check_metadata(package, REAL_ROOT, allow_existing_tag=True)
    assert failed(checks) == [], [c.line() for c in checks]
    assert failed([tool.check_sync(package)]) == []


@pytest.mark.skipif(os.environ.get("RELEASE_CHECK_FULL") != "1", reason="builds and installs the real package (set RELEASE_CHECK_FULL=1)")
def test_full_release_check_of_the_real_package(monkeypatch):
    monkeypatch.setattr("sys.argv", ["tool.py", "--skip-tests", "--allow-existing-tag"])
    assert tool.main() == 0
