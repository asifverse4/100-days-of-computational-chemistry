"""Tests for the spectra-kit package. Run after `pip install ./spectra-kit`: pytest spectra-kit"""
import importlib.metadata
import importlib.util
import re
import shutil
import subprocess
import sys
import tomllib
from pathlib import Path

import pytest
from spectra_kit import __version__, cli
from spectra_kit.registry import BY_KEY, TOOLS, command

KIT = Path(__file__).resolve().parents[1]
ROOT = KIT.parent


def after_docstring(text: str) -> str:
    """Everything after the module docstring."""
    return text[re.search(r'"""(.*?)"""', text, re.DOTALL).end():]


def load_sync():
    spec = importlib.util.spec_from_file_location("spectra_kit_sync", KIT / "sync.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# ---------- registry and sync ----------
def test_registry_is_consistent():
    assert len({t.key for t in TOOLS}) == len(TOOLS) == 9
    assert len({t.module for t in TOOLS}) == 9 and [t.day for t in TOOLS] == list(range(11, 20))
    for tool in TOOLS:
        assert (ROOT / tool.source / "tool.py").is_file(), tool.source
        for name in tool.examples:
            assert (ROOT / tool.source / "example" / name).is_file(), name
        assert command(tool) == f"spectra-{tool.key}"


def test_package_copies_match_the_daily_tools():
    """The copies in src/ must be exactly what sync.py generates from days/ (this is the drift check)."""
    sync = load_sync()
    assert sync.stale_files(sync.expected_files()) == []


def test_sync_detects_a_stale_copy(tmp_path, monkeypatch):
    sync = load_sync()
    files = sync.expected_files()
    victim = next(p for p in files if p.name == "uvvis.py")
    original = victim.read_bytes()
    try:
        victim.write_bytes(original + b"# edited by hand\n")
        problems = sync.stale_files(files)
        assert problems == [f"out of date: {victim.relative_to(ROOT).as_posix()}"]
    finally:
        victim.write_bytes(original)
    assert sync.stale_files(files) == []


def test_sync_ignores_windows_line_endings():
    sync = load_sync()
    assert sync.same(b"a\r\nb\r\n", b"a\nb\n") and not sync.same(b"a\n", b"b\n")


def test_rendered_modules_only_change_the_docstring():
    sync = load_sync()
    for tool in TOOLS:
        original = (ROOT / tool.source / "tool.py").read_text()
        copy = sync.render_module(tool)
        assert after_docstring(copy) == after_docstring(original)  # all the code is untouched
        assert copy.startswith("# Generated from") and "python tool.py" not in copy.split('"""')[1]
        assert command(tool) in copy.split('"""')[1]


def test_version_matches_pyproject():
    project = tomllib.loads((KIT / "pyproject.toml").read_text())["project"]
    assert project["version"] == __version__
    scripts = tomllib.loads((KIT / "pyproject.toml").read_text())["project"]["scripts"]
    assert set(scripts) == {"spectra", *(command(t) for t in TOOLS)}
    for tool in TOOLS:
        assert scripts[command(tool)] == f"spectra_kit.{tool.module}:main"


def test_license_copy_matches_the_root_license():
    assert (KIT / "LICENSE").read_bytes().replace(b"\r\n", b"\n") == (ROOT / "LICENSE").read_bytes().replace(b"\r\n", b"\n")


# ---------- the dispatcher ----------
def test_list_version_and_help(capsys):
    assert cli.main(["list"]) == 0
    listing = capsys.readouterr().out
    assert all(t.key in listing and command(t) in listing for t in TOOLS)
    assert cli.main(["--version"]) == 0 and capsys.readouterr().out.strip() == f"spectra-kit {__version__}"
    assert cli.main([]) == 0 and "usage: spectra <tool>" in capsys.readouterr().out


def test_unknown_tool(capsys):
    assert cli.main(["nonsense"]) == 2
    assert "unknown tool 'nonsense'" in capsys.readouterr().err


def test_examples_copy_and_overwrite_protection(tmp_path, capsys):
    dest = tmp_path / "ex"
    assert cli.main(["examples", str(dest)]) == 0
    assert "Copied 14 example files" in capsys.readouterr().out
    for key, name in cli.example_files():
        assert (dest / key / name).is_file()
    assert cli.main(["examples", str(dest)]) == 2  # refuses to overwrite
    assert "already exist" in capsys.readouterr().err
    assert cli.main(["examples", str(dest), "--force"]) == 0
    assert cli.main(["examples", "a", "b"]) == 2


def test_every_tool_help_shows_the_installed_command_name(capsys):
    for tool in TOOLS:
        with pytest.raises(SystemExit) as exit_info:
            cli.main([tool.key, "--help"])
        assert exit_info.value.code == 0
        out = capsys.readouterr().out
        assert f"spectra {tool.key}" in out  # argparse's usage line
        assert command(tool) in out and "python tool.py" not in out  # the docstring was rewritten


def test_dispatcher_restores_sys_argv():
    before = list(sys.argv)
    cli.main(["list"])
    assert sys.argv == before
    cli.run_tool("stokes", ["345", "450"])
    assert sys.argv == before


# ---------- every tool, on its packaged examples ----------
@pytest.fixture(scope="module")
def examples(tmp_path_factory):
    dest = tmp_path_factory.mktemp("examples") / "spectra-examples"
    assert cli.main(["examples", str(dest)]) == 0
    return dest


@pytest.mark.parametrize("key, args, expected", [
    ("uvvis", ["{ex}/uvvis/synthetic_cqd.csv", "-o", "{out}/uvvis.png", "--dpi", "60"], "[ok] 2 spectra"),
    ("tauc", ["{ex}/tauc/synthetic_tauc.csv", "-o", "{out}/tauc.png", "--dpi", "60"], "[ok] 2 spectra"),
    ("qy", ["--ref", "{ex}/qy/synthetic_reference.csv", "--sample", "{ex}/qy/synthetic_sample.csv", "--standard",
            "fluorescein", "--sample-solvent", "ethanol"], "Relative quantum yield (gradient method)"),
    ("pl", ["{ex}/pl/synthetic_pl.csv", "-o", "{out}/pl.png", "--mask-scatter", "15", "--dpi", "60"], "[ok] 8 spectra"),
    ("stokes", ["--abs", "{ex}/stokes/synthetic_absorption.csv", "--em", "{ex}/stokes/synthetic_emission.csv",
                "--abs-peak", "longest"], "105.4"),
    ("ftir", ["{ex}/ftir/synthetic_ftir.csv"], "8 peaks"),
    ("xrd", ["{ex}/xrd/synthetic_xrd.csv", "--inst-fwhm", "0.08", "--williamson-hall"], "Williamson-Hall (6 peaks)"),
    ("tem", ["{ex}/tem/synthetic_sizes.csv", "--column", "diameter_nm"], "n = 250 particles"),
    ("baseline", ["{ex}/baseline/synthetic_spectrum.csv", "--baseline", "asls", "--baseline-lambda", "1e7", "--truth",
                  "{ex}/baseline/synthetic_baseline_truth.csv"], "against the true baseline"),
])
def test_each_tool_runs_on_its_packaged_examples(key, args, expected, examples, tmp_path, capsys):
    argv = [a.format(ex=examples, out=tmp_path) for a in args]
    assert cli.main([key, *argv]) == 0
    assert expected in capsys.readouterr().out


def test_every_registered_tool_is_covered_by_the_example_test():
    covered = {"uvvis", "tauc", "qy", "pl", "stokes", "ftir", "xrd", "tem", "baseline"}
    assert covered == {t.key for t in TOOLS}


# ---------- the installed package ----------
def installed_scripts():
    try:
        dist = importlib.metadata.distribution("spectra-kit")
    except importlib.metadata.PackageNotFoundError:
        pytest.skip("spectra-kit is not installed (run: pip install ./spectra-kit)")
    return {ep.name: ep.value for ep in dist.entry_points if ep.group == "console_scripts"}


def test_installed_entry_points():
    scripts = installed_scripts()
    assert set(scripts) == {"spectra", *(command(t) for t in TOOLS)}
    assert importlib.metadata.version("spectra-kit") == __version__


def test_installed_package_ships_the_examples():
    installed_scripts()
    for key, name in cli.example_files():
        from importlib import resources
        assert (resources.files("spectra_kit") / "examples" / key / name).is_file()


@pytest.mark.parametrize("name", ["spectra", *(command(t) for t in TOOLS)])
def test_console_scripts_run(name):
    installed_scripts()
    exe = shutil.which(name)
    if exe is None:
        pytest.skip(f"{name} is not on the PATH")
    result = subprocess.run([exe, "--help"], check=False, capture_output=True, text=True, timeout=120)
    assert result.returncode == 0 and re.search(r"usage:", result.stdout)


def test_python_dash_m_works():
    result = subprocess.run([sys.executable, "-m", "spectra_kit", "--version"], check=False, capture_output=True,
                            text=True, timeout=120)
    assert result.returncode == 0 and result.stdout.strip() == f"spectra-kit {__version__}"


def test_baseline_to_ftir_chain(examples, tmp_path, capsys):
    """The output of one tool is readable by another: baseline-correct a spectrum, then pick its peaks."""
    corrected = tmp_path / "corrected.csv"
    assert cli.main(["baseline", str(examples / "baseline" / "synthetic_spectrum.csv"), "--baseline", "asls",
                     "--baseline-lambda", "1e7", "--csv", str(corrected), "--corrected-only"]) == 0
    capsys.readouterr()
    assert cli.main(["ftir", str(corrected), "--signal", "absorbance"]) == 0
    assert "6 peaks" in capsys.readouterr().out
    assert BY_KEY["ftir"].day == 16
