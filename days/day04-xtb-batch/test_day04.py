import subprocess
from pathlib import Path

import pytest
import tool

# Excerpt in the format xtb prints at the end of an optimization
CONVERGED = """
          | TOTAL ENERGY               -8.123456789012 Eh   |
          | GRADIENT NORM               0.000812345678 Eh/a0|
          | HOMO-LUMO GAP              12.345678901234 eV   |
 *** GEOMETRY OPTIMIZATION CONVERGED AFTER 7 ITERATIONS ***
         :::::::::::::::::::::::::::::::::::::::::::::::::::::
         ::                     SUMMARY                     ::
         :::::::::::::::::::::::::::::::::::::::::::::::::::::
         :: total energy              -8.130000000000 Eh    ::
         :: gradient norm              0.000400000000 Eh/a0 ::
         :: HOMO-LUMO gap             12.400000000000 eV    ::
"""
NOT_CONVERGED = """
          | TOTAL ENERGY               -5.000000000000 Eh   |
 *** FAILED TO CONVERGE GEOMETRY OPTIMIZATION IN 200 ITERATIONS ***
"""
CRASHED = "abnormal termination of xtb\n"
XYZ = "3\nwater | energy = 1.0 kcal/mol\nO 0.0 0.0 0.0\nH 0.0 0.0 1.0\nH 0.0 1.0 0.0\n"


def test_parse_converged_uses_last_values():
    r = tool.parse_xtb_output(CONVERGED)
    assert r["status"] == "converged"
    assert r["energy_eh"] == pytest.approx(-8.13)
    assert r["gap_ev"] == pytest.approx(12.4)
    assert r["grad_norm"] == pytest.approx(0.0004)
    assert r["iterations"] == 7


def test_parse_not_converged():
    r = tool.parse_xtb_output(NOT_CONVERGED)
    assert r["status"] == "not_converged"
    assert r["gap_ev"] is None


def test_parse_failed():
    r = tool.parse_xtb_output(CRASHED)
    assert r["status"] == "failed"
    assert r["energy_eh"] is None


def test_clean_xyz_replaces_comment():
    out = tool.clean_xyz(XYZ, "water")
    assert out.splitlines()[:2] == ["3", "water"]
    assert len(out.splitlines()) == 5


def test_clean_xyz_rejects_bad_files():
    with pytest.raises(ValueError):
        tool.clean_xyz("hello\nworld\n", "x")
    with pytest.raises(ValueError):
        tool.clean_xyz("5\ncomment\nH 0 0 0\n", "x")


def test_collect_entries_from_folder_with_duplicates_and_bad_file(tmp_path):
    (tmp_path / "a.xyz").write_text(XYZ)
    (tmp_path / "b b.xyz").write_text(XYZ)
    (tmp_path / "bad.xyz").write_text("nope")
    entries = tool.collect_entries(str(tmp_path), charge=-1, uhf=0)
    assert [e["name"] for e in entries] == ["a", "b_b", "bad"]
    assert entries[0]["charge"] == -1
    assert entries[2]["error"] is not None


def test_collect_entries_errors(tmp_path):
    with pytest.raises(ValueError):
        tool.collect_entries(str(tmp_path / "missing"))
    with pytest.raises(ValueError):
        tool.collect_entries(str(tmp_path))  # empty folder


def test_unique_names():
    seen = set()
    assert [tool.unique("x", seen) for _ in range(3)] == ["x", "x_2", "x_3"]


def test_run_xtb_builds_command_and_parses(tmp_path, monkeypatch):
    calls = {}

    def fake_run(cmd, cwd, **kw):
        calls["cmd"] = cmd
        (Path(cwd) / "xtbopt.xyz").write_text(XYZ)
        return subprocess.CompletedProcess(cmd, 0, stdout=CONVERGED, stderr="")

    monkeypatch.setattr(tool.subprocess, "run", fake_run)
    entry = {"name": "water", "xyz": tool.clean_xyz(XYZ, "water"), "charge": 0, "uhf": 0, "error": None}
    row = tool.run_xtb(entry, tmp_path, solvent="water", opt_level="tight")
    assert row["status"] == "converged"
    assert row["energy_eh"] == pytest.approx(-8.13)
    assert row["opt_xyz"].endswith("xtbopt.xyz")
    assert "--alpb" in calls["cmd"] and "water" in calls["cmd"]
    assert calls["cmd"][calls["cmd"].index("--opt") + 1] == "tight"
    assert (tmp_path / "water" / "xtb.out").is_file()


def test_run_xtb_timeout(tmp_path, monkeypatch):
    def fake_run(cmd, **kw):
        raise subprocess.TimeoutExpired(cmd, 1)

    monkeypatch.setattr(tool.subprocess, "run", fake_run)
    entry = {"name": "m", "xyz": tool.clean_xyz(XYZ, "m"), "charge": 0, "uhf": 0, "error": None}
    assert tool.run_xtb(entry, tmp_path, timeout=1)["status"] == "timeout"


def test_main_without_xtb_returns_2(monkeypatch, capsys):
    monkeypatch.setattr(tool.shutil, "which", lambda _: None)
    monkeypatch.setattr("sys.argv", ["tool.py", "whatever"])
    assert tool.main() == 2
    assert "xtb executable not found" in capsys.readouterr().err


def test_smiles_input_sets_charge_and_atoms(tmp_path):
    pytest.importorskip("rdkit")
    smi = tmp_path / "in.smi"
    smi.write_text("# c\nCC(=O)[O-] acetate\nCCO ethanol\nnot_smiles bad\n")
    entries = tool.collect_entries(str(smi))
    assert entries[0]["charge"] == -1
    assert entries[1]["xyz"].startswith("9\n")  # ethanol: 9 atoms with H
    assert entries[2]["error"] is not None
