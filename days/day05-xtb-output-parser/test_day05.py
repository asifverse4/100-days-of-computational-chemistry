from pathlib import Path

import pytest
import tool

EXAMPLE = Path(__file__).parent / "example" / "synthetic"

# Hand-written excerpt in xtb's output format (values are made up for testing)
SAMPLE = """
 * xtb version 6.7.1 (demo) compiled by 'demo' on 2026-01-01
          :  # atoms                       3                      :
          :  Hamiltonian                  GFN2-xTB                :
   #   Z          covCN         q      C6AA      α(0)
   1   8 O        2.000    -0.500     10.000     5.000
   2   1 H        1.000     0.250      3.000     2.000
   3   1 H        1.000     0.250      3.000     2.000

    4        2.0000           -0.5000000 (HOMO)    -13.6057
    5                         -0.0500000 (LUMO)     -1.3606
          | TOTAL ENERGY               -4.900000000000 Eh   |
          | HOMO-LUMO GAP               9.000000000000 eV   |
 *** GEOMETRY OPTIMIZATION CONVERGED AFTER 4 ITERATIONS ***
         :: total energy              -5.000000000000 Eh    ::
         :: gradient norm              0.000500000000 Eh/a0 ::
         :: HOMO-LUMO gap             12.245100000000 eV    ::
         :: SCC energy                -5.100000000000 Eh    ::
         :: -> dispersion             -0.000100000000 Eh    ::
         :: repulsion energy           0.100000000000 Eh    ::
         :: total charge               0.000000000000 e     ::
normal termination of xtb
"""


def test_parse_text_uses_last_values():
    r = tool.parse_xtb_text(SAMPLE)
    assert r["status"] == "converged"
    assert r["total_energy_eh"] == pytest.approx(-5.0)
    assert r["gap_ev"] == pytest.approx(12.2451)
    assert r["gradient_norm"] == pytest.approx(0.0005)
    assert r["homo_ev"] == pytest.approx(-13.6057)
    assert r["lumo_ev"] == pytest.approx(-1.3606)
    assert r["scc_energy_eh"] == pytest.approx(-5.1)
    assert r["dispersion_eh"] == pytest.approx(-0.0001)
    assert r["total_charge"] == 0.0
    assert (r["iterations"], r["n_atoms"]) == (4, 3)
    assert (r["method"], r["xtb_version"]) == ("GFN2-xTB", "6.7.1")


def test_gap_falls_back_to_orbital_energies():
    text = "xtb version 6.7.1\n 4  2.0  -0.5 (HOMO)  -13.6000\n 5  -0.05 (LUMO)  -1.6000\n"
    assert tool.parse_xtb_text(text)["gap_ev"] == pytest.approx(12.0)


def test_status_singlepoint_failed_and_not_converged():
    sp = "xtb version 6\n         :: total energy   -1.000000 Eh    ::\n"
    assert tool.parse_xtb_text(sp)["status"] == "singlepoint"
    assert tool.parse_xtb_text("xtb version 6\nabnormal termination of xtb\n")["status"] == "failed"
    nc = sp + " *** FAILED TO CONVERGE GEOMETRY OPTIMIZATION IN 200 ITERATIONS ***\n"
    assert tool.parse_xtb_text(nc)["status"] == "not_converged"


def test_parse_charge_table():
    table = tool.parse_charge_table(SAMPLE)
    assert [(r["element"], r["charge"]) for r in table] == [("O", -0.5), ("H", 0.25), ("H", 0.25)]


def test_is_xtb_output():
    assert tool.is_xtb_output(SAMPLE)
    assert not tool.is_xtb_output("just some log text")


def test_charges_file_preferred_over_table(tmp_path):
    (tmp_path / "xtb.out").write_text(SAMPLE)
    (tmp_path / "charges").write_text("-0.6\n0.3\n0.3\n")
    (tmp_path / "xtbopt.xyz").write_text("3\nw\nO 0 0 0\nH 0 1 0\nH 1 0 0\n")
    row, atoms = tool.parse_run(tmp_path / "xtb.out")
    assert row["charge_source"].startswith("charges_file")
    assert [a["charge"] for a in atoms] == [-0.6, 0.3, 0.3]
    assert [a["element"] for a in atoms] == ["O", "H", "H"]


def test_table_used_when_no_charges_file(tmp_path):
    (tmp_path / "xtb.out").write_text(SAMPLE)
    row, atoms = tool.parse_run(tmp_path / "xtb.out")
    assert row["charge_source"].startswith("output_table")
    assert len(atoms) == 3


def test_run_name_from_folder_or_file(tmp_path):
    assert tool.run_name(tmp_path / "ethanol" / "xtb.out") == "ethanol"
    assert tool.run_name(tmp_path / "conf1.out") == "conf1"


def test_find_outputs_skips_non_xtb_files_in_folders(tmp_path):
    (tmp_path / "a").mkdir()
    (tmp_path / "a" / "xtb.out").write_text(SAMPLE)
    (tmp_path / "notes.log").write_text("not xtb")
    paths = [p.name for p, _ in tool.find_outputs([str(tmp_path)])]
    assert sorted(paths) == ["notes.log", "xtb.out"]
    with pytest.raises(tool.NotXtbOutput):
        tool.parse_run(tmp_path / "notes.log")
    with pytest.raises(ValueError):
        tool.find_outputs([str(tmp_path / "missing")])


def test_relative_energies():
    rows = [{"total_energy_eh": -5.0, "n_atoms": 9}, {"total_energy_eh": -4.999, "n_atoms": 9},
            {"total_energy_eh": None, "n_atoms": 9}]
    tool.add_relative_energies(rows)
    assert rows[0]["rel_energy_kcal"] == 0.0
    assert rows[1]["rel_energy_kcal"] == pytest.approx(0.001 * tool.EH_TO_KCAL, abs=1e-3)
    assert rows[2]["rel_energy_kcal"] is None


def test_cli_end_to_end_on_example(tmp_path, monkeypatch):
    out, ch = tmp_path / "s.csv", tmp_path / "c.csv"
    monkeypatch.setattr("sys.argv", ["tool.py", str(EXAMPLE), "-o", str(out), "--charges-csv", str(ch), "--relative"])
    assert tool.main() == 0
    lines = out.read_text().splitlines()
    assert lines[0].split(",")[:7] == ["name", "status", "method", "xtb_version", "n_atoms",
                                      "total_energy_eh", "rel_energy_kcal"]
    assert lines[1].startswith("water,converged,GFN2-xTB,6.7.1,3,-5.0,0.0")
    assert len(ch.read_text().splitlines()) == 4


def test_cli_no_files_returns_2(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr("sys.argv", ["tool.py", str(tmp_path)])
    assert tool.main() == 2
    assert "no xtb output files" in capsys.readouterr().err
