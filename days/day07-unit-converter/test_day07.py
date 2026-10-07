from itertools import permutations

import pytest
import tool


@pytest.mark.parametrize("value, src, dst, expected", [
    (1, "eh", "ev", 27.211386245988),
    (1, "eh", "kcal/mol", 627.5094740631),
    (1, "eh", "kj/mol", 2625.4996394799),
    (1, "eh", "cm-1", 219474.6313632),
    (1, "ev", "cm-1", 8065.543937349),
    (1, "ev", "nm", 1239.841984332),
    (1, "kcal/mol", "kj/mol", 4.184),
    (1, "kcal/mol", "cm-1", 349.7550882),
])
def test_known_values(value, src, dst, expected):
    assert tool.convert(value, src, dst) == pytest.approx(expected, rel=1e-8)


def test_thermal_energy_at_room_temperature():
    assert tool.convert(298.15, "k", "kcal/mol") == pytest.approx(0.5925, abs=1e-4)


def test_optical_gap_from_absorption_edge():
    assert tool.convert(450, "nm", "ev") == pytest.approx(2.7552, abs=1e-4)
    assert tool.convert(2.7552, "ev", "nm") == pytest.approx(450, abs=0.05)


def test_wavenumber_to_frequency():
    assert tool.convert(1000, "cm-1", "thz") == pytest.approx(29.9792458, rel=1e-9)


@pytest.mark.parametrize("src, dst", list(permutations(tool.ORDER, 2)))
def test_round_trip_all_pairs(src, dst):
    back = tool.convert(tool.convert(7.5, src, dst), dst, src)
    assert back == pytest.approx(7.5, rel=1e-12)


def test_negative_values_allowed_for_energies_and_imaginary_frequencies():
    assert tool.convert(-0.5, "eh", "ev") == pytest.approx(-13.6056931, rel=1e-6)
    assert tool.convert(-100, "cm-1", "kj/mol") < 0


@pytest.mark.parametrize("value, src, dst", [(0, "nm", "ev"), (-5, "nm", "ev"), (0, "ev", "nm"), (-1, "ev", "nm")])
def test_inverse_unit_needs_positive_values(value, src, dst):
    with pytest.raises(ValueError):
        tool.convert(value, src, dst)


@pytest.mark.parametrize("name, key", [
    ("Hartree", "eh"), ("Eh", "eh"), ("eV", "ev"), ("kcal/mol", "kcal/mol"), ("kcal mol-1", "kcal/mol"),
    ("kJ/mol", "kj/mol"), ("cm^-1", "cm-1"), ("cm\u207b\u00b9", "cm-1"), ("1/cm", "cm-1"),
    ("nm", "nm"), ("THz", "thz"), ("K", "k"),
])
def test_unit_aliases(name, key):
    assert tool.canonical(name) == key


def test_unknown_unit_lists_supported_units():
    with pytest.raises(ValueError, match="Supported: Eh, eV"):
        tool.canonical("furlong")


def run(monkeypatch, *argv):
    monkeypatch.setattr("sys.argv", ["tool.py", *map(str, argv)])
    return tool.main()


def test_cli_single_and_all_units(monkeypatch, capsys):
    assert run(monkeypatch, 1, "eh", "ev") == 0
    assert capsys.readouterr().out.strip() == "1 Eh = 27.2114 eV"
    assert run(monkeypatch, 2.5, "ev") == 0
    out = capsys.readouterr().out
    assert "2.5 eV = 495.937 nm" in out
    assert "2.5 eV = 0.0918733 Eh" in out
    assert out.count("\n") == 7  # every unit except the source


def test_cli_negative_energy_to_nm_is_marked_not_available(monkeypatch, capsys):
    assert run(monkeypatch, -1, "eh", "nm") == 0
    assert "n/a" in capsys.readouterr().out


def test_cli_errors_return_2(monkeypatch, capsys):
    assert run(monkeypatch, "abc", "ev", "eh") == 2
    assert "not a number" in capsys.readouterr().err
    assert run(monkeypatch, 1, "foo", "ev") == 2
    assert run(monkeypatch, "--csv", "x.csv") == 2


def test_csv_mode(tmp_path, monkeypatch):
    src = tmp_path / "e.csv"
    src.write_text("name,energy_eh\na,-0.5\nb,\nc,-1.0\n")
    out = tmp_path / "out.csv"
    assert run(monkeypatch, "--csv", src, "--column", "energy_eh", "--from", "eh",
               "--to", "ev", "kcal/mol", "-o", out) == 0
    lines = out.read_text().splitlines()
    assert lines[0] == "name,energy_eh,energy_eh_ev,energy_eh_kcal_mol"
    assert float(lines[1].split(",")[2]) == pytest.approx(-13.6056931, rel=1e-6)
    assert lines[2] == "b,,,"
    assert float(lines[3].split(",")[3]) == pytest.approx(-627.5094740631, rel=1e-9)


def test_csv_errors(tmp_path, monkeypatch, capsys):
    src = tmp_path / "e.csv"
    src.write_text("name,energy_eh\na,oops\n")
    assert run(monkeypatch, "--csv", src, "--column", "energy_eh", "--from", "eh", "--to", "ev") == 2
    assert "line 2" in capsys.readouterr().err
    assert run(monkeypatch, "--csv", src, "--column", "missing", "--from", "eh", "--to", "ev") == 2
    assert "not found" in capsys.readouterr().err
    assert run(monkeypatch, "--csv", tmp_path / "nope.csv", "--column", "a", "--from", "eh", "--to", "ev") == 2
