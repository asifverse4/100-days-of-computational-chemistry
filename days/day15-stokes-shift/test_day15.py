import csv
from pathlib import Path

import numpy as np
import pytest
import tool

EXAMPLE = Path(__file__).parent / "example"


def gaussian(x, center, sigma, height=1.0):
    return height * np.exp(-0.5 * ((x - center) / sigma) ** 2)


def two_band_absorption(noise=0.0, seed=0):
    x = np.arange(220.0, 521.0, 1.0)
    y = gaussian(x, 240, 18, 1.0) + gaussian(x, 345, 26, 0.42)
    if noise:
        y = y + np.random.default_rng(seed).normal(0, noise, x.size)
    return tool.Spectrum("abs", x, y)


def run(monkeypatch, *argv):
    monkeypatch.setattr("sys.argv", ["tool.py", *map(str, argv)])
    return tool.main()


# ---------- the calculation ----------
def test_exact_values():
    r = tool.stokes(400, 500)
    assert (r["shift_nm"], r["shift_cm1"]) == (100, pytest.approx(5000.0))
    assert r["shift_ev"] == pytest.approx(1239.841984 * (1 / 400 - 1 / 500))
    r = tool.stokes(340, 450)
    assert r["shift_cm1"] == pytest.approx(1e7 / 340 - 1e7 / 450) and r["shift_cm1"] == pytest.approx(7189.5, abs=0.1)


def test_zero_and_negative_shift_and_validation():
    assert tool.stokes(450, 450)["shift_cm1"] == 0
    assert tool.stokes(500, 400)["shift_nm"] < 0
    with pytest.raises(ValueError, match="positive"):
        tool.stokes(0, 400)


def test_same_shift_in_all_three_units_is_consistent():
    r = tool.stokes(350, 470)
    assert r["shift_ev"] == pytest.approx(r["shift_cm1"] * 1239.841984 / 1e7)  # 1 eV = 8065.54 cm-1


def test_uncertainty_matches_a_simulation():
    rng = np.random.default_rng(1)
    a = 350 + rng.normal(0, 2.0, 20000)
    e = 470 + rng.normal(0, 2.0, 20000)
    simulated = np.std(1e7 / a - 1e7 / e)
    assert tool.stokes(350, 470, 2.0, 2.0)["shift_cm1_unc"] == pytest.approx(simulated, rel=0.03)
    assert tool.stokes(350, 470)["shift_cm1_unc"] == 0
    assert tool.stokes(350, 470, 3.0, 4.0)["shift_nm_unc"] == pytest.approx(5.0)


def test_unit_conversion_to_nm():
    assert tool.to_nm(400, "nm") == 400
    assert tool.to_nm(25000, "cm-1") == pytest.approx(400)
    assert tool.to_nm(3.0, "ev") == pytest.approx(1239.841984 / 3.0)
    with pytest.raises(ValueError, match="positive"):
        tool.to_nm(-1, "ev")


def test_warnings():
    assert tool.shift_warnings("s", tool.stokes(500, 400))[0].startswith("s: emission (400.0 nm) is at a shorter")
    assert "under 5 nm" in tool.shift_warnings("s", tool.stokes(450, 452))[0]
    assert tool.shift_warnings("s", tool.stokes(350, 450)) == []


# ---------- picking the bands ----------
def test_max_versus_longest_absorption_band():
    spec = two_band_absorption(noise=0.003, seed=2)
    assert tool.find_band(spec, "max") == pytest.approx(240, abs=1)
    assert tool.find_band(spec, "longest") == pytest.approx(345, abs=1)


def test_range_restricts_the_search():
    spec = two_band_absorption()
    assert tool.find_band(spec, "max", (300, 450)) == pytest.approx(345, abs=0.5)
    with pytest.raises(ValueError, match="fewer than 5 points"):
        tool.find_band(spec, "max", (300, 302))


def test_prominence_controls_which_bands_count():
    spec = two_band_absorption()
    assert tool.find_band(spec, "longest", prominence=0.1) == pytest.approx(345, abs=1)
    with pytest.raises(ValueError, match="no band found"):
        tool.find_band(tool.Spectrum("flat", spec.x, np.full_like(spec.x, 0.3)), "longest")


def test_smoothing_validation_and_noisy_band():
    spec = two_band_absorption(noise=0.02, seed=3)
    assert tool.find_band(spec, "longest", smooth=21, prominence=0.1) == pytest.approx(345, abs=2)
    with pytest.raises(ValueError, match="odd window"):
        tool.find_band(spec, "max", smooth=10)


def test_refined_position_beats_the_grid_on_a_coarse_scan():
    x = np.arange(400.0, 501.0, 5.0)
    spec = tool.Spectrum("coarse", x, gaussian(x, 452.3, 25))
    assert tool.find_band(spec, "max") == pytest.approx(452.3, abs=0.2)


# ---------- reading ----------
def write(tmp_path, text, name="s.csv"):
    p = tmp_path / name
    p.write_text(text)
    return p


def test_read_spectra_and_pick_column(tmp_path):
    body = "".join(f"{300 + 10 * i},{i / 10},{i / 5}\n" for i in range(12))
    spectra = tool.read_spectra(write(tmp_path, "nm,A,B\n" + body))
    assert [s.name for s in spectra] == ["A", "B"]
    assert tool.pick_column(spectra, None, "f").name == "A"
    assert tool.pick_column(spectra, "B", "f").name == "B"
    assert tool.pick_column(spectra, "2", "f").name == "B"
    with pytest.raises(ValueError, match="column 'C' not found"):
        tool.pick_column(spectra, "C", "f")
    with pytest.raises(ValueError, match="no numeric data"):
        tool.read_spectra(write(tmp_path, "text only\n"))
    with pytest.raises(ValueError, match="cannot read"):
        tool.read_spectra(tmp_path / "missing.csv")


def test_read_table_variants(tmp_path):
    rows = tool.read_table(EXAMPLE / "samples.csv", "nm")
    assert rows[0] == ("sample_A", 345.0, 450.0) and len(rows) == 4
    p = write(tmp_path, "345;450\n350;470\n")  # no header, no names: columns are positional
    with pytest.raises(ValueError):
        tool.read_table(p, "nm")
    p = write(tmp_path, "Sample,Absorption max (eV),Emission max (eV)\nQD1,3.6,2.8\n")
    (name, a, e), = tool.read_table(p, "ev")
    assert (name, a, e) == ("QD1", pytest.approx(1239.841984 / 3.6), pytest.approx(1239.841984 / 2.8))


def test_read_table_errors(tmp_path):
    with pytest.raises(ValueError, match="no rows found"):
        tool.read_table(write(tmp_path, "name,abs,em\n"), "nm")
    with pytest.raises(ValueError, match="cannot read"):
        tool.read_table(tmp_path / "missing.csv", "nm")
    with pytest.raises(ValueError, match="line 3"):
        tool.read_table(write(tmp_path, "name,abs,em\nA,300,400\nB,oops\n"), "nm")


# ---------- output ----------
@pytest.mark.parametrize("name, magic", [("p.png", b"\x89PNG"), ("p.pdf", b"%PDF"), ("p.svg", b"<?xm")])
def test_plot_formats(tmp_path, name, magic):
    abs_spec = two_band_absorption()
    x = np.arange(380.0, 701.0)
    em_spec = tool.Spectrum("em", x, gaussian(x, 450, 38))
    r = tool.stokes(345, 450)
    tool.plot_spectra(abs_spec, em_spec, r, tmp_path / name, 60)
    assert (tmp_path / name).read_bytes()[:4] == magic
    tool.plot_table([{"name": "A", **r}, {"name": "B", **tool.stokes(340, 430)}], tmp_path / ("t" + name), 60)
    assert (tmp_path / ("t" + name)).is_file()


def test_result_table_and_csv(tmp_path):
    rows = [{"name": "A", **tool.stokes(345, 450, 1.0, 1.0)}, {"name": "B", **tool.stokes(340, 430)}]
    assert "shift (cm-1)" in tool.result_table(rows) and "6763" in tool.result_table(rows)
    tool.write_csv(tmp_path / "r.csv", rows)
    parsed = list(csv.DictReader((tmp_path / "r.csv").open()))
    assert list(parsed[0]) == tool.RESULT_COLUMNS
    assert float(parsed[0]["shift_cm1"]) == pytest.approx(6763.3, abs=0.1) and parsed[0]["shift_cm1_unc"] != ""
    assert parsed[1]["shift_cm1_unc"] == ""


# ---------- CLI ----------
def test_cli_two_values(monkeypatch, capsys):
    assert run(monkeypatch, 345, 450) == 0
    out = capsys.readouterr().out
    assert "105.0" in out and "6763" in out and "0.839" in out
    assert run(monkeypatch, 400, 500, "--unc", 2) == 0
    assert "+/- 2 nm" in capsys.readouterr().out
    assert run(monkeypatch, 25000, 20000, "--input-unit", "cm-1") == 0
    assert "100.0" in capsys.readouterr().out


def test_cli_spectra_max_vs_longest(tmp_path, monkeypatch, capsys):
    base = ["--abs", EXAMPLE / "synthetic_absorption.csv", "--em", EXAMPLE / "synthetic_emission.csv"]
    assert run(monkeypatch, *base, "--csv", tmp_path / "r.csv") == 0
    assert "239.8" in capsys.readouterr().out  # the tallest band is the 240 nm one
    assert run(monkeypatch, *base, "--abs-peak", "longest", "-o", tmp_path / "s.png", "--dpi", 60) == 0
    assert "344.6" in capsys.readouterr().out and (tmp_path / "s.png").read_bytes()[:4] == b"\x89PNG"
    assert run(monkeypatch, *base, "--abs-range", 300, 450, "--em-range", 400, 600, "--smooth", 9) == 0


def test_cli_table_with_plot_and_csv(tmp_path, monkeypatch, capsys):
    assert run(monkeypatch, "--table", EXAMPLE / "samples.csv", "-o", tmp_path / "t.png", "--csv", tmp_path / "r.csv",
               "--dpi", 60) == 0
    out = capsys.readouterr().out
    assert "sample_D" in out and "7778" in out
    assert len(list(csv.DictReader((tmp_path / "r.csv").open()))) == 4


def test_cli_warns_about_anti_stokes(monkeypatch, capsys):
    assert run(monkeypatch, 500, 400) == 0
    assert "shorter wavelength" in capsys.readouterr().err


def test_cli_errors(tmp_path, monkeypatch, capsys):
    assert run(monkeypatch) == 2  # nothing given
    assert run(monkeypatch, 345) == 2  # one value
    assert run(monkeypatch, 345, 450, "-o", tmp_path / "x.png") == 2  # plot needs spectra or a table
    assert run(monkeypatch, 345, 450, "--table", EXAMPLE / "samples.csv") == 2  # two inputs at once
    assert run(monkeypatch, "--abs", EXAMPLE / "synthetic_absorption.csv") == 2  # --em missing
    assert run(monkeypatch, "--abs", EXAMPLE / "synthetic_absorption.csv", "--em", EXAMPLE / "synthetic_emission.csv",
               "--abs-range", 400, 300) == 2
    assert run(monkeypatch, "--table", EXAMPLE / "samples.csv", "-o", tmp_path / "x.jpg") == 2
    assert run(monkeypatch, 345, 450, "--unc", -1) == 2
    err = capsys.readouterr().err
    assert "choose one input" in err and "exactly two values" in err and "-o needs spectra" in err
    assert "--em must be given" in err or "must be given together" in err
    assert "LOW below HIGH" in err and "must end in" in err and "must not be negative" in err
