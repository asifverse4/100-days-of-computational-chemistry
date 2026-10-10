import csv
from pathlib import Path

import numpy as np
import pytest
import tool

EXAMPLE = Path(__file__).parent / "example"
K = 1.0e7
A_VALUES = np.array([0.012, 0.025, 0.041, 0.058, 0.077, 0.095])


def measurement(a, qy, n, noise=0.0, seed=0, offset=0.0):
    """Intensities that follow F = K * QY * f / n^2 (f = 1 - 10^-A), as the IUPAC equation implies."""
    a = np.asarray(a, dtype=float)
    f = 1 - 10 ** (-a)
    rng = np.random.default_rng(seed)
    return tool.Measurement(a, K * qy * f / n**2 * (1 + (rng.normal(0, noise, a.size) if noise else 0)) + offset)


def run(monkeypatch, *argv):
    monkeypatch.setattr("sys.argv", ["tool.py", *map(str, argv)])
    return tool.main()


# ---------- reference data ----------
@pytest.mark.parametrize("key, qy, unc, solvent", [
    ("quinine-sulfate", 0.546, None, "water"), ("quinine-hclo4", 0.60, 0.02, "water"),
    ("fluorescein", 0.925, 0.015, "water"), ("rhodamine-6g", 0.95, 0.015, "ethanol"),
    ("rhodamine-6g-water", 0.92, 0.02, "water"),
])
def test_standards_match_the_iupac_table(key, qy, unc, solvent):
    s = tool.STANDARDS[key]
    assert (s.qy, s.qy_unc, s.solvent) == (qy, unc, solvent)
    assert "Brouwer" in s.reference and s.excitation and s.usable_excitation and s.emission


def test_every_standard_solvent_has_a_refractive_index():
    assert all(s.solvent in tool.SOLVENTS for s in tool.STANDARDS.values())
    assert tool.SOLVENTS["water"] == 1.333 and tool.SOLVENTS["ethanol"] == 1.361


def test_resolve_n():
    assert tool.resolve_n("Ethanol", None, "sample") == 1.361
    assert tool.resolve_n("anything", 1.5, "sample") == 1.5
    with pytest.raises(ValueError, match="unknown solvent"):
        tool.resolve_n("hexanol", None, "sample")
    with pytest.raises(ValueError, match="--sample-solvent"):
        tool.resolve_n(None, None, "sample")
    with pytest.raises(ValueError, match="at least 1"):
        tool.resolve_n(None, 0.9, "ref")


def test_absorption_factor():
    a = np.array([0.1, 1.0])
    assert tool.absorption_factor(a, False) == pytest.approx([0.2056717653, 0.9])
    assert tool.absorption_factor(a, True).tolist() == [0.1, 1.0]


# ---------- the calculation ----------
def single(ref_qy=0.5, ref=(0.05, 1000.0), sample=(0.05, 1000.0), n_ref=1.333, n_sample=1.333, **kw):
    r = tool.Measurement(np.array([ref[0]]), np.array([ref[1]]))
    s = tool.Measurement(np.array([sample[0]]), np.array([sample[1]]))
    return tool.compute(r, s, ref_qy=ref_qy, ref_qy_unc=kw.pop("ref_qy_unc", None), n_ref=n_ref, n_sample=n_sample, **kw)


def test_single_point_basic_ratios():
    assert single()["qy"] == pytest.approx(0.5)
    assert single(sample=(0.05, 400.0))["qy"] == pytest.approx(0.2)  # fewer photons -> lower yield
    assert single(n_ref=1.0, n_sample=1.5)["qy"] == pytest.approx(0.5 * 2.25)  # refractive index squared
    # half the light absorbed by the sample with the same emission means twice the yield
    f_ref = 1 - 10 ** -0.1
    a_half = -np.log10(1 - f_ref / 2)
    assert single(ref=(0.1, 1000.0), sample=(a_half, 1000.0))["qy"] == pytest.approx(1.0, rel=1e-9)


def test_single_point_absorbance_flag_uses_a_directly():
    assert single(ref=(0.1, 1000.0), sample=(0.05, 1000.0), use_absorbance=True)["qy"] == pytest.approx(1.0)
    assert single(ref=(0.1, 1000.0), sample=(0.05, 1000.0))["qy"] != pytest.approx(1.0, abs=0.01)


def test_single_point_has_no_statistical_uncertainty_but_carries_the_reference_one():
    assert single()["uncertainty"] is None
    assert single(ref_qy=0.6, ref_qy_unc=0.03)["uncertainty"] == pytest.approx(0.03)


def test_gradient_recovers_the_true_yield():
    ref = measurement(A_VALUES, 0.925, 1.333)
    sample = measurement(A_VALUES * 0.9, 0.35, 1.361)
    result = tool.compute(ref, sample, ref_qy=0.925, ref_qy_unc=None, n_ref=1.333, n_sample=1.361)
    assert result["method"] == "gradient" and result["qy"] == pytest.approx(0.35, rel=1e-6)
    assert result["uncertainty"] == pytest.approx(0, abs=1e-6)


def test_gradient_ignores_a_constant_background_with_free_intercept():
    ref = measurement(A_VALUES, 0.5, 1.333, offset=2000)
    sample = measurement(A_VALUES, 0.25, 1.333, offset=2000)
    result = tool.compute(ref, sample, ref_qy=0.5, ref_qy_unc=None, n_ref=1.333, n_sample=1.333)
    assert result["qy"] == pytest.approx(0.25, rel=1e-6)


def test_through_origin_fit():
    f = np.array([0.02, 0.05, 0.1])
    fit = tool.fit_gradient(f, 100 * f + np.array([1.0, -1.0, 0.5]), through_origin=True)
    manual = np.sum(f * (100 * f + np.array([1.0, -1.0, 0.5]))) / np.sum(f * f)
    assert fit.slope == pytest.approx(manual) and fit.intercept == 0.0
    with pytest.raises(ValueError, match="not positive"):
        tool.fit_gradient(f, -f)


def test_uncertainty_matches_the_scatter_of_repeated_noisy_experiments():
    qys, sigmas = [], []
    for seed in range(300):
        ref = measurement(A_VALUES, 0.9, 1.333, noise=0.02, seed=2 * seed)
        sample = measurement(A_VALUES, 0.3, 1.333, noise=0.02, seed=2 * seed + 1)
        result = tool.compute(ref, sample, ref_qy=0.9, ref_qy_unc=None, n_ref=1.333, n_sample=1.333)
        qys.append(result["qy"])
        sigmas.append(result["uncertainty"])
    assert np.mean(qys) == pytest.approx(0.3, abs=0.005)
    assert np.mean(sigmas) == pytest.approx(np.std(qys), rel=0.25)


def test_reference_uncertainty_adds_to_the_result():
    ref, sample = measurement(A_VALUES, 0.9, 1.333, 0.01, 1), measurement(A_VALUES, 0.3, 1.333, 0.01, 2)
    kw = {"ref_qy": 0.9, "n_ref": 1.333, "n_sample": 1.333}
    plain = tool.compute(ref, sample, ref_qy_unc=None, **{k: v for k, v in kw.items()})
    with_unc = tool.compute(ref, sample, ref_qy_unc=0.05, **kw)
    assert with_unc["uncertainty"] > plain["uncertainty"]


@pytest.mark.parametrize("nr, ns", [(1, 4), (4, 1), (2, 2), (2, 5)])
def test_wrong_number_of_points(nr, ns):
    ref, sample = measurement(A_VALUES[:nr], 0.9, 1.333), measurement(A_VALUES[:ns], 0.3, 1.333)
    with pytest.raises(ValueError, match="either one point each"):
        tool.compute(ref, sample, ref_qy=0.9, ref_qy_unc=None, n_ref=1.333, n_sample=1.333)


def test_reference_qy_must_be_a_fraction():
    for bad in (0, -0.1, 1.5):
        with pytest.raises(ValueError, match="between 0 and 1"):
            single(ref_qy=bad)


# ---------- warnings ----------
def notes(ref, sample, **kw):
    result = tool.compute(ref, sample, ref_qy=0.5, ref_qy_unc=None, n_ref=1.333, n_sample=1.333, **kw)
    return tool.warnings_for(ref, sample, result)


def test_no_warnings_for_a_clean_measurement():
    assert notes(measurement(A_VALUES, 0.5, 1.333), measurement(A_VALUES, 0.4, 1.333)) == []


def test_high_absorbance_warning():
    high = np.array([0.02, 0.06, 0.12, 0.2])
    assert any("above 0.1" in n for n in notes(measurement(high, 0.5, 1.333), measurement(high, 0.4, 1.333)))


def test_nonlinear_and_intercept_and_range_warnings():
    ref = measurement(A_VALUES, 0.5, 1.333)
    curved = tool.Measurement(A_VALUES, ref.intensity * np.array([1, 1, 1, 0.9, 0.8, 0.6]))
    assert any("not very linear" in n for n in notes(ref, curved))
    offset = measurement(A_VALUES, 0.4, 1.333, offset=5e5)
    assert any("large intercept" in n for n in notes(ref, offset))
    narrow = measurement(np.array([0.05, 0.06, 0.07]), 0.4, 1.333)
    assert any("narrow range" in n for n in notes(ref, narrow))


def test_quantum_yield_above_one_warning():
    ref = measurement(A_VALUES, 0.5, 1.333)
    sample = measurement(A_VALUES, 1.5, 1.333)  # computed yield 1.5 with a 0.5 reference
    assert any("above 1" in n for n in notes(ref, sample))


# ---------- reading ----------
def test_read_by_header_names_and_position(tmp_path):
    p = tmp_path / "m.csv"
    p.write_text("integrated area,absorbance (A)\n1000,0.02\n2000,0.04\n3000,0.06\n")
    m = tool.read_measurement(p)
    assert m.a.tolist() == [0.02, 0.04, 0.06] and m.intensity.tolist() == [1000, 2000, 3000]
    p.write_text("0.02;1000,5\n0.04;2000,5\n")
    m = tool.read_measurement(p)
    assert m.intensity.tolist() == [1000.5, 2000.5]


def test_read_errors(tmp_path):
    p = tmp_path / "m.csv"
    p.write_text("a,b\nx,y\n")
    with pytest.raises(ValueError, match="no numeric data"):
        tool.read_measurement(p)
    p.write_text("0.0,100\n0.1,200\n")
    with pytest.raises(ValueError, match="above 0"):
        tool.read_measurement(p)
    p.write_text("0.1,-5\n0.2,10\n")
    with pytest.raises(ValueError, match="not be negative"):
        tool.read_measurement(p)
    with pytest.raises(ValueError, match="cannot read"):
        tool.read_measurement(tmp_path / "missing.csv")


# ---------- output ----------
@pytest.mark.parametrize("name, magic", [("p.png", b"\x89PNG"), ("p.pdf", b"%PDF"), ("p.svg", b"<?xm")])
def test_plot_formats_gradient_and_single_point(tmp_path, name, magic):
    ref, sample = measurement(A_VALUES, 0.9, 1.333), measurement(A_VALUES, 0.3, 1.333)
    result = tool.compute(ref, sample, ref_qy=0.9, ref_qy_unc=None, n_ref=1.333, n_sample=1.333)
    tool.plot_fit(ref, sample, result, tmp_path / name, "x", 60)
    assert (tmp_path / name).read_bytes()[:4] == magic
    one = single()
    tool.plot_fit(tool.Measurement(one["fr"], np.array([1000.0])), tool.Measurement(one["fs"], np.array([1000.0])),
                  one, tmp_path / "one.png", "x", 60)
    assert (tmp_path / "one.png").is_file()


def test_standards_table_lists_every_standard():
    text = tool.standards_table()
    assert all(key in text for key in tool.STANDARDS) and "Brouwer" in text and "0.925" in text


# ---------- CLI ----------
def test_cli_example_recovers_the_demo_yield(tmp_path, monkeypatch, capsys):
    out_csv = tmp_path / "qy.csv"
    assert run(monkeypatch, "--ref", EXAMPLE / "synthetic_reference.csv", "--sample", EXAMPLE / "synthetic_sample.csv",
               "--standard", "fluorescein", "--sample-solvent", "ethanol", "--plot", tmp_path / "p.png",
               "--csv", out_csv, "--dpi", 60) == 0
    out = capsys.readouterr().out
    assert "gradient method" in out and "fluorescein in 0.1 M NaOH" in out
    row = next(csv.DictReader(out_csv.open()))
    assert float(row["quantum_yield"]) == pytest.approx(0.35, abs=0.02)
    assert float(row["uncertainty"]) > 0 and row["method"] == "gradient"
    assert (tmp_path / "p.png").read_bytes()[:4] == b"\x89PNG"


def test_cli_single_point_and_custom_reference(monkeypatch, capsys):
    assert run(monkeypatch, "--ref-point", 0.05, 1.2e6, "--sample-point", 0.05, 4.0e5, "--ref-qy", 0.6,
               "--ref-solvent", "water", "--sample-solvent", "water") == 0
    out = capsys.readouterr().out
    assert "single-point" in out and "QY sample = 0.200" in out and "custom reference" in out


def test_cli_lists(monkeypatch, capsys):
    assert run(monkeypatch, "--list-standards") == 0
    assert "rhodamine-6g" in capsys.readouterr().out
    assert run(monkeypatch, "--list-solvents") == 0
    assert "ethanol" in capsys.readouterr().out


def test_cli_warnings_go_to_stderr(monkeypatch, capsys):
    assert run(monkeypatch, "--ref-point", 0.3, 1.0e6, "--sample-point", 0.3, 9.0e6, "--ref-qy", 0.9,
               "--ref-n", 1.333, "--sample-n", 1.333) == 0
    err = capsys.readouterr().err
    assert "above 0.1" in err and "above 1" in err


def test_cli_errors(tmp_path, monkeypatch, capsys):
    ex = ["--ref", EXAMPLE / "synthetic_reference.csv", "--sample", EXAMPLE / "synthetic_sample.csv"]
    assert run(monkeypatch, *ex, "--sample-solvent", "ethanol") == 2  # no reference value
    assert run(monkeypatch, *ex, "--standard", "fluorescein") == 2  # no sample solvent
    assert run(monkeypatch, *ex, "--standard", "fluorescein", "--sample-solvent", "ethanol", "--plot", tmp_path / "p.jpg") == 2
    assert run(monkeypatch, "--sample", EXAMPLE / "synthetic_sample.csv", "--standard", "fluorescein",
               "--sample-solvent", "ethanol") == 2  # no reference data
    assert run(monkeypatch, "--ref", tmp_path / "missing.csv", "--ref-point", 0.1, 5.0, "--sample-point", 0.1, 5.0,
               "--ref-qy", 0.5, "--ref-n", 1.3, "--sample-n", 1.3) == 2  # both file and point
    err = capsys.readouterr().err
    assert "reference quantum yield" in err and "--sample-solvent" in err and "--plot must end in" in err
    assert "not both" in err
    with pytest.raises(SystemExit):
        run(monkeypatch, *ex, "--standard", "unknown-dye")
