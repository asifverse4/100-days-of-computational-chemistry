import csv
from pathlib import Path

import numpy as np
import pytest
import tool

EXAMPLE = Path(__file__).parent / "example" / "synthetic_tauc.csv"
WL = np.arange(250.0, 701.0, 1.0)
E = tool.HC_EV_NM / WL


def absorbance_for(transition, gap, b=170.0, offset=0.0, noise=0.0, seed=0):
    """Absorbance (1 cm path) that follows the Tauc law of the given transition above `gap`."""
    above = np.clip(E - gap, 0, None)
    alpha_e = {"direct": np.sqrt(b * above), "indirect": b * above**2,
               "direct-forbidden": b * above**1.5, "indirect-forbidden": b * above**3}[transition]
    a = alpha_e / E / tool.LN10 + offset
    return a + np.random.default_rng(seed).normal(0, noise, a.size) if noise else a


def spectrum(transition="direct", gap=3.0, **kw):
    return tool.Spectrum("s", WL.copy(), absorbance_for(transition, gap, **kw))


def run_analyze(spec, **overrides):
    opts = {"transition": "direct", "signal": "absorbance", "path_length_cm": 1.0, "x_unit": "nm",
            "zero_at_nm": None, "fit_range": None, "emin": None, "emax": None, "width_ev": 0.25, "r2_min": 0.99}
    opts.update(overrides)
    return tool.analyze(spec, **opts)


def run(monkeypatch, *argv):
    monkeypatch.setattr("sys.argv", ["tool.py", *map(str, argv)])
    return tool.main()


# ---------- conversions ----------
def test_to_energy_nm_and_ev():
    out = tool.to_energy(tool.Spectrum("s", np.array([400.0, 800.0, 200.0]), np.array([1.0, 2.0, 3.0])), "nm")
    assert out.x.tolist() == pytest.approx([1239.841984 / 800, 1239.841984 / 400, 1239.841984 / 200])
    assert out.y.tolist() == [2.0, 1.0, 3.0]  # y follows its own x
    ev = tool.to_energy(tool.Spectrum("s", np.array([3.0, 2.0]), np.array([1.0, 2.0])), "ev")
    assert ev.x.tolist() == [2.0, 3.0]
    with pytest.raises(ValueError, match="positive"):
        tool.to_energy(tool.Spectrum("s", np.array([0.0, 1.0]), np.array([1.0, 1.0])), "nm")


def test_absorption_proxy():
    assert tool.absorption_proxy(np.array([1.0]), "absorbance", 1.0)[0] == pytest.approx(2.302585)
    assert tool.absorption_proxy(np.array([1.0]), "absorbance", 0.5)[0] == pytest.approx(4.60517)
    with pytest.raises(ValueError, match="path-length"):
        tool.absorption_proxy(np.array([1.0]), "absorbance", 0.0)


def test_kubelka_munk():
    assert tool.absorption_proxy(np.array([0.5]), "reflectance", 1.0)[0] == pytest.approx(0.25)
    assert tool.absorption_proxy(np.array([50.0]), "reflectance-percent", 1.0)[0] == pytest.approx(0.25)
    for bad in (0.0, 1.0, 1.5):
        with pytest.raises(ValueError, match="between 0 and 1"):
            tool.absorption_proxy(np.array([bad]), "reflectance", 1.0)


def test_tauc_exponents_and_clipping():
    e, a = np.array([2.0, 2.0, 2.0]), np.array([2.0, -1.0, 2.0])
    assert tool.tauc_y(e, a, "direct").tolist() == [16.0, 0.0, 16.0]
    assert tool.tauc_y(e[:1], a[:1], "indirect")[0] == pytest.approx(2.0)
    assert tool.tauc_y(e[:1], a[:1], "direct-forbidden")[0] == pytest.approx(4 ** (2 / 3))
    assert tool.tauc_y(e[:1], a[:1], "indirect-forbidden")[0] == pytest.approx(4 ** (1 / 3))


# ---------- the linear fit ----------
def test_linear_fit_exact_line():
    e = np.linspace(3.2, 3.6, 20)
    fit = tool.linear_fit(e, 2 * e - 6)
    assert fit.gap == pytest.approx(3.0) and fit.sigma == pytest.approx(0, abs=1e-9)
    assert fit.r2 == pytest.approx(1.0) and fit.n == 20


def test_linear_fit_rejects_non_positive_slope():
    e = np.linspace(3, 4, 10)
    with pytest.raises(ValueError, match="not positive"):
        tool.linear_fit(e, -e)


def test_uncertainty_matches_repeated_noisy_fits():
    """The propagated sigma should match the scatter of the intercept over many noisy repeats."""
    e = np.linspace(3.2, 3.6, 40)
    rng = np.random.default_rng(5)
    gaps, sigmas = [], []
    for _ in range(400):
        fit = tool.linear_fit(e, 100 * (e - 3.0) + rng.normal(0, 2.0, e.size))
        gaps.append(fit.gap)
        sigmas.append(fit.sigma)
    assert np.mean(gaps) == pytest.approx(3.0, abs=0.01)
    assert np.mean(sigmas) == pytest.approx(np.std(gaps), rel=0.2)


# ---------- finding the band gap ----------
@pytest.mark.parametrize("transition, gap", [("direct", 3.0), ("indirect", 2.4), ("direct-forbidden", 3.3),
                                             ("indirect-forbidden", 2.9)])
def test_recovers_known_gap_for_each_transition(transition, gap):
    row, *_ = run_analyze(spectrum(transition, gap, noise=0.001, seed=4), transition=transition)
    assert row["band_gap_ev"] == pytest.approx(gap, abs=0.04)
    assert row["method"] == "auto" and row["r_squared"] >= 0.99 and row["n_points"] >= 5
    assert row["warning"] == ""


def test_zero_at_removes_baseline_offset():
    spec = spectrum("direct", 3.0, offset=0.05, noise=0.001, seed=2)
    without, *_ = run_analyze(spec)
    with_zero, *_ = run_analyze(spec, zero_at_nm=700)
    assert abs(with_zero["band_gap_ev"] - 3.0) < 0.04
    assert abs(with_zero["band_gap_ev"] - 3.0) < abs(without["band_gap_ev"] - 3.0)


def test_diffuse_reflectance_with_kubelka_munk():
    f = np.maximum(np.sqrt(170 * np.clip(E - 3.2, 0, None)) / E, 1e-4)
    reflectance = 1 + f - np.sqrt(f**2 + 2 * f)  # inverse of F(R) = (1-R)^2 / (2R)
    spec = tool.Spectrum("powder", WL.copy(), reflectance * 100)
    row, *_ = run_analyze(spec, signal="reflectance-percent")
    assert row["band_gap_ev"] == pytest.approx(3.2, abs=0.04)


def test_manual_range_is_used_and_validated():
    spec = spectrum("direct", 3.0)
    row, *_ = run_analyze(spec, fit_range=(3.3, 3.7))
    assert row["method"] == "manual" and row["band_gap_ev"] == pytest.approx(3.0, abs=0.01)
    assert 3.29 <= row["fit_emin_ev"] and row["fit_emax_ev"] <= 3.71
    with pytest.raises(ValueError, match="fewer than 5 points"):
        run_analyze(spec, fit_range=(3.3, 3.31))


def test_warning_when_intercept_is_far_from_fit_range():
    row, *_ = run_analyze(spectrum("direct", 3.0), fit_range=(4.4, 4.9))
    assert "more than 1 eV below" in row["warning"]


def test_search_range_limits_are_respected():
    row, *_ = run_analyze(spectrum("direct", 3.0), emin=3.4, emax=4.0)
    assert row["fit_emin_ev"] >= 3.4 and row["fit_emax_ev"] <= 4.0


def test_auto_fit_errors_and_threshold_flag():
    e = np.linspace(2.0, 4.0, 200)
    with pytest.raises(ValueError, match="narrower than --window"):
        tool.auto_fit(e[:5], e[:5], width_ev=0.25)
    with pytest.raises(ValueError, match="no linear region"):
        tool.auto_fit(e, -e)
    fit, met = tool.auto_fit(e, e**4, r2_min=0.99999)  # strongly curved, cannot be linear
    assert met is False and fit.r2 < 0.99999


# ---------- reading ----------
def test_read_spectra_formats(tmp_path):
    rows = "".join(f"{300 + 10 * i},{i / 10}\n" for i in range(12))
    p = tmp_path / "a.csv"
    p.write_text("Wavelength (nm),Abs\n" + rows)
    (spec,) = tool.read_spectra(p)
    assert spec.name == "Abs" and spec.x.size == 12
    p.write_text("1,1\n2,2\n3,3\n")
    with pytest.raises(ValueError, match="only 3 data points"):
        tool.read_spectra(p)
    with pytest.raises(ValueError, match="cannot read"):
        tool.read_spectra(tmp_path / "missing.csv")


# ---------- output ----------
@pytest.mark.parametrize("name, magic", [("t.png", b"\x89PNG"), ("t.pdf", b"%PDF"), ("t.svg", b"<?xm")])
def test_plot_formats_and_three_panels(tmp_path, name, magic):
    items = [run_analyze(spectrum("direct", g)) for g in (2.8, 3.0, 3.2)]
    tool.plot_tauc(items, tmp_path / name, title="t", dpi=60, figsize=None)
    assert (tmp_path / name).read_bytes()[:4] == magic


def test_result_table_and_csv(tmp_path):
    row, *_ = run_analyze(spectrum("direct", 3.0))
    table = tool.result_table([row])
    assert "Eg (eV)" in table and "direct" in table and f"{row['band_gap_ev']:.3f}" in table
    tool.write_csv(tmp_path / "r.csv", [row])
    parsed = list(csv.DictReader((tmp_path / "r.csv").open()))
    assert list(parsed[0]) == tool.RESULT_COLUMNS and float(parsed[0]["band_gap_ev"]) == pytest.approx(3.0, abs=0.04)


# ---------- CLI ----------
def test_cli_example_recovers_the_demo_gaps(tmp_path, monkeypatch, capsys):
    out_csv = tmp_path / "gaps.csv"
    assert run(monkeypatch, EXAMPLE, "-o", tmp_path / "t.png", "--csv", out_csv, "--dpi", 60) == 0
    assert (tmp_path / "t.png").read_bytes()[:4] == b"\x89PNG"
    rows = {r["sample"]: r for r in csv.DictReader(out_csv.open())}
    assert float(rows["sample_A"]["band_gap_ev"]) == pytest.approx(3.10, abs=0.02)
    assert float(rows["sample_B"]["band_gap_ev"]) == pytest.approx(2.85, abs=0.02)
    assert "[ok] 2 spectra" in capsys.readouterr().out


def test_cli_options(tmp_path, monkeypatch, capsys):
    assert run(monkeypatch, EXAMPLE, "-o", tmp_path / "t.png", "--columns", "sample_A", "--fit-range", 3.4, 3.8,
               "--zero-at", 700, "--dpi", 60) == 0
    out = capsys.readouterr().out
    assert "sample_A" in out and "sample_B" not in out
    assert run(monkeypatch, EXAMPLE, "-o", tmp_path / "t.svg", "--transition", "indirect", "--window", 0.4) in (0, 2)


def test_cli_errors(tmp_path, monkeypatch, capsys):
    assert run(monkeypatch, EXAMPLE, "-o", tmp_path / "t.jpg") == 2
    assert run(monkeypatch, tmp_path / "missing.csv", "-o", tmp_path / "t.png") == 2
    assert run(monkeypatch, EXAMPLE, "-o", tmp_path / "t.png", "--columns", "nope") == 2
    assert run(monkeypatch, EXAMPLE, "-o", tmp_path / "t.png", "--signal", "reflectance") == 2
    assert run(monkeypatch, EXAMPLE, "-o", tmp_path / "t.png", "--x-unit", "ev", "--zero-at", 700) == 2
    err = capsys.readouterr().err
    assert "must end in" in err and "column(s) not found: nope" in err and "between 0 and 1" in err
    with pytest.raises(SystemExit):
        run(monkeypatch, EXAMPLE, "--transition", "sideways")
