import csv
from pathlib import Path

import numpy as np
import pytest
import tool

EXAMPLE = Path(__file__).parent / "example" / "synthetic_xrd.csv"
LAMBDA, K = 1.54056, 0.9
FCC_A = 4.05


def reflections(size_nm=20.0, strain=0.002, inst=0.08, heights=(1000, 450, 250, 280, 110, 90)):
    out = []
    for (hkl, h) in zip([(1, 1, 1), (2, 0, 0), (2, 2, 0), (3, 1, 1), (2, 2, 2), (4, 0, 0)], heights):
        d = FCC_A / np.sqrt(sum(i * i for i in hkl))
        theta = np.arcsin(LAMBDA / (2 * d))
        beta = (K * LAMBDA / (size_nm * 10) + 4 * strain * np.sin(theta)) / np.cos(theta)
        out.append((float(np.degrees(2 * theta)), float(np.degrees(beta)), float(np.hypot(np.degrees(beta), inst)), h))
    return out  # (2theta, true corrected FWHM, observed FWHM, height)


def pattern(noise=True, seed=3, **kw):
    x = np.arange(20.0, 100.0, 0.02)
    y = 60 + 0.4 * (100 - x)
    for tt, _, fwhm, h in reflections(**kw):
        y = y + h * tool.pseudo_voigt(x, 1.0, tt, fwhm, 0.5, 0.0, 0.0)
    return tool.Pattern(x, np.random.default_rng(seed).poisson(y).astype(float) if noise else y)


def run(monkeypatch, *argv):
    monkeypatch.setattr("sys.argv", ["tool.py", *map(str, argv)])
    return tool.main()


def write(tmp_path, text, name="s.csv"):
    p = tmp_path / name
    p.write_text(text)
    return p


# ---------- constants and the physics ----------
def test_wavelengths_match_the_published_tables():
    assert tool.WAVELENGTHS == {"cu-ka1": 1.54056, "cu-ka": 1.54184, "co-ka1": 1.78897, "co-ka": 1.79026,
                                "cr-ka1": 2.28970, "cr-ka": 2.29100, "mo-ka1": 0.70930, "mo-ka": 0.71073}


def test_bragg_d():
    assert tool.bragg_d(60.0, 1.54056) == pytest.approx(1.54056)  # sin(30 degrees) = 0.5
    assert tool.bragg_d(38.47, 1.54056) == pytest.approx(1.54056 / (2 * np.sin(np.radians(19.235))))
    assert tool.bragg_d(38.47, LAMBDA) == pytest.approx(2.3382, abs=1e-3)


def test_correct_fwhm_methods():
    assert tool.correct_fwhm(0.45, 0.10) == pytest.approx(np.sqrt(0.45**2 - 0.10**2))
    assert tool.correct_fwhm(0.45, 0.10, "lorentzian") == pytest.approx(0.35)
    assert tool.correct_fwhm(0.45, 0.10, "none") == 0.45
    assert tool.correct_fwhm(0.45, 0.0) == 0.45
    assert tool.correct_fwhm(0.10, 0.10) is None and tool.correct_fwhm(0.05, 0.10, "lorentzian") is None


def test_scherrer_size_hand_calculation():
    beta, theta = np.radians(0.45), np.radians(38.5) / 2
    expected = K * LAMBDA / (beta * np.cos(theta)) / 10
    assert tool.scherrer_size_nm(38.5, 0.45, LAMBDA, K) == pytest.approx(expected)
    assert tool.scherrer_size_nm(38.5, 0.45, LAMBDA, K) == pytest.approx(18.69, abs=0.02)


def test_scherrer_scaling():
    base = tool.scherrer_size_nm(40.0, 0.5, LAMBDA)
    assert tool.scherrer_size_nm(40.0, 1.0, LAMBDA) == pytest.approx(base / 2)
    assert tool.scherrer_size_nm(40.0, 0.5, LAMBDA, 0.94) == pytest.approx(base * 0.94 / 0.9)
    assert tool.scherrer_size_nm(40.0, 0.5, 2 * LAMBDA) == pytest.approx(2 * base)


def test_analyze_peak_with_and_without_instrument():
    plain = tool.analyze_peak(38.5, 0.45, LAMBDA, k=K, inst=0.0, method="gaussian")
    corrected = tool.analyze_peak(38.5, 0.45, LAMBDA, k=K, inst=0.2, method="gaussian")
    assert corrected["size_nm"] > plain["size_nm"]  # removing instrument width makes the peak narrower, the size larger
    assert tool.analyze_peak(38.5, 0.15, LAMBDA, k=K, inst=0.2, method="gaussian")["size_nm"] is None


def test_size_uncertainty_matches_a_simulation():
    rng = np.random.default_rng(1)
    sizes = [tool.scherrer_size_nm(40.0, tool.correct_fwhm(w, 0.2), LAMBDA) for w in rng.normal(0.6, 0.01, 20000)]
    row = tool.analyze_peak(40.0, 0.6, LAMBDA, k=K, inst=0.2, method="gaussian", fwhm_unc=0.01)
    assert row["size_unc_nm"] == pytest.approx(np.std(sizes), rel=0.05)


# ---------- Williamson-Hall ----------
def exact_rows(size_nm=20.0, strain=0.002):
    return [{"two_theta_deg": tt, "fwhm_corrected_deg": beta} for tt, beta, _, _ in reflections(size_nm, strain)]


def test_williamson_hall_recovers_exact_size_and_strain():
    wh = tool.williamson_hall(exact_rows(), LAMBDA, K)
    assert wh["size_nm"] == pytest.approx(20.0, rel=1e-6) and wh["strain"] == pytest.approx(0.002, rel=1e-6)
    assert wh["r2"] == pytest.approx(1.0) and wh["n"] == 6


def test_williamson_hall_edge_cases():
    with pytest.raises(ValueError, match="at least 3 peaks"):
        tool.williamson_hall(exact_rows()[:2], LAMBDA)
    rows = exact_rows()
    rows[0]["fwhm_corrected_deg"] = None  # not wider than the instrument: left out
    assert tool.williamson_hall(rows, LAMBDA)["n"] == 5
    flat = [{"two_theta_deg": tt, "fwhm_corrected_deg": 0.3 * (1 - 0.0004 * tt)} for tt in (30, 50, 70, 90)]
    assert tool.williamson_hall(flat, LAMBDA)["size_nm"] is None or tool.williamson_hall(flat, LAMBDA)["strain"] < 0


# ---------- peak shape and fitting ----------
def test_pseudo_voigt_height_and_fwhm():
    assert tool.pseudo_voigt(np.array([10.0]), 5.0, 10.0, 0.4, 0.3, 2.0, 0.0)[0] == pytest.approx(7.0)
    half = tool.pseudo_voigt(np.array([10.2]), 5.0, 10.0, 0.4, 0.3, 0.0, 0.0)[0]
    assert half == pytest.approx(2.5)  # half the height at x0 + FWHM/2 for any mixing


def test_fit_peak_recovers_position_and_width():
    x = np.arange(35.0, 42.0, 0.02)
    y = 100 + 800 * tool.pseudo_voigt(x, 1.0, 38.47, 0.6, 0.5, 0.0, 0.0)
    y = np.random.default_rng(2).poisson(y).astype(float)
    fit = tool.fit_peak(tool.Pattern(x, y), 38.5, 36.0, 41.0)
    assert fit["x0"] == pytest.approx(38.47, abs=0.01) and fit["fwhm"] == pytest.approx(0.6, rel=0.03)
    assert fit["r2"] > 0.95 and fit["fwhm_unc"] > 0
    assert fit["window"][0] >= 35.0 and fit["window"][1] <= 42.0  # clipped to the data


def test_fit_peak_needs_enough_points():
    with pytest.raises(ValueError, match="too few points"):
        tool.fit_peak(tool.Pattern(np.arange(35.0, 42.0, 0.02), np.ones(350)), 38.0, 38.0, 38.05)


def test_fit_all_windows_stop_halfway_between_peaks():
    p = pattern(noise=False)
    fits = tool.fit_all(p, [38.5, 44.7])
    assert fits[0]["window"][1] <= 41.6 + 0.02 and fits[1]["window"][0] >= 41.6 - 0.02


# ---------- detecting peaks ----------
def test_detects_the_true_peaks_and_no_noise_peaks():
    found = tool.detect_peaks(pattern())
    expected = [tt for tt, *_ in reflections()]
    assert len(found) == len(expected) and found == pytest.approx(expected, abs=0.25)  # detection is approximate
    fitted = [f["x0"] for f in tool.fit_all(pattern(), found)]
    assert fitted == pytest.approx(expected, abs=0.1)  # the fit refines the positions


def test_detection_controls():
    p = pattern()
    assert len(tool.detect_peaks(p, prominence=0.3)) == 2  # only the two strongest reflections pass 30% of the range
    assert tool.detect_peaks(tool.Pattern(p.x, np.full_like(p.x, 5.0))) == []
    assert len(tool.detect_peaks(p, min_separation=10.0)) < 6


def test_auto_window_is_odd_and_scales_with_the_step():
    assert tool.auto_window(tool.Pattern(np.arange(0, 10, 0.02), np.zeros(500))) == 13
    assert tool.auto_window(tool.Pattern(np.arange(0, 10, 0.05), np.zeros(200))) == 5
    assert tool.auto_window(tool.Pattern(np.arange(0, 10, 0.01), np.zeros(1000))) % 2 == 1


# ---------- reading ----------
def test_read_pattern_variants(tmp_path):
    rows = "".join(f"{20 + 0.1 * i:.1f},{100 + i}\n" for i in range(100))
    p = write(tmp_path, "Instrument: demo\ntwo theta,counts\n" + rows)
    pat = tool.read_pattern(p)
    assert pat.x.size == 100 and pat.y[0] == 100
    p = write(tmp_path, "tt,a,b\n" + "".join(f"{20 + 0.1 * i:.1f},{i},{2 * i}\n" for i in range(100)))
    assert tool.read_pattern(p, "b").y[-1] == 198 and tool.read_pattern(p, "1").y[-1] == 99
    with pytest.raises(ValueError, match="column 'z' not found"):
        tool.read_pattern(p, "z")
    with pytest.raises(ValueError, match="no numeric data"):
        tool.read_pattern(write(tmp_path, "text\n"))
    with pytest.raises(ValueError, match="at least 50"):
        tool.read_pattern(write(tmp_path, "1,1\n2,2\n3,3\n"))
    with pytest.raises(ValueError, match="cannot read"):
        tool.read_pattern(tmp_path / "missing.csv")


# ---------- warnings ----------
def test_size_warnings():
    big = tool.analyze_peak(40.0, 0.07, LAMBDA, k=K, inst=0.0, method="none")  # about 120 nm
    assert any("unreliable above" in n for n in tool.size_warnings([big], 0.0, None))
    near = tool.analyze_peak(40.0, 0.11, LAMBDA, k=K, inst=0.10, method="gaussian")
    assert any("times the instrumental width" in n for n in tool.size_warnings([near], 0.10, None))
    none = tool.analyze_peak(40.0, 0.08, LAMBDA, k=K, inst=0.10, method="gaussian")
    assert any("not wider than the instrumental" in n for n in tool.size_warnings([none], 0.10, None))
    assert tool.size_warnings([tool.analyze_peak(40.0, 0.8, LAMBDA, k=K, inst=0.1, method="gaussian")], 0.1, None) == []
    assert any("poor" in n for n in tool.size_warnings([], 0.0, [{"x0": 30.0, "r2": 0.1}]))


# ---------- output ----------
@pytest.mark.parametrize("name, magic", [("p.png", b"\x89PNG"), ("p.pdf", b"%PDF"), ("p.svg", b"<?xm")])
def test_plot_formats_with_and_without_williamson_hall(tmp_path, name, magic):
    p = pattern()
    fits = tool.fit_all(p, tool.detect_peaks(p))
    rows = [tool.analyze_peak(f["x0"], f["fwhm"], LAMBDA, k=K, inst=0.08, method="gaussian") for f in fits]
    tool.plot_xrd(p, fits, rows, tool.williamson_hall(rows, LAMBDA), tmp_path / name, title="t", dpi=60, figsize=None)
    assert (tmp_path / name).read_bytes()[:4] == magic
    tool.plot_xrd(p, fits, rows, None, tmp_path / ("n" + name), title=None, dpi=60, figsize=(5, 3))
    assert (tmp_path / ("n" + name)).is_file()


def test_table_and_csv(tmp_path):
    rows = [tool.analyze_peak(38.5, 0.5, LAMBDA, k=K, inst=0.1, method="gaussian", fwhm_unc=0.01),
            tool.analyze_peak(44.7, 0.05, LAMBDA, k=K, inst=0.1, method="gaussian")]
    assert "n/a" in tool.result_table(rows)
    tool.write_csv(tmp_path / "r.csv", rows, None)
    parsed = list(csv.DictReader((tmp_path / "r.csv").open()))
    assert list(parsed[0]) == tool.RESULT_COLUMNS and parsed[1]["size_nm"] == "" and parsed[0]["size_unc_nm"] != ""


# ---------- CLI ----------
def test_cli_example_pattern(tmp_path, monkeypatch, capsys):
    out_csv = tmp_path / "r.csv"
    assert run(monkeypatch, EXAMPLE, "--inst-fwhm", 0.08, "--williamson-hall", "-o", tmp_path / "p.png",
               "--csv", out_csv, "--dpi", 60) == 0
    out = capsys.readouterr().out
    assert "Scherrer size over 6 peaks" in out and "Williamson-Hall (6 peaks)" in out
    rows = list(csv.DictReader(out_csv.open()))
    assert [round(float(r["two_theta_deg"]), 1) for r in rows] == pytest.approx([38.5, 44.7, 65.1, 78.2, 82.4, 99.1], abs=0.1)
    wh_size = float(out.split("size ")[-1].split(" +/-")[0])
    assert wh_size == pytest.approx(20.0, abs=3.0)
    assert (tmp_path / "p.png").read_bytes()[:4] == b"\x89PNG"


def test_cli_options(tmp_path, monkeypatch, capsys):
    assert run(monkeypatch, EXAMPLE, "--source", "cu-ka", "--k", 0.94, "--correction", "lorentzian", "--inst-fwhm", 0.08,
               "--range", 30, 70, "--smooth", 0, "--prominence", 0.1, "--min-separation", 1) == 0
    assert "wavelength 1.54184" in capsys.readouterr().out
    assert run(monkeypatch, EXAMPLE, "--wavelength", 1.7889, "--peak-positions", 38.5, 44.7) == 0
    assert "wavelength 1.7889" in capsys.readouterr().out


def test_cli_direct_peaks_and_williamson_hall(monkeypatch, capsys):
    args = []
    for tt, beta, observed, _ in reflections():
        args += ["--peak", round(tt, 3), round(observed, 5)]
    assert run(monkeypatch, *args, "--inst-fwhm", 0.08, "--williamson-hall") == 0
    out = capsys.readouterr().out
    assert "Williamson-Hall (6 peaks): size 20.0" in out
    assert run(monkeypatch, "--peak", 38.5, 0.45, "--peak", 44.7, 0.47) == 0


def test_cli_errors(tmp_path, monkeypatch, capsys):
    assert run(monkeypatch) == 2  # nothing given
    assert run(monkeypatch, EXAMPLE, "--peak", 38.5, 0.4) == 2  # both
    assert run(monkeypatch, "--peak", 38.5, 0.4, "-o", tmp_path / "x.png") == 2
    assert run(monkeypatch, EXAMPLE, "-o", tmp_path / "x.jpg") == 2
    assert run(monkeypatch, EXAMPLE, "--range", 80, 30) == 2
    assert run(monkeypatch, EXAMPLE, "--range", 30, 30.2) == 2
    assert run(monkeypatch, EXAMPLE, "--smooth", 4) == 2
    assert run(monkeypatch, EXAMPLE, "--k", 0) == 2
    assert run(monkeypatch, "--peak", 200, 0.4) == 2
    assert run(monkeypatch, "--peak", 38.5, 0.4, "--williamson-hall") == 2  # needs 3 peaks
    flat = write(tmp_path, "".join(f"{20 + 0.1 * i:.1f},50\n" for i in range(100)), "flat.csv")
    assert run(monkeypatch, flat) == 2
    err = capsys.readouterr().err
    for text in ("either a pattern file", "needs a pattern file", "must end in", "LOW below HIGH", "fewer than 50 points",
                 "odd window", "must be positive", "between 0 and 180", "at least 3 peaks", "no peaks found"):
        assert text in err, text


def test_wide_window_for_a_very_broad_peak(tmp_path, monkeypatch, capsys):
    x = np.arange(10.0, 60.0, 0.05)
    y = 50 + 400 * tool.pseudo_voigt(x, 1.0, 26.0, 6.0, 0.5, 0.0, 0.0)
    src = write(tmp_path, "".join(f"{a:.2f},{b:.3f}\n" for a, b in zip(x, np.random.default_rng(4).poisson(y))))
    assert run(monkeypatch, src, "--window", 14, "--smooth", 0, "--peak-positions", 26.0) == 0
    out = capsys.readouterr().out
    row = next(ln.split() for ln in out.splitlines() if ln.split() and ln.split()[0].replace(".", "", 1).isdigit())
    assert float(row[2]) == pytest.approx(6.0, rel=0.06) and float(row[4]) < 3.0  # size of a few nm
    assert run(monkeypatch, src, "--window", 0) == 2
    assert "--window must be positive" in capsys.readouterr().err
