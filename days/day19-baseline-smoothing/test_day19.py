import csv
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import tool

EXAMPLE = Path(__file__).parent / "example"
SPECTRUM, TRUTH = EXAMPLE / "synthetic_spectrum.csv", EXAMPLE / "synthetic_baseline_truth.csv"
X = np.arange(400.0, 1801.0, 1.0)
PEAKS = [(520, 22, 0.55), (760, 18, 0.30), (1004, 14, 1.00), (1250, 28, 0.45), (1450, 35, 0.70), (1620, 24, 0.50)]


def peaks_signal(x=X, peaks=PEAKS):
    return sum(h * np.exp(-0.5 * ((x - c) / w) ** 2) for c, w, h in peaks)


def curved_baseline(x=X):
    u = (x - 400.0) / 1400.0
    return 0.15 + 0.55 * u - 0.35 * u**2 + 0.30 * np.exp(-0.5 * ((x - 1200.0) / 380.0) ** 2)


def spectrum(noise=0.006, seed=0):
    y = peaks_signal() + curved_baseline() + np.random.default_rng(seed).normal(0, noise, X.size)
    return tool.Spectrum("s", X.copy(), y)


def args(**kw):
    base = {"smooth": "none", "window": 11, "order": 3, "sigma": None, "smooth_lambda": 100.0, "baseline": "none",
            "baseline_from": "smoothed", "degree": 4, "baseline_lambda": 1e5, "asymmetry": 0.01, "width": None,
            "no_lls": False, "anchors": None, "anchor_width": 10.0, "anchor_kind": "linear"}
    base.update(kw)
    return SimpleNamespace(**base)


def run(monkeypatch, *argv):
    monkeypatch.setattr("sys.argv", ["tool.py", *map(str, argv)])
    return tool.main()


def write(tmp_path, text, name="s.csv"):
    p = tmp_path / name
    p.write_text(text)
    return p


# ---------- reading and noise ----------
def test_read_spectrum_variants_and_errors(tmp_path):
    rows = "".join(f"{1000 - 10 * i},{i / 10}\n" for i in range(30))  # descending, like many FTIR exports
    s = tool.read_spectrum(write(tmp_path, "Instrument: demo\ncm-1,abs\n" + rows))
    assert s.x[0] < s.x[-1] and s.name == "abs" and s.y.size == 30
    wide = write(tmp_path, "x,a,b\n" + "".join(f"{i},{i},{2 * i}\n" for i in range(30)))
    assert tool.read_spectrum(wide, "b").y[-1] == 58 and tool.read_spectrum(wide, "1").y[-1] == 29
    with pytest.raises(ValueError, match="column 'z' not found"):
        tool.read_spectrum(wide, "z")
    with pytest.raises(ValueError, match="no numeric data"):
        tool.read_spectrum(write(tmp_path, "text\n"))
    with pytest.raises(ValueError, match="only 3 data points"):
        tool.read_spectrum(write(tmp_path, "1,1\n2,2\n3,3\n"))
    with pytest.raises(ValueError, match="cannot read"):
        tool.read_spectrum(tmp_path / "missing.csv")


def test_noise_sigma_ignores_baseline_and_peaks():
    white = np.random.default_rng(1).normal(0, 0.05, 5000)
    assert tool.noise_sigma(white) == pytest.approx(0.05, rel=0.05)
    assert tool.noise_sigma(white + 10 * np.sin(np.linspace(0, 3, 5000)) + 3 * np.exp(-0.5 * ((np.arange(5000) - 2500) / 200) ** 2)) \
        == pytest.approx(0.05, rel=0.08)


def test_step_and_uniformity():
    assert tool.step_of(X) == 1.0 and tool.spacing_uniformity(X) == 0.0
    assert tool.spacing_uniformity(np.array([0.0, 1.0, 2.0, 4.0, 5.0])) == pytest.approx(1.0)


# ---------- smoothing ----------
def test_penalty_matrix_ignores_straight_lines():
    p = tool.penalty_matrix(8, 2)
    line = 3.0 + 2.0 * np.arange(8)
    assert p.shape == (8, 8) and np.allclose(p @ line, 0) and np.allclose(p @ np.ones(8), 0)
    assert np.linalg.norm(p @ np.arange(8.0) ** 2) > 1  # a curve is penalized


def test_whittaker_limits_and_errors():
    y = np.random.default_rng(2).normal(0, 1, 50)
    assert tool.whittaker_smooth(y, 0.0) == pytest.approx(y)
    line = 1.0 + 0.5 * np.arange(50)
    assert tool.whittaker_smooth(line, 1e6) == pytest.approx(line, abs=1e-6)  # a straight line is not penalized
    assert np.std(tool.whittaker_smooth(y, 1e4)) < 0.2 * np.std(y)
    with pytest.raises(ValueError, match="must not be negative"):
        tool.whittaker_smooth(y, -1)
    with pytest.raises(ValueError, match="at least 4 points"):
        tool.whittaker_smooth(np.ones(3), 1.0)


def test_savgol_keeps_a_polynomial_and_validates():
    x = np.linspace(0, 1, 101)
    parabola = 2 + 3 * x - 4 * x**2
    assert tool.smooth(parabola, "savgol", window=11, order=3) == pytest.approx(parabola, abs=1e-9)
    for window in (10, 1):
        with pytest.raises(ValueError, match="odd and at least 3"):
            tool.smooth(parabola, "savgol", window=window)
    with pytest.raises(ValueError, match="larger than --order"):
        tool.smooth(parabola, "savgol", window=3, order=3)
    with pytest.raises(ValueError, match="not longer than the data"):
        tool.smooth(parabola, "savgol", window=201)


def test_moving_average_gaussian_median_basics():
    flat = np.full(60, 4.0)
    for method in ("moving-average", "gaussian", "median", "whittaker", "savgol"):
        assert tool.smooth(flat, method) == pytest.approx(flat)
    spike = np.zeros(61)
    spike[30] = 10.0
    assert tool.smooth(spike, "median", window=5)[30] == 0.0  # the spike is removed
    assert tool.smooth(spike, "moving-average", window=5)[30] == pytest.approx(2.0)
    area = tool.smooth(spike, "gaussian", sigma_points=3.0).sum()
    assert area == pytest.approx(10.0, rel=1e-6)  # smoothing keeps the area
    with pytest.raises(ValueError, match="--sigma must be positive"):
        tool.smooth(spike, "gaussian", sigma_points=0)
    with pytest.raises(ValueError, match="between 1 and the number"):
        tool.smooth(spike, "moving-average", window=500)
    with pytest.raises(ValueError, match="unknown smoothing"):
        tool.smooth(spike, "magic")
    assert tool.smooth(spike, "none").tolist() == spike.tolist()


@pytest.mark.parametrize("method, kw", [("savgol", {"window": 15}), ("moving-average", {"window": 9}),
                                        ("gaussian", {"sigma_points": 3.0}), ("median", {"window": 9}),
                                        ("whittaker", {"lam": 50.0})])
def test_every_smoother_reduces_noise_and_keeps_the_peaks(method, kw):
    s = spectrum(noise=0.02, seed=3)
    out = tool.smooth(s.y, method, **kw)
    assert tool.noise_sigma(out) < 0.5 * tool.noise_sigma(s.y)
    assert out.max() > 0.93 * s.y.max() - 0.05  # the tallest peak survives a sensible window


def test_a_wide_window_flattens_the_peaks():
    s = spectrum()
    assert tool.smooth(s.y, "moving-average", window=9).max() > 0.97 * s.y.max()
    assert tool.smooth(s.y, "moving-average", window=61).max() < 0.85 * s.y.max()  # measured: 77% of the true height


# ---------- baselines ----------
def rms(a, b):
    return float(np.sqrt(np.mean((a - b) ** 2)))


def test_modpoly_recovers_a_polynomial_exactly_without_peaks():
    base = 1 + 0.002 * X - 1.5e-6 * X**2
    fit, n = tool.modpoly(X, base, degree=2)
    assert fit == pytest.approx(base, abs=1e-8) and n <= 3


def test_modpoly_stays_under_the_peaks():
    s = spectrum()
    fit, n = tool.modpoly(X, s.y, degree=5)
    assert rms(fit, curved_baseline()) < 0.035 and 1 <= n <= 100
    assert np.mean(s.y >= fit - 0.03) > 0.99  # the fit lies at or below the data
    with pytest.raises(ValueError, match="--degree"):
        tool.modpoly(X, s.y, degree=-1)


def test_asls_recovers_the_baseline_and_validates():
    s = spectrum()
    z, n = tool.asls(s.y, 1e7, 0.01)
    assert rms(z, curved_baseline()) < 0.025 and n < 20
    assert np.mean(s.y > z) > 0.5  # the baseline sits under the data: most points are above it (measured: 78%)
    with pytest.raises(ValueError, match="between 0 and 1"):
        tool.asls(s.y, 1e5, 1.5)
    with pytest.raises(ValueError, match="must be positive"):
        tool.asls(s.y, 0.0, 0.01)


def test_asls_stiffness_matters():
    s = spectrum()
    stiff = rms(tool.asls(s.y, 1e7, 0.01)[0], curved_baseline())
    loose = rms(tool.asls(s.y, 1e5, 0.01)[0], curved_baseline())
    assert loose > 2 * stiff  # a baseline that is too flexible starts to follow the peaks


def test_snip_basics():
    assert tool.snip(np.full(200, 3.0), 10) == pytest.approx(3.0)
    line = 2.0 + 0.01 * X
    for lls in (True, False):  # clipping never rises above the data and stays within 5% of a sloping line's total rise
        b = tool.snip(line, 50, lls=lls)
        assert np.all(b <= line + 1e-9) and np.max(line - b) < 0.05 * np.ptp(line)
    s = spectrum()
    base = tool.snip(s.y, 50)
    assert rms(base, curved_baseline()) < 0.03 and np.mean(s.y >= base - 0.03) > 0.99
    assert rms(tool.snip(s.y, 50, lls=False), curved_baseline()) < 0.03
    with pytest.raises(ValueError, match="at least 1 point"):
        tool.snip(s.y, 0)


def test_snip_window_too_small_follows_the_peaks():
    s = spectrum()
    assert rms(tool.snip(s.y, 5), curved_baseline()) > 3 * rms(tool.snip(s.y, 50), curved_baseline())


def test_anchor_baseline_linear_cubic_and_errors():
    y = 1.0 + 0.001 * X
    base = tool.anchor_baseline(X, y, [450, 900, 1500], 10)
    between = (X >= 450) & (X <= 1500)
    assert base[between] == pytest.approx(y[between], abs=1e-9)  # exact between the anchors for a straight line
    assert base[0] == pytest.approx(base[np.argmin(np.abs(X - 450))])  # constant before the first anchor
    s = spectrum()
    cubic = tool.anchor_baseline(X, s.y, [420, 650, 880, 1130, 1360, 1540, 1790], 20, "cubic")
    assert rms(cubic, curved_baseline()) < 0.025
    with pytest.raises(ValueError, match="at least 2 positions"):
        tool.anchor_baseline(X, y, [500], 10)
    with pytest.raises(ValueError, match="4 for a cubic"):
        tool.anchor_baseline(X, y, [500, 900, 1200], 10, "cubic")
    with pytest.raises(ValueError, match="must be positive"):
        tool.anchor_baseline(X, y, [500, 900], 0)
    with pytest.raises(ValueError, match="no data within"):
        tool.anchor_baseline(X, y, [500, 5000], 10)


def test_estimate_baseline_needs_its_parameters():
    s = spectrum()
    with pytest.raises(ValueError, match="needs --width"):
        tool.estimate_baseline(X, s.y, "snip", args())
    with pytest.raises(ValueError, match="needs --anchors"):
        tool.estimate_baseline(X, s.y, "anchors", args())
    base, text = tool.estimate_baseline(X, s.y, "none", args())
    assert base.max() == 0 and text == "none"
    assert "degree 4" in tool.estimate_baseline(X, s.y, "poly", args())[1]
    with pytest.raises(ValueError, match="unknown baseline"):
        tool.estimate_baseline(X, s.y, "magic", args())


# ---------- the pipeline and its checks ----------
def test_process_subtracts_the_baseline_from_the_smoothed_spectrum():
    s = spectrum()
    r = tool.process(s, args(smooth="savgol", baseline="asls", baseline_lambda=1e7))
    assert r["corrected"] == pytest.approx(r["smoothed"] - r["baseline"])
    assert r["noise_smoothed"] < r["noise_raw"] and 0.97 < r["peak_retained"] <= 1.01
    raw_based = tool.process(s, args(smooth="savgol", baseline="asls", baseline_lambda=1e7, baseline_from="raw"))
    assert raw_based["baseline"] != pytest.approx(r["baseline"], abs=1e-9)


def test_peak_prominence_check_tells_good_baselines_from_bad_ones():
    s = spectrum()
    good = tool.process(s, args(baseline="asls", baseline_lambda=1e7))["prominence_kept"]
    bad = tool.process(s, args(baseline="asls", baseline_lambda=1e3, asymmetry=0.5))["prominence_kept"]
    assert good[1] == 6 and 0.95 < good[0] < 1.1
    assert bad[0] < 0.3
    assert tool.process(s, args())["prominence_kept"] is None
    assert tool.prominence_kept(np.zeros(100), np.zeros(100), 0.01) is None  # no peaks
    assert tool.prominence_kept(s.y, s.y, 0.0) is None


def test_baseline_error_measures_the_difference():
    truth = curved_baseline()
    err = tool.baseline_error(X, truth + 0.1, X, truth)
    assert err["rms"] == pytest.approx(0.1) and err["max"] == pytest.approx(0.1) and err["mean"] == pytest.approx(0.1)
    coarse = tool.baseline_error(X, truth, X[::50], truth[::50])  # the truth is interpolated onto the spectrum
    assert coarse["rms"] < 0.01


def test_warnings():
    s = spectrum()
    clean = tool.process(s, args(smooth="savgol", baseline="asls", baseline_lambda=1e7))
    assert tool.warnings_for(s, clean, args(smooth="savgol", baseline="asls")) == []
    wide = tool.process(s, args(smooth="moving-average", window=61, baseline="asls", baseline_lambda=1e7))
    assert any("lowered the highest peak" in n for n in tool.warnings_for(s, wide, args(smooth="moving-average", baseline="asls")))
    flexible = tool.process(s, args(baseline="asls", baseline_lambda=1e3, asymmetry=0.5))
    assert any("following the peaks" in n for n in tool.warnings_for(s, flexible, args(baseline="asls")))
    assert any("high polynomial degree" in n for n in tool.warnings_for(s, clean, args(baseline="poly", degree=12)))
    uneven = tool.Spectrum("u", np.concatenate([np.arange(0, 100, 1.0), np.arange(100, 300, 5.0)]), np.ones(140))
    assert any("not uniform" in n for n in tool.warnings_for(uneven, clean, args()))
    below = dict(clean, corrected=np.full(X.size, -1.0))
    assert any("below zero" in n for n in tool.warnings_for(s, below, args(baseline="asls")))


# ---------- output ----------
@pytest.mark.parametrize("name, magic", [("p.png", b"\x89PNG"), ("p.pdf", b"%PDF"), ("p.svg", b"<?xm")])
def test_plot_formats_and_options(tmp_path, name, magic):
    s = spectrum()
    r = tool.process(s, args(smooth="savgol", baseline="asls", baseline_lambda=1e7))
    tool.plot_result(s, r, tmp_path / name, truth=(X, curved_baseline()), smoothed=True, reverse_x=True, xlabel="x",
                     ylabel="y", title="t", figsize=None, dpi=60)
    assert (tmp_path / name).read_bytes()[:4] == magic
    tool.plot_result(s, r, tmp_path / ("n" + name), truth=None, smoothed=False, reverse_x=False, xlabel="x", ylabel="y",
                     title=None, figsize=(5, 4), dpi=60)
    assert (tmp_path / ("n" + name)).is_file()


def test_csv_columns_and_values(tmp_path):
    s = spectrum()
    r = tool.process(s, args(smooth="gaussian", sigma=3.0, baseline="snip", width=100))
    tool.write_csv(tmp_path / "c.csv", s, r)
    rows = list(csv.DictReader((tmp_path / "c.csv").open()))
    assert list(rows[0]) == tool.CSV_COLUMNS and len(rows) == X.size
    row = rows[700]
    assert float(row["corrected"]) == pytest.approx(float(row["smoothed"]) - float(row["baseline"]), abs=1e-6)


def test_corrected_only_csv(tmp_path):
    s = spectrum()
    r = tool.process(s, args(baseline="snip", width=100))
    tool.write_csv(tmp_path / "c.csv", s, r, corrected_only=True)
    rows = list(csv.DictReader((tmp_path / "c.csv").open()))
    assert list(rows[0]) == ["x", "corrected"] and len(rows) == X.size
    assert float(rows[700]["corrected"]) == pytest.approx(r["corrected"][700], abs=1e-6)


# ---------- CLI ----------
def test_cli_recommended_settings_recover_the_example_baseline(tmp_path, monkeypatch, capsys):
    out_csv = tmp_path / "c.csv"
    assert run(monkeypatch, SPECTRUM, "--smooth", "savgol", "--window", 11, "--baseline", "asls", "--baseline-lambda", 1e7,
               "--truth", TRUTH, "-o", tmp_path / "p.png", "--csv", out_csv, "--dpi", 60, "--reverse-x") == 0
    out = capsys.readouterr().out
    assert "highest peak kept at 99" in out and "102% of the peak prominence kept" in out
    error = float(out.split("RMS error ")[1].split()[0])
    assert error < 0.025
    assert (tmp_path / "p.png").read_bytes()[:4] == b"\x89PNG" and out_csv.is_file()
    assert run(monkeypatch, SPECTRUM, "--baseline", "poly", "--csv", tmp_path / "two.csv", "--corrected-only") == 0
    assert next(iter(csv.reader((tmp_path / "two.csv").open()))) == ["x", "corrected"]


@pytest.mark.parametrize("flags", [
    ["--baseline", "poly", "--degree", 5],
    ["--baseline", "asls", "--baseline-lambda", 1e7],
    ["--baseline", "snip", "--width", 100],
    ["--baseline", "snip", "--width", 100, "--no-lls"],
    ["--baseline", "anchors", "--anchors", 420, 650, 880, 1130, 1360, 1540, 1790, "--anchor-width", 20, "--anchor-kind", "cubic"],
])
def test_every_baseline_method_recovers_the_known_baseline(monkeypatch, capsys, flags):
    assert run(monkeypatch, SPECTRUM, *flags, "--truth", TRUTH) == 0
    out = capsys.readouterr().out
    assert float(out.split("RMS error ")[1].split()[0]) < 0.03
    assert float(out.split("peaks: ")[1].split("%")[0]) > 90


def test_cli_warns_about_bad_settings(monkeypatch, capsys):
    assert run(monkeypatch, SPECTRUM, "--baseline", "asls", "--baseline-lambda", 1e3, "--asymmetry", 0.5) == 0
    assert "following the peaks" in capsys.readouterr().err
    assert run(monkeypatch, SPECTRUM, "--smooth", "moving-average", "--window", 61) == 0
    assert "lowered the highest peak" in capsys.readouterr().err


def test_cli_other_smoothers_and_range(tmp_path, monkeypatch, capsys):
    for flags in (["--smooth", "gaussian", "--sigma", 3], ["--smooth", "median", "--window", 7],
                  ["--smooth", "whittaker", "--smooth-lambda", 50], ["--smooth", "moving-average", "--window", 5]):
        assert run(monkeypatch, SPECTRUM, *flags, "--xmin", 450, "--xmax", 1700) == 0
    assert "x 450 to 1700" in capsys.readouterr().out


def test_cli_errors(tmp_path, monkeypatch, capsys):
    assert run(monkeypatch, SPECTRUM, "-o", tmp_path / "x.jpg") == 2
    assert run(monkeypatch, tmp_path / "missing.csv") == 2
    assert run(monkeypatch, SPECTRUM, "--column", "nope") == 2
    assert run(monkeypatch, SPECTRUM, "--baseline", "snip") == 2
    assert run(monkeypatch, SPECTRUM, "--baseline", "anchors") == 2
    assert run(monkeypatch, SPECTRUM, "--smooth", "savgol", "--window", 10) == 2
    assert run(monkeypatch, SPECTRUM, "--xmin", 1000, "--xmax", 1003) == 2
    assert run(monkeypatch, SPECTRUM, "--baseline", "asls", "--asymmetry", 2) == 2
    assert run(monkeypatch, SPECTRUM, "--truth", tmp_path / "missing.csv", "--baseline", "poly") == 2
    err = capsys.readouterr().err
    for text in ("must end in", "column 'nope' not found", "needs --width", "needs --anchors", "odd and at least 3",
                 "fewer than 10 points", "between 0 and 1", "cannot read"):
        assert text in err, text
