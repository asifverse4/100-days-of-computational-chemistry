import csv
from pathlib import Path

import numpy as np
import pytest
import tool

EXAMPLE = Path(__file__).parent / "example" / "synthetic_cqd.csv"


def gaussians(centers, heights, sigma=15.0, lo=200.0, hi=700.0, step=1.0, noise=0.0, seed=0, name="s"):
    x = np.arange(lo, hi + step, step)
    y = sum(h * np.exp(-0.5 * ((x - c) / sigma) ** 2) for c, h in zip(centers, heights))
    if noise:
        y = y + np.random.default_rng(seed).normal(0, noise, x.size)
    return tool.Spectrum(name, x, y)


def run(monkeypatch, *argv):
    monkeypatch.setattr("sys.argv", ["tool.py", *map(str, argv)])
    return tool.main()


# ---------- peak detection ----------
def test_peak_positions_and_heights():
    peaks = tool.detect_peaks(gaussians([250, 400], [1.0, 0.5]))
    assert [round(p["wavelength_nm"]) for p in peaks] == [250, 400]
    assert [p["absorbance"] for p in peaks] == pytest.approx([1.0, 0.5], rel=0.01)


def test_fwhm_of_isolated_gaussian():
    sigma = 20.0
    (peak,) = tool.detect_peaks(gaussians([400], [1.0], sigma=sigma))
    assert peak["fwhm_nm"] == pytest.approx(2.3548 * sigma, rel=0.03)
    assert peak["prominence"] == pytest.approx(1.0, rel=0.01)


def test_prominence_threshold():
    spec = gaussians([250, 450], [1.0, 0.03])
    assert len(tool.detect_peaks(spec)) == 1  # default threshold is 5% of the range
    assert len(tool.detect_peaks(spec, prominence=0.01)) == 2


def test_min_separation_keeps_the_stronger_peak():
    spec = gaussians([300, 306], [1.0, 0.8], sigma=2.0)
    close = tool.detect_peaks(spec, prominence=0.01, min_separation_nm=10)
    far = tool.detect_peaks(spec, prominence=0.01, min_separation_nm=2)
    assert len(close) == 1 and round(close[0]["wavelength_nm"]) == 300
    assert len(far) == 2


def test_no_peaks_in_flat_or_monotonic_data():
    x = np.arange(200, 300.0)
    assert tool.detect_peaks(tool.Spectrum("flat", x, np.full_like(x, 0.2))) == []
    assert tool.detect_peaks(tool.Spectrum("ramp", x, x / 300)) == []


def test_smoothing_removes_noise_peaks():
    spec = gaussians([300], [1.0], noise=0.02, seed=3)
    smoothed = tool.smooth(spec, 21)
    peaks = tool.detect_peaks(smoothed, prominence=0.1)
    assert len(peaks) == 1 and abs(peaks[0]["wavelength_nm"] - 300) <= 2
    assert np.std(smoothed.y - gaussians([300], [1.0]).y) < np.std(spec.y - gaussians([300], [1.0]).y)


# ---------- processing ----------
def test_normalize_and_errors():
    spec = tool.normalize(gaussians([300], [4.0]))
    assert spec.y.max() == pytest.approx(1.0)
    with pytest.raises(ValueError, match="cannot normalize"):
        tool.normalize(tool.Spectrum("neg", np.arange(10.0), -np.ones(10)))


def test_zero_at_uses_nearest_wavelength():
    spec = tool.Spectrum("s", np.arange(700.0, 801.0, 10), np.linspace(0.1, 0.3, 11))
    out = tool.zero_at(spec, 803)
    assert out.y[-1] == 0.0 and out.y[0] == pytest.approx(-0.2)


def test_crop_and_too_few_points():
    spec = gaussians([300], [1.0])
    out = tool.crop(spec, 250, 350)
    assert out.x.min() == 250 and out.x.max() == 350
    with pytest.raises(ValueError, match="fewer than 5 points"):
        tool.crop(spec, 300, 302)


def test_percent_t_conversion():
    assert tool.percent_t_to_absorbance(np.array([100.0, 10.0, 50.0])) == pytest.approx([0, 1, 0.30103], abs=1e-5)
    with pytest.raises(ValueError, match="above 0"):
        tool.percent_t_to_absorbance(np.array([50.0, 0.0]))


@pytest.mark.parametrize("window", [4, 3, 6])
def test_smooth_rejects_bad_windows(window):
    with pytest.raises(ValueError, match="odd window"):
        tool.smooth(gaussians([300], [1.0]), window)


def test_smooth_window_longer_than_data():
    with pytest.raises(ValueError, match="longer than the data"):
        tool.smooth(tool.Spectrum("s", np.arange(7.0), np.arange(7.0)), 9)


# ---------- reading files ----------
ROWS = "".join(f"{w},{a}\n" for w, a in [(300, 0.1), (310, 0.2), (320, 0.5), (330, 0.2), (340, 0.1), (350, 0.05)])


def write(tmp_path, text, name="s.csv"):
    p = tmp_path / name
    p.write_text(text)
    return p


def test_read_headerless_two_columns(tmp_path):
    (spec,) = tool.read_spectra(write(tmp_path, ROWS, "sample1.csv"))
    assert spec.name == "sample1" and spec.x.tolist() == [300, 310, 320, 330, 340, 350]
    assert spec.y[2] == 0.5


def test_read_header_metadata_and_footer(tmp_path):
    text = "Instrument: demo\nDate: today\nWavelength (nm),Abs\n" + ROWS + "End of file\n"
    (spec,) = tool.read_spectra(write(tmp_path, text))
    assert spec.name == "Abs" and spec.x.size == 6


def test_read_tab_separated_multi_column(tmp_path):
    text = "nm\tA\tB\n" + "".join(f"{300 + 10 * i}\t{i / 10}\t{i / 5}\n" for i in range(6))
    a, b = tool.read_spectra(write(tmp_path, text))
    assert (a.name, b.name) == ("A", "B") and b.y[5] == pytest.approx(1.0)


def test_read_semicolon_with_decimal_comma(tmp_path):
    text = "wl;abs\n" + "".join(f"{300 + 10 * i};0,{i}5\n" for i in range(6))
    (spec,) = tool.read_spectra(write(tmp_path, text))
    assert spec.y.tolist() == pytest.approx([0.05, 0.15, 0.25, 0.35, 0.45, 0.55])


def test_read_whitespace_separated_descending_and_duplicates(tmp_path):
    rows = [(350, 0.05), (340, 0.1), (340, 0.1), (330, 0.2), (320, 0.5), (310, 0.2), (300, 0.1)]
    (spec,) = tool.read_spectra(write(tmp_path, "".join(f"{w} {a}\n" for w, a in rows)))
    assert spec.x.tolist() == [300, 310, 320, 330, 340, 350]  # sorted, duplicate removed


def test_read_errors(tmp_path):
    with pytest.raises(ValueError, match="no numeric data"):
        tool.read_spectra(write(tmp_path, "just text\nmore text\n"))
    with pytest.raises(ValueError, match="only 3 data points"):
        tool.read_spectra(write(tmp_path, "1,1\n2,2\n3,3\n"))
    with pytest.raises(ValueError, match="cannot read"):
        tool.read_spectra(tmp_path / "missing.csv")


def test_unique_names_across_files(tmp_path):
    a = tool.read_spectra(write(tmp_path, "nm,Abs\n" + ROWS, "a.csv"))
    b = tool.read_spectra(write(tmp_path, "nm,Abs\n" + ROWS, "b.csv"))
    names = [s.name for s in tool.unique_names([(tmp_path / "a.csv", a), (tmp_path / "b.csv", b)])]
    assert names == ["a Abs", "b Abs"]


# ---------- output ----------
@pytest.mark.parametrize("name, magic", [("p.png", b"\x89PNG"), ("p.pdf", b"%PDF"), ("p.svg", b"<?xml")])
def test_plot_formats(tmp_path, name, magic):
    spec = gaussians([300], [1.0])
    tool.plot([spec], tool.detect_peaks(spec), tmp_path / name, dpi=60)
    assert (tmp_path / name).read_bytes()[:4] == magic[:4]


def test_peak_table_and_csv(tmp_path):
    peaks = tool.detect_peaks(gaussians([300], [1.0], name="sampleX"))
    assert "sampleX" in tool.peak_table(peaks) and "300.0" in tool.peak_table(peaks)
    tool.write_peaks_csv(tmp_path / "p.csv", peaks)
    rows = list(csv.DictReader((tmp_path / "p.csv").open()))
    assert list(rows[0]) == tool.PEAK_COLUMNS and float(rows[0]["wavelength_nm"]) == pytest.approx(300, abs=1)


# ---------- CLI ----------
def test_cli_example(tmp_path, monkeypatch, capsys):
    out, peaks_csv = tmp_path / "plot.png", tmp_path / "peaks.csv"
    assert run(monkeypatch, EXAMPLE, "-o", out, "--peaks-csv", peaks_csv, "--dpi", 80) == 0
    assert out.read_bytes()[:4] == b"\x89PNG"
    assert "[ok] 2 spectra" in capsys.readouterr().out
    rows = list(csv.DictReader(peaks_csv.open()))
    assert [r["sample"] for r in rows] == ["sample_A", "sample_A", "sample_B", "sample_B"]
    found = [float(r["wavelength_nm"]) for r in rows]
    assert found[0] == pytest.approx(240, abs=5) and found[1] == pytest.approx(340, abs=8)


def test_cli_processing_options(tmp_path, monkeypatch):
    assert run(monkeypatch, EXAMPLE, EXAMPLE, "-o", tmp_path / "p.png", "--normalize", "--stack", 0.3, "--smooth", 11,
               "--zero-at", 600, "--xmin", 230, "--xmax", 600, "--no-labels", "--title", "t", "--dpi", 60) == 0


def test_cli_percent_t(tmp_path, monkeypatch, capsys):
    spec = gaussians([300], [0.8])
    t_values = 100 * 10 ** (-spec.y)
    src = write(tmp_path, "nm,T\n" + "".join(f"{x},{t}\n" for x, t in zip(spec.x, t_values)))
    assert run(monkeypatch, src, "-o", tmp_path / "p.png", "--percent-t", "--dpi", 60) == 0
    assert "300.0" in capsys.readouterr().out


def test_cli_no_peaks_is_not_an_error(tmp_path, monkeypatch, capsys):
    src = write(tmp_path, "".join(f"{300 + i},{0.1 + i / 1000}\n" for i in range(50)))
    assert run(monkeypatch, src, "-o", tmp_path / "p.png", "--dpi", 60) == 0
    assert "no peaks found" in capsys.readouterr().out


def test_cli_errors(tmp_path, monkeypatch, capsys):
    assert run(monkeypatch, EXAMPLE, "-o", tmp_path / "x.jpg") == 2
    assert run(monkeypatch, tmp_path / "missing.csv", "-o", tmp_path / "x.png") == 2
    assert run(monkeypatch, EXAMPLE, "-o", tmp_path / "x.png", "--smooth", 4) == 2
    assert run(monkeypatch, EXAMPLE, "-o", tmp_path / "x.png", "--xmin", 599, "--xmax", 600) == 2
    err = capsys.readouterr().err
    assert "must end in" in err and "odd window" in err and "fewer than 5 points" in err
