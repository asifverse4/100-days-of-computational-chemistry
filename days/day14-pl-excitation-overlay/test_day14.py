import csv
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import tool

EXAMPLE = Path(__file__).parent / "example" / "synthetic_pl.csv"
WL = np.arange(320.0, 701.0, 1.0)


def band(center, sigma=38.0, amplitude=1.0, noise=0.0, seed=0, x=WL):
    y = amplitude * np.exp(-0.5 * ((x - center) / sigma) ** 2)
    return y + np.random.default_rng(seed).normal(0, noise, x.size) if noise else y


def spec(ex, y=None, x=WL):
    return tool.Spectrum(float(ex), x.copy(), band(430 + 0.45 * (ex - 300)) if y is None else y, f"{ex:g} nm")


def args(**kw):
    base = {"xmin": None, "xmax": None, "smooth": 0, "mask_scatter": 0.0, "normalize": False}
    base.update(kw)
    return SimpleNamespace(**base)


def run(monkeypatch, *argv):
    monkeypatch.setattr("sys.argv", ["tool.py", *map(str, argv)])
    return tool.main()


def write(tmp_path, text, name="s.csv"):
    p = tmp_path / name
    p.write_text(text)
    return p


# ---------- reading and excitation wavelengths ----------
def rows(n=12, cols=1):
    return "".join(f"{400 + 10 * i}," + ",".join(f"{0.1 * (c + 1) + i / 100}" for c in range(cols)) + "\n"
                   for i in range(n))


@pytest.mark.parametrize("header", ["340", "340 nm", "Ex 340", "ex=340nm", "PL (lambda_ex 340)"])
def test_excitation_from_header_variants(tmp_path, header):
    raws = tool.read_raw(write(tmp_path, f"emission,{header}\n" + rows()))
    assert tool.assign_excitation(raws, None)[0].ex == 340.0


def test_wide_file_sorted_by_excitation(tmp_path):
    raws = tool.read_raw(write(tmp_path, "nm,ex380,ex340,ex360\n" + rows(cols=3)))
    assert [s.ex for s in tool.assign_excitation(raws, None)] == [340, 360, 380]


def test_excitation_from_file_name_when_header_has_no_number(tmp_path):
    raws = tool.read_raw(write(tmp_path, "nm,Intensity\n" + rows(), "ex365.csv"))
    assert tool.assign_excitation(raws, None)[0].ex == 365.0
    raws = tool.read_raw(write(tmp_path, rows(), "sample_ex420nm.txt"))  # headerless
    assert tool.assign_excitation(raws, None)[0].ex == 420.0


def test_ex_override_and_errors(tmp_path):
    raws = tool.read_raw(write(tmp_path, "nm,a,b\n" + rows(cols=2), "data.csv"))
    assert [s.ex for s in tool.assign_excitation(raws, [300, 320])] == [300, 320]
    with pytest.raises(ValueError, match="gives 3 wavelengths but there are 2"):
        tool.assign_excitation(raws, [300, 320, 340])
    with pytest.raises(ValueError, match="same excitation"):
        tool.assign_excitation(raws, [300, 300])
    with pytest.raises(ValueError, match="cannot find the excitation"):
        tool.assign_excitation(raws, None)


def test_read_errors(tmp_path):
    with pytest.raises(ValueError, match="no numeric data"):
        tool.read_raw(write(tmp_path, "just text\n"))
    with pytest.raises(ValueError, match="only 3 data points"):
        tool.read_raw(write(tmp_path, "1,1\n2,2\n3,3\n"))
    with pytest.raises(ValueError, match="cannot read"):
        tool.read_raw(tmp_path / "missing.csv")


# ---------- processing ----------
def test_crop_smooth_normalize():
    s = spec(340)
    assert tool.crop(s, 400, 500).x.min() == 400
    with pytest.raises(ValueError, match="fewer than 5 points"):
        tool.crop(s, 400, 402)
    with pytest.raises(ValueError, match="odd window"):
        tool.smooth(s, 10)
    with pytest.raises(ValueError, match="longer than the data"):
        tool.smooth(tool.Spectrum(1, np.arange(7.0), np.arange(7.0), "s"), 9)
    assert tool.normalize(spec(340, y=band(450, amplitude=4.0))).y.max() == pytest.approx(1.0)
    with pytest.raises(ValueError, match="cannot normalize"):
        tool.normalize(spec(340, y=-np.ones(WL.size)))


def test_mask_scatter_blanks_both_orders():
    masked = tool.mask_scatter(spec(340), 15)
    nan = np.isnan(masked.y)
    assert nan[(WL >= 325) & (WL <= 355)].all() and nan[(WL >= 665) & (WL <= 695)].all()
    assert not nan[(WL >= 400) & (WL <= 600)].any()
    assert tool.normalize(masked).y[~nan].max() == pytest.approx(1.0)  # NaN-safe


# ---------- finding the emission maximum ----------
@pytest.mark.parametrize("center", [431.3, 457.0, 489.7])
def test_noisy_broad_band_maximum_is_accurate(center):
    s = spec(300, y=band(center, noise=0.005, seed=3))
    assert tool.find_maximum(s)["emission_max_nm"] == pytest.approx(center, abs=0.6)


def test_refine_falls_back_for_a_very_narrow_peak():
    x = np.arange(400.0, 420.0, 1.0)
    y = np.zeros_like(x)
    y[10] = 1.0
    assert tool.refine_peak(x, y, 10) == (410.0, 1.0)


def test_fwhm_of_a_gaussian():
    peak = tool.find_maximum(spec(300, y=band(500, sigma=30)))
    assert peak["fwhm_nm"] == pytest.approx(2.3548 * 30, rel=0.02)
    assert peak["intensity"] == pytest.approx(1.0, abs=0.01)


def test_fwhm_is_none_when_the_band_is_cut_off():
    assert tool.find_maximum(spec(300, y=band(690, sigma=40)))["fwhm_nm"] is None


def test_search_offset_ignores_rayleigh_scatter():
    y = band(450, amplitude=0.3) + 0.9 * np.exp(-0.5 * ((WL - 340) / 4) ** 2)  # strong scatter at the excitation
    assert tool.find_maximum(spec(340, y=y), search_offset=10)["emission_max_nm"] == pytest.approx(450, abs=1)
    assert tool.find_maximum(spec(340, y=y), search_offset=-50)["emission_max_nm"] == pytest.approx(340, abs=1)


def test_second_order_scatter_needs_masking():
    y = band(430, amplitude=0.3) + 0.35 * np.exp(-0.5 * ((WL - 600) / 6) ** 2)  # ex 300 nm, second order at 600
    raw = tool.find_maximum(spec(300, y=y))
    assert raw["emission_max_nm"] == pytest.approx(600, abs=1)
    masked = tool.find_maximum(tool.mask_scatter(spec(300, y=y), 15))
    assert masked["emission_max_nm"] == pytest.approx(430, abs=1)


def test_no_maximum_when_nothing_is_above_the_excitation():
    assert tool.find_maximum(spec(695)) is None


# ---------- the shift fit ----------
def test_fit_shift_recovers_the_slope():
    peaks = [{"excitation_nm": ex, "emission_max_nm": 430 + 0.45 * (ex - 300)} for ex in (300, 340, 380, 420)]
    fit = tool.fit_shift(peaks)
    assert fit["slope"] == pytest.approx(0.45) and fit["r2"] == pytest.approx(1.0)
    assert fit["shift_nm"] == pytest.approx(54.0) and fit["ex_range_nm"] == 120.0
    assert tool.fit_shift(peaks[:2]) is None


def test_no_shift_gives_zero_slope():
    peaks = [{"excitation_nm": ex, "emission_max_nm": 450.0} for ex in (300, 340, 380)]
    fit = tool.fit_shift(peaks)
    assert fit["slope"] == pytest.approx(0, abs=1e-9) and fit["r2"] is None


def test_scatter_warnings():
    peaks = [{"excitation_nm": 300, "emission_max_nm": 598.0}, {"excitation_nm": 340, "emission_max_nm": 345.0},
             {"excitation_nm": 360, "emission_max_nm": 460.0}]
    notes = tool.scatter_warnings(peaks)
    assert len(notes) == 2 and "second-order" in notes[0] and "Rayleigh" in notes[1]


# ---------- output ----------
@pytest.mark.parametrize("name, magic", [("p.png", b"\x89PNG"), ("p.pdf", b"%PDF"), ("p.svg", b"<?xm")])
def test_plot_formats(tmp_path, name, magic):
    spectra = [spec(ex) for ex in (300, 340, 380)]
    peaks = [tool.find_maximum(s) for s in spectra]
    tool.plot_pl(spectra, peaks, tool.fit_shift(peaks), tmp_path / name, cmap="viridis", stack=0.0,
                 overlay_only=False, normalized=False, title="t", figsize=None, dpi=60)
    assert (tmp_path / name).read_bytes()[:4] == magic


def test_plot_single_spectrum_and_stack(tmp_path):
    s = spec(340)
    tool.plot_pl([s], [tool.find_maximum(s)], None, tmp_path / "one.png", cmap="plasma", stack=0.0,
                 overlay_only=False, normalized=True, title=None, figsize=(4, 3), dpi=60)
    spectra = [spec(300), spec(340)]
    peaks = [tool.find_maximum(x) for x in spectra]
    tool.plot_pl(spectra, peaks, None, tmp_path / "stack.png", cmap="viridis", stack=0.5, overlay_only=True,
                 normalized=False, title=None, figsize=None, dpi=60)
    assert (tmp_path / "one.png").is_file() and (tmp_path / "stack.png").is_file()


def test_peak_table_and_csv(tmp_path):
    peaks = [tool.find_maximum(spec(300, y=band(500, sigma=30))), tool.find_maximum(spec(300, y=band(690, sigma=40)))]
    table = tool.peak_table(peaks)
    assert "n/a" in table and "em max" in table
    tool.write_csv(tmp_path / "m.csv", peaks)
    parsed = list(csv.DictReader((tmp_path / "m.csv").open()))
    assert list(parsed[0]) == tool.PEAK_COLUMNS and parsed[1]["fwhm_nm"] == ""


# ---------- CLI ----------
def test_cli_example_with_masking_recovers_the_shift(tmp_path, monkeypatch, capsys):
    out_csv = tmp_path / "m.csv"
    assert run(monkeypatch, EXAMPLE, "-o", tmp_path / "p.png", "--mask-scatter", 15, "--csv", out_csv, "--dpi", 60) == 0
    captured = capsys.readouterr()
    assert "slope 0.449" in captured.out and captured.err == ""
    parsed = list(csv.DictReader(out_csv.open()))
    assert len(parsed) == 8
    for row in parsed:
        truth = 430 + 0.45 * (float(row["excitation_nm"]) - 300)
        assert float(row["emission_max_nm"]) == pytest.approx(truth, abs=0.7)
    assert (tmp_path / "p.png").read_bytes()[:4] == b"\x89PNG"


def test_cli_without_masking_warns_about_second_order_scatter(tmp_path, monkeypatch, capsys):
    assert run(monkeypatch, EXAMPLE, "-o", tmp_path / "p.png", "--dpi", 60) == 0
    captured = capsys.readouterr()
    assert "[warn] ex 300 nm" in captured.err and "second-order" in captured.err
    assert "599" in captured.out or "600" in captured.out


def test_cli_separate_files_and_options(tmp_path, monkeypatch, capsys):
    files = []
    for ex in (320, 360, 400):
        text = "".join(f"{x},{y}\n" for x, y in zip(WL, band(430 + 0.45 * (ex - 300), noise=0.003, seed=ex)))
        files.append(write(tmp_path, text, f"ex{ex}.csv"))
    assert run(monkeypatch, *files, "-o", tmp_path / "p.png", "--normalize", "--smooth", 11, "--stack", 0.3,
               "--xmin", 340, "--xmax", 680, "--cmap", "plasma", "--overlay-only", "--title", "t", "--dpi", 60) == 0
    assert "[ok] 3 spectra" in capsys.readouterr().out


def test_cli_two_spectra_no_fit_message(tmp_path, monkeypatch, capsys):
    text = "nm,ex320,ex360\n" + "".join(f"{x},{a},{b}\n" for x, a, b in zip(WL, band(439), band(457)))
    assert run(monkeypatch, write(tmp_path, text), "-o", tmp_path / "p.png", "--dpi", 60) == 0
    assert "fewer than 3 maxima" in capsys.readouterr().out


def test_cli_errors(tmp_path, monkeypatch, capsys):
    assert run(monkeypatch, EXAMPLE, "-o", tmp_path / "x.jpg") == 2
    assert run(monkeypatch, EXAMPLE, "-o", tmp_path / "x.png", "--cmap", "nonsense") == 2
    assert run(monkeypatch, tmp_path / "missing.csv", "-o", tmp_path / "x.png") == 2
    assert run(monkeypatch, EXAMPLE, "-o", tmp_path / "x.png", "--smooth", 4) == 2
    assert run(monkeypatch, EXAMPLE, "-o", tmp_path / "x.png", "--ex", 300, 320) == 2
    err = capsys.readouterr().err
    assert "must end in" in err and "unknown colormap" in err and "odd window" in err and "gives 2 wavelengths" in err
