import csv
from pathlib import Path

import numpy as np
import pytest
import tool

EXAMPLE = Path(__file__).parent / "example" / "synthetic_ftir.csv"
X = np.arange(400.0, 4001.0, 2.0)


def bands(*specs, noise=0.0, seed=0, baseline=0.02):
    """Absorbance from (center, FWHM, height) Gaussian bands."""
    y = np.full_like(X, baseline)
    for center, fwhm, height in specs:
        y += height * np.exp(-0.5 * ((X - center) / (fwhm / 2.3548)) ** 2)
    return y + np.random.default_rng(seed).normal(0, noise, X.size) if noise else y


def run(monkeypatch, *argv):
    monkeypatch.setattr("sys.argv", ["tool.py", *map(str, argv)])
    return tool.main()


def write(tmp_path, text, name="s.csv"):
    p = tmp_path / name
    p.write_text(text)
    return p


# ---------- the reference table ----------
UCLA_RANGES = {
    "Alkyl C-H stretch": (2850, 2950), "Alkenyl C-H stretch": (3010, 3100), "Alkenyl C=C stretch": (1620, 1680),
    "Alcohol/phenol O-H stretch (hydrogen-bonded)": (3200, 3550), "Carboxylic acid O-H stretch": (2500, 3000),
    "Amine N-H stretch": (3300, 3500), "Amide N-H stretch": (3500, 3700), "Nitrile C\u2261N stretch": (2220, 2260),
    "Alkynyl C\u2261C stretch": (2100, 2260), "Aldehyde C=O stretch": (1690, 1740), "Ketone C=O stretch": (1680, 1750),
    "Ester C=O stretch": (1735, 1750), "Carboxylic acid C=O stretch": (1710, 1780),
    "Amide C=O stretch (amide I)": (1630, 1690), "Aromatic C=C stretch": (1500, 1700),
    "Aromatic C-H bend (out of plane)": (680, 860),
}


def test_table_is_well_formed():
    groups = [b.group for b in tool.TABLE]
    assert len(groups) == len(set(groups))
    for b in tool.TABLE:
        assert 400 <= b.lo < b.hi <= 4000, b.group
        assert b.intensity in tool.INTENSITY_WORDS and b.shape in ("", "broad", "sharp") and b.source in tool.SOURCES


def test_ranges_taken_from_the_ucla_table_are_unchanged():
    by_group = {b.group: b for b in tool.TABLE}
    for group, (lo, hi) in UCLA_RANGES.items():
        assert (by_group[group].lo, by_group[group].hi) == (lo, hi), group
        assert by_group[group].source == "UCLA"


def test_every_source_is_described_and_unchecked_entries_are_labelled():
    assert {b.source for b in tool.TABLE} <= set(tool.SOURCES)
    assert "not checked" in tool.SOURCES["text"] and "ucla" in tool.SOURCES["UCLA"].lower()
    nh_bend = next(b for b in tool.TABLE if b.group == "Amine N-H bend")
    assert (nh_bend.lo, nh_bend.hi, nh_bend.source) == (1550, 1650, "PMC")


# ---------- signal handling ----------
def test_detect_signal_and_conversion():
    assert tool.detect_signal(np.array([90.0, 40.0])) == "transmittance"
    assert tool.detect_signal(np.array([0.9, 0.4])) == "absorbance"
    assert tool.detect_signal(np.array([90.0]), "absorbance") == "absorbance"
    assert tool.to_absorbance(np.array([10.0, 100.0]), "transmittance").tolist() == pytest.approx([1.0, 0.0])
    assert tool.to_absorbance(np.array([0.1]), "transmittance-fraction")[0] == pytest.approx(1.0)
    assert tool.to_absorbance(np.array([0.3]), "absorbance")[0] == 0.3
    with pytest.raises(ValueError, match="above 0"):
        tool.to_absorbance(np.array([50.0, 0.0]), "transmittance")


# ---------- reading ----------
def test_read_descending_with_header_and_metadata(tmp_path):
    rows = "".join(f"{4000 - 10 * i},{i / 100}\n" for i in range(30))
    spec = tool.read_spectrum(write(tmp_path, "Instrument: demo\nwavenumber,absorbance\n" + rows))
    assert spec.x[0] < spec.x[-1] and spec.name == "absorbance" and spec.x.size == 30


def test_read_column_choice_and_errors(tmp_path):
    rows = "".join(f"{400 + 10 * i},{i / 100},{i / 50}\n" for i in range(30))
    p = write(tmp_path, "cm-1,A,B\n" + rows)
    assert tool.read_spectrum(p, "B").y[-1] == pytest.approx(29 / 50)
    assert tool.read_spectrum(p, "1").name == "A"
    with pytest.raises(ValueError, match="column 'C' not found"):
        tool.read_spectrum(p, "C")
    with pytest.raises(ValueError, match="no numeric data"):
        tool.read_spectrum(write(tmp_path, "text\n"))
    with pytest.raises(ValueError, match="only 3 data points"):
        tool.read_spectrum(write(tmp_path, "1,1\n2,2\n3,3\n"))
    with pytest.raises(ValueError, match="cannot read"):
        tool.read_spectrum(tmp_path / "missing.csv")


# ---------- finding peaks ----------
def test_peak_positions_and_widths():
    a = bands((3400, 380, 0.55), (2925, 40, 0.18), (1710, 45, 0.28), (1100, 120, 0.45), noise=0.002, seed=1)
    peaks = tool.find_peaks_in(X, a)
    assert [round(p["wavenumber_cm1"]) for p in peaks] == pytest.approx([3400, 2925, 1710, 1100], abs=3)
    assert peaks[0]["wavenumber_cm1"] > peaks[-1]["wavenumber_cm1"]  # high to low wavenumber
    assert peaks[0]["fwhm_cm1"] == pytest.approx(380, rel=0.1) and peaks[1]["fwhm_cm1"] == pytest.approx(40, rel=0.15)


def test_intensity_classes_are_relative_to_the_strongest_peak():
    peaks = tool.find_peaks_in(X, bands((3000, 40, 1.0), (2000, 40, 0.4), (1000, 40, 0.1)), prominence=0.05)
    assert [p["intensity_class"] for p in peaks] == ["strong", "medium", "weak"]
    assert tool.classify(0.6) == "strong" and tool.classify(0.25) == "medium" and tool.classify(0.24) == "weak"


def test_prominence_and_separation_limits():
    a = bands((2000, 30, 1.0), (1500, 30, 0.03))
    assert len(tool.find_peaks_in(X, a)) == 1  # the small band is under 5% of the range
    assert len(tool.find_peaks_in(X, a, prominence=0.01)) == 2
    close = bands((2000, 8, 1.0), (2014, 8, 0.8))
    assert len(tool.find_peaks_in(X, close, prominence=0.05, min_separation=10)) == 2
    assert len(tool.find_peaks_in(X, close, prominence=0.05, min_separation=40)) == 1


def test_no_peaks_in_a_flat_spectrum():
    assert tool.find_peaks_in(X, np.full_like(X, 0.1)) == []


def test_refined_position_is_finer_than_the_grid():
    a = bands((1711.3, 40, 0.5), baseline=0.0)
    assert tool.find_peaks_in(X, a)[0]["wavenumber_cm1"] == pytest.approx(1711.3, abs=0.3)


# ---------- suggestions ----------
def groups(position, cls="medium", fwhm=None, **kw):
    return [m["group"] for m in tool.suggest(position, cls, fwhm, **kw)]


def test_typical_assignments():
    assert groups(3400, "strong", 376)[0] == "Alcohol/phenol O-H stretch (hydrogen-bonded)"
    assert groups(2925, "medium", 38)[0] == "Alkyl C-H stretch"
    assert groups(2855, "weak", 27)[0] == "Alkyl C-H stretch"  # weak peaks still go to the range that contains them
    assert groups(2240, "medium", 15)[0] == "Nitrile C\u2261N stretch"
    assert groups(1100, "strong", 120)[0].startswith("C-O stretch")


def test_carbonyl_region_lists_the_carbonyl_groups():
    names = groups(1710, "strong", 34, top=6)
    assert {"Aldehyde C=O stretch", "Ketone C=O stretch", "Carboxylic acid C=O stretch"} <= set(names)
    assert not any("O-H stretch" in n for n in names)


def test_edge_positions_use_a_smooth_penalty():
    just_outside = tool.suggest(1709.6, "strong", 34, top=10)
    assert "Carboxylic acid C=O stretch" in [m["group"] for m in just_outside]
    scores = {m["group"]: m["score"] for m in tool.suggest(1710.0, "strong", 34, top=10)}
    shifted = {m["group"]: m["score"] for m in tool.suggest(1700.0, "strong", 34, top=10)}
    assert scores["Carboxylic acid C=O stretch"] > shifted["Carboxylic acid C=O stretch"]


def test_broad_and_sharp_shapes_change_the_ranking():
    broad = groups(3000, "medium", 300, tolerance=60)
    narrow = groups(3000, "medium", 20, tolerance=60)
    assert broad.index("Carboxylic acid O-H stretch") < narrow.index("Carboxylic acid O-H stretch")


def test_outside_the_table_gives_nothing():
    assert tool.suggest(2350, "weak", 10) == []  # where atmospheric CO2 absorbs
    assert tool.suggest(2262, "medium", 10, tolerance=0) == []
    assert tool.suggest(2270, "medium", 10, tolerance=10) != []


def test_top_limits_the_list_and_scores_are_sorted():
    result = tool.suggest(1710, "strong", 34, top=2)
    assert len(result) == 2 and result[0]["score"] >= result[1]["score"]


def test_describe_marks_fingerprint_and_empty_matches():
    peaks = [{"wavenumber_cm1": 1100.0, "height": 0.5, "prominence": 0.5, "relative": 1.0, "intensity_class": "strong",
              "fwhm_cm1": 120.0},
             {"wavenumber_cm1": 2350.0, "height": 0.1, "prominence": 0.1, "relative": 0.2, "intensity_class": "weak",
              "fwhm_cm1": 10.0}]
    text = "\n".join(tool.describe(peaks, 3, 10.0))
    assert "fingerprint region" in text and "no band in the table" in text


# ---------- output ----------
@pytest.mark.parametrize("signal, magic, name", [("transmittance", b"\x89PNG", "p.png"), ("absorbance", b"%PDF", "p.pdf"),
                                                 ("transmittance", b"<?xm", "p.svg")])
def test_plot_formats(tmp_path, signal, magic, name):
    a = bands((3400, 380, 0.55), (1710, 45, 0.28))
    y = 100 * 10 ** (-a) if signal == "transmittance" else a
    spec = tool.Spectrum("s", X, y)
    peaks = tool.find_peaks_in(X, a)
    tool.plot_ftir(spec, signal, peaks, tmp_path / name, regions=True, max_labels=5, title="t", figsize=(6, 3), dpi=60)
    assert (tmp_path / name).read_bytes()[:4] == magic


def test_csv_columns_and_content(tmp_path):
    peaks = tool.find_peaks_in(X, bands((2925, 40, 0.3), (1710, 45, 0.28)))
    tool.write_csv(tmp_path / "p.csv", peaks, 2, 10.0)
    rows = list(csv.DictReader((tmp_path / "p.csv").open()))
    assert list(rows[0]) == tool.PEAK_COLUMNS and len(rows) == 2
    assert "Alkyl C-H stretch" in rows[0]["suggestions"] and " | " in rows[1]["suggestions"]


# ---------- CLI ----------
def test_cli_example_spectrum(tmp_path, monkeypatch, capsys):
    out_csv = tmp_path / "p.csv"
    assert run(monkeypatch, EXAMPLE, "-o", tmp_path / "p.png", "--csv", out_csv, "--dpi", 60) == 0
    out = capsys.readouterr().out
    assert "8 peaks" in out and "percent transmittance" in out
    rows = list(csv.DictReader(out_csv.open()))
    found = [float(r["wavenumber_cm1"]) for r in rows]
    assert found == pytest.approx([3400, 2925, 2855, 1710, 1620, 1400, 1260, 1100], abs=3)
    assert "Alkyl C-H stretch" in rows[1]["suggestions"] and "Carboxylic acid C=O stretch" in rows[3]["suggestions"]
    assert (tmp_path / "p.png").read_bytes()[:4] == b"\x89PNG"


def test_cli_range_smoothing_and_options(tmp_path, monkeypatch, capsys):
    assert run(monkeypatch, EXAMPLE, "--xmin", 1000, "--xmax", 1800, "--smooth", 11, "--top", 1, "--prominence", 0.05) == 0
    out = capsys.readouterr().out
    assert "peaks in 1000-1800" in out and "2925" not in out
    assert run(monkeypatch, EXAMPLE, "-o", tmp_path / "p.svg", "--no-regions", "--max-labels", 3, "--title", "t") == 0


def test_cli_absorbance_input_and_warning(tmp_path, monkeypatch, capsys):
    a = bands((2925, 40, 0.3), (1710, 45, 0.28))
    p = write(tmp_path, "".join(f"{x},{v}\n" for x, v in zip(X, a)))
    assert run(monkeypatch, p) == 0
    out = capsys.readouterr()
    assert "signal read as absorbance" in out.out and "2 peaks" in out.out
    assert run(monkeypatch, EXAMPLE, "--signal", "absorbance") == 0
    assert "look like percent transmittance" in capsys.readouterr().err


def test_cli_errors(tmp_path, monkeypatch, capsys):
    assert run(monkeypatch, EXAMPLE, "-o", tmp_path / "x.jpg") == 2
    assert run(monkeypatch, tmp_path / "missing.csv") == 2
    assert run(monkeypatch, EXAMPLE, "--smooth", 4) == 2
    assert run(monkeypatch, EXAMPLE, "--xmin", 1000, "--xmax", 1004) == 2
    assert run(monkeypatch, EXAMPLE, "--top", 0) == 2
    assert run(monkeypatch, EXAMPLE, "--column", "nope") == 2
    err = capsys.readouterr().err
    assert "must end in" in err and "odd window" in err and "fewer than 10 points" in err and "--top must be" in err
