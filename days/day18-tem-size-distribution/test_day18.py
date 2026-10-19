import csv
from pathlib import Path

import numpy as np
import pytest
import tool

EXAMPLE = Path(__file__).parent / "example"
SIZES, IMAGEJ = EXAMPLE / "synthetic_sizes.csv", EXAMPLE / "synthetic_imagej.csv"


def lognormal_sample(n, median=4.2, gsd=1.32, seed=0):
    return np.random.default_rng(seed).lognormal(np.log(median), np.log(gsd), n)


def write(tmp_path, text, name="s.csv"):
    p = tmp_path / name
    p.write_text(text)
    return p


def run(monkeypatch, *argv):
    monkeypatch.setattr("sys.argv", ["tool.py", *map(str, argv)])
    return tool.main()


# ---------- reading ----------
def test_single_column_with_and_without_header(tmp_path):
    names, data = tool.read_columns(write(tmp_path, "diameter\n3.1\n4.2\n5.3\n6.4\n7.5\n"))
    assert names == ["diameter"] and data.shape == (5, 1)
    names, data = tool.read_columns(write(tmp_path, "3.1\n4.2\n5.3\n6.4\n7.5\n"))
    assert names == ["column 1"] and data[:, 0].tolist() == [3.1, 4.2, 5.3, 6.4, 7.5]


def test_pick_column_by_name_number_and_errors(tmp_path):
    path = write(tmp_path, "id,Feret,Area\n1,3.0,7\n2,4.0,12\n3,5.0,20\n4,6.0,28\n5,7.0,38\n")
    names, data = tool.read_columns(path)
    assert tool.pick_column(names, data, "feret", path)[1].tolist() == [3, 4, 5, 6, 7]  # not case sensitive
    assert tool.pick_column(names, data, "3", path)[0] == "Area"
    with pytest.raises(ValueError, match="choose one with --column"):
        tool.pick_column(names, data, None, path)
    with pytest.raises(ValueError, match="column 'nope' not found. Available: id, Feret, Area"):
        tool.pick_column(names, data, "nope", path)


def test_imagej_style_blank_first_header_and_decimal_comma(tmp_path):
    names, data = tool.read_columns(write(tmp_path, " ,Area,Mean\n1,100,5\n2,200,6\n"))
    assert names == ["column 1", "Area", "Mean"] and data.shape == (2, 3)
    names, data = tool.read_columns(write(tmp_path, "d;n\n3,5;1\n4,5;2\n"))
    assert data[:, 0].tolist() == [3.5, 4.5]


def test_read_errors(tmp_path):
    with pytest.raises(ValueError, match="no numbers found"):
        tool.read_columns(write(tmp_path, "text only\n"))
    with pytest.raises(ValueError, match="cannot read"):
        tool.read_columns(tmp_path / "missing.csv")


# ---------- converting to diameters ----------
def test_diameters_from_area_and_scale():
    radii = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    areas = np.pi * radii**2
    assert tool.to_diameters(areas, from_area=True).tolist() == pytest.approx((2 * radii).tolist())
    assert tool.to_diameters(areas, from_area=True, scale=0.25).tolist() == pytest.approx((0.5 * radii).tolist())
    assert tool.to_diameters(radii, scale=0.5).tolist() == pytest.approx((0.5 * radii).tolist())


def test_diameter_validation():
    with pytest.raises(ValueError, match="at least 5"):
        tool.to_diameters(np.array([1.0, 2.0, np.nan, 3.0]))
    with pytest.raises(ValueError, match="2 values are 0 or negative"):
        tool.to_diameters(np.array([1.0, 0.0, -2.0, 3.0, 4.0, 5.0]))
    with pytest.raises(ValueError, match="--unit-scale must be positive"):
        tool.to_diameters(np.ones(6), scale=0)
    assert tool.to_diameters(np.array([1.0, 2.0, np.nan, 3.0, 4.0, 5.0])).size == 5  # blanks are dropped


# ---------- the statistics ----------
def test_exact_values_on_a_small_sample():
    d = np.array([1.0, 2.0, 4.0, 8.0, 16.0])
    s = tool.describe(d)
    assert s["mean"] == pytest.approx(6.2) and s["sd"] == pytest.approx(np.std(d, ddof=1)) and s["median"] == 4.0
    assert s["geometric_mean"] == pytest.approx(4.0)  # ln d is 0, ln2, 2 ln2, 3 ln2, 4 ln2
    assert s["sigma_ln"] == pytest.approx(np.log(2) * np.std([0, 1, 2, 3, 4], ddof=1))
    assert s["geometric_sd"] == pytest.approx(np.exp(s["sigma_ln"]))
    assert s["lognormal_mean"] == pytest.approx(np.exp(s["mu_ln"] + s["sigma_ln"] ** 2 / 2))
    assert s["lognormal_mode"] == pytest.approx(np.exp(s["mu_ln"] - s["sigma_ln"] ** 2))
    assert s["pdi"] == pytest.approx((s["sd"] / s["mean"]) ** 2) and s["cv_percent"] == pytest.approx(100 * s["sd"] / s["mean"])
    t = 2.7764451  # Student t, 4 degrees of freedom, 97.5%
    assert s["mean_ci95_high"] - s["mean"] == pytest.approx(t * s["sd"] / np.sqrt(5), rel=1e-6)
    assert s["sigma_ln_se"] == pytest.approx(s["sigma_ln"] / np.sqrt(8)) and s["mu_ln_se"] == pytest.approx(s["sigma_ln"] / np.sqrt(5))


def test_identical_values_are_rejected():
    with pytest.raises(ValueError, match="all the values are identical"):
        tool.describe(np.full(10, 5.0))


def test_parameters_are_recovered_from_a_large_sample():
    s = tool.describe(lognormal_sample(50000, 4.2, 1.32, seed=1))
    assert s["geometric_mean"] == pytest.approx(4.2, rel=0.005) and s["geometric_sd"] == pytest.approx(1.32, rel=0.005)
    assert s["shapiro_p_normal"] is None  # Shapiro-Wilk is only used up to 5000 values


def test_confidence_intervals_have_the_stated_coverage():
    """About 95% of the 95% intervals must contain the true value (simulated samples of 100 particles)."""
    mu, sigma = np.log(4.2), np.log(1.32)
    rng = np.random.default_rng(5)
    inside_gm = inside_mean = inside_sigma = 0
    runs = 600
    for _ in range(runs):
        s = tool.describe(rng.lognormal(mu, sigma, 100))
        inside_gm += s["geometric_mean_ci95_low"] <= 4.2 <= s["geometric_mean_ci95_high"]
        inside_mean += s["mean_ci95_low"] <= np.exp(mu + sigma**2 / 2) <= s["mean_ci95_high"]
        inside_sigma += abs(s["sigma_ln"] - sigma) <= 1.96 * s["sigma_ln_se"]
    for hits in (inside_gm, inside_mean, inside_sigma):
        assert 0.92 <= hits / runs <= 0.98


def test_aic_prefers_the_right_model():
    assert tool.describe(lognormal_sample(300, 4.0, 1.6, seed=2))["delta_aic"] > 10
    normal = np.random.default_rng(3).normal(50, 3, 300)
    assert tool.describe(normal)["delta_aic"] < 3  # no strong preference for the lognormal


def test_shapiro_on_log_sizes_flags_two_populations():
    rng = np.random.default_rng(4)
    two = np.concatenate([rng.lognormal(np.log(3), 0.1, 150), rng.lognormal(np.log(9), 0.1, 150)])
    assert tool.describe(two)["shapiro_p_lognormal"] < 0.001
    assert tool.describe(lognormal_sample(300, seed=6))["shapiro_p_lognormal"] > 0.01


# ---------- histogram bins ----------
def test_bin_edges_options():
    d = lognormal_sample(200)
    assert len(tool.bin_edges(d, "12")) == 13
    for rule in ("auto", "fd", "sturges", "scott", "doane", "rice", "sqrt"):
        assert len(tool.bin_edges(d, rule)) >= 2
    edges = tool.bin_edges(d, "auto", width=0.5)
    assert np.allclose(np.diff(edges), 0.5) and edges[0] <= d.min() and edges[-1] >= d.max()
    with pytest.raises(ValueError, match="unknown --bins"):
        tool.bin_edges(d, "magic")
    with pytest.raises(ValueError, match="at least 1"):
        tool.bin_edges(d, "0")
    with pytest.raises(ValueError, match="--bin-width must be positive"):
        tool.bin_edges(d, "auto", width=-1)


# ---------- warnings ----------
def notes(d, resolution=None, bins="auto"):
    return tool.warnings_for(d, tool.describe(d), tool.bin_edges(d, bins), resolution)


def test_clean_sample_gives_no_warnings():
    assert notes(lognormal_sample(300, seed=8)) == []


def test_small_samples_are_flagged_with_different_strength():
    assert any("very uncertain" in n for n in notes(lognormal_sample(20)))
    small = notes(lognormal_sample(60, seed=9))
    assert any("small sample" in n for n in small) and not any("very uncertain" in n for n in small)


def test_resolution_tail_and_model_warnings():
    d = lognormal_sample(300, seed=10)
    assert any("narrower than the measurement resolution" in n for n in notes(d, resolution=0.5, bins="60"))
    with_aggregate = np.append(d, 40.0)
    assert any("more than 5 times the median" in n for n in notes(with_aggregate))
    two = np.concatenate([np.random.default_rng(11).lognormal(np.log(3), 0.1, 150),
                          np.random.default_rng(12).lognormal(np.log(9), 0.1, 150)])
    assert any("not normal at the 5% level" in n for n in notes(two))
    normal = np.random.default_rng(13).normal(50, 2, 400)
    s = tool.describe(normal)
    s["delta_aic"] = -4.0
    assert any("normal distribution fits better" in n for n in tool.warnings_for(normal, s, tool.bin_edges(normal), None))


# ---------- output ----------
@pytest.mark.parametrize("name, magic", [("p.png", b"\x89PNG"), ("p.pdf", b"%PDF"), ("p.svg", b"<?xm")])
def test_plot_formats_and_options(tmp_path, name, magic):
    d = lognormal_sample(200)
    s, edges = tool.describe(d), tool.bin_edges(d)
    tool.plot_distribution(d, s, edges, tmp_path / name, density=False, also_normal=True, qq=True, unit="nm", title="t",
                           figsize=None, dpi=60)
    assert (tmp_path / name).read_bytes()[:4] == magic
    tool.plot_distribution(d, s, edges, tmp_path / ("d" + name), density=True, also_normal=False, qq=False, unit="um",
                           title=None, figsize=(4, 3), dpi=60)
    assert (tmp_path / ("d" + name)).is_file()


def test_summary_and_histogram_csv(tmp_path):
    d = lognormal_sample(200)
    s, edges = tool.describe(d), tool.bin_edges(d, "10")
    tool.write_summary_csv(tmp_path / "s.csv", s)
    row = next(iter(csv.DictReader((tmp_path / "s.csv").open())))
    assert list(row) == tool.SUMMARY_COLUMNS and row["n"] == "200"
    assert float(row["geometric_mean"]) == pytest.approx(s["geometric_mean"], rel=1e-5)
    tool.write_hist_csv(tmp_path / "h.csv", d, edges)
    bins = list(csv.DictReader((tmp_path / "h.csv").open()))
    assert len(bins) == 10 and sum(int(b["count"]) for b in bins) == 200
    assert sum(float(b["fraction"]) for b in bins) == pytest.approx(1.0, abs=1e-4)


def test_summary_csv_leaves_unavailable_values_blank(tmp_path):
    s = tool.describe(lognormal_sample(50000, seed=14))
    tool.write_summary_csv(tmp_path / "s.csv", s)
    assert next(iter(csv.DictReader((tmp_path / "s.csv").open())))["shapiro_p_normal"] == ""


# ---------- CLI ----------
def test_cli_synthetic_sizes(tmp_path, monkeypatch, capsys):
    out_csv = tmp_path / "s.csv"
    assert run(monkeypatch, SIZES, "--column", "diameter_nm", "-o", tmp_path / "p.png", "--csv", out_csv,
               "--hist-csv", tmp_path / "h.csv", "--dpi", 60) == 0
    out = capsys.readouterr().out
    assert "n = 250 particles" in out and "lognormal fit" in out
    row = next(iter(csv.DictReader(out_csv.open())))
    assert float(row["geometric_mean_ci95_low"]) <= 4.2 <= float(row["geometric_mean_ci95_high"])
    assert float(row["geometric_sd"]) == pytest.approx(1.32, abs=0.05)
    assert (tmp_path / "p.png").read_bytes()[:4] == b"\x89PNG" and (tmp_path / "h.csv").is_file()


def test_cli_imagej_areas_with_scale(monkeypatch, capsys):
    assert run(monkeypatch, IMAGEJ, "--column", "Area", "--from-area", "--unit-scale", 0.25) == 0
    out = capsys.readouterr().out
    assert "equivalent circle diameter, scaled by 0.25" in out
    geometric_mean = float(out.split("geometric mean")[1].split()[0])
    assert geometric_mean == pytest.approx(6.0, abs=0.35)


def test_cli_options_and_warnings(tmp_path, monkeypatch, capsys):
    assert run(monkeypatch, SIZES, "--column", "2", "--bin-width", 0.1, "--unit-scale", 1.0, "--density", "--also-normal",
               "--no-qq", "--unit-label", "nm", "--title", "t", "--figsize", 5, 3, "-o", tmp_path / "p.svg") == 0
    few = write(tmp_path, "".join(f"{x}\n" for x in (3.1, 3.5, 4.0, 4.4, 5.1, 6.0, 7.2)), "few.csv")
    assert run(monkeypatch, few) == 0
    assert "very uncertain" in capsys.readouterr().err


def test_cli_errors(tmp_path, monkeypatch, capsys):
    assert run(monkeypatch, IMAGEJ) == 2  # several columns, none chosen
    assert run(monkeypatch, SIZES, "--column", "nope") == 2
    assert run(monkeypatch, tmp_path / "missing.csv") == 2
    assert run(monkeypatch, SIZES, "--column", "diameter_nm", "-o", tmp_path / "x.jpg") == 2
    assert run(monkeypatch, SIZES, "--column", "diameter_nm", "--bins", "magic") == 2
    assert run(monkeypatch, SIZES, "--column", "diameter_nm", "--unit-scale", 0) == 2
    same = write(tmp_path, "5\n5\n5\n5\n5\n5\n", "same.csv")
    assert run(monkeypatch, same) == 2
    err = capsys.readouterr().err
    for text in ("choose one with --column", "column 'nope' not found", "must end in", "unknown --bins",
                 "--unit-scale must be positive", "all the values are identical"):
        assert text in err, text
