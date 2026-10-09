from pathlib import Path

import pytest
import tool
from PIL import Image
from rdkit import Chem

EXAMPLE = Path(__file__).parent / "example" / "input.smi"
PNG_MAGIC = b"\x89PNG\r\n\x1a\n"


def run(monkeypatch, *argv):
    monkeypatch.setattr("sys.argv", ["tool.py", *map(str, argv)])
    return tool.main()


def test_read_entries(tmp_path):
    f = tmp_path / "in.smi"
    f.write_text("# comment\nCCO ethanol\nc1ccccc1\n\n")
    assert tool.read_entries(str(f)) == [("CCO", "ethanol"), ("c1ccccc1", "mol3")]
    assert tool.read_entries("CCO") == [("CCO", "molecule")]
    assert tool.read_entries("C" * 400) == [("C" * 400, "molecule")]


@pytest.mark.parametrize("kind, expected", [
    ("name", "ethanol"), ("smiles", "CCO"), ("formula", "C2H6O"), ("mw", "ethanol  46.1"), ("none", ""),
])
def test_make_legend(kind, expected):
    assert tool.make_legend(Chem.MolFromSmiles("CCO"), "CCO", "ethanol", kind) == expected


def test_grid_size():
    assert tool.grid_size(3, 300, None, 300) == 300
    assert tool.grid_size(3, 300, 6.5, 300) == 650
    assert tool.grid_size(4, 300, 3.3, 600) == 495
    assert tool.grid_size(3, 300, 0.01, 72) == 50  # never absurdly small
    with pytest.raises(ValueError):
        tool.grid_size(3, 300, -1, 300)


def test_highlight_lists_and_hits():
    mols = [Chem.MolFromSmiles(s) for s in ("c1ccccc1C", "CCO")]
    patt = tool.pattern_from_smarts("c1ccccc1")
    atoms, bonds, hits = tool.highlight_lists(mols, patt)
    assert hits == 1
    assert len(atoms[0]) == 6 and len(bonds[0]) == 6  # six ring atoms, six ring bonds (not the methyl bond)
    assert atoms[1] == [] and bonds[1] == []


def test_bad_smarts():
    with pytest.raises(ValueError, match="invalid SMARTS"):
        tool.pattern_from_smarts("c1cc(")


def test_align_counts_only_matching_molecules():
    mols = [Chem.MolFromSmiles(s) for s in ("c1ccccc1C", "Cc1ccccc1O", "CCO")]
    assert tool.align_to_pattern(mols, tool.pattern_from_smarts("c1ccccc1")) == 2


def test_build_image_dimensions():
    mols = [Chem.MolFromSmiles("CCO")] * 4
    img = tool.build_image(mols, ["a", "b", "c", "d"], 2, 200)
    assert img.size == (400, 400)


def test_page_paths():
    assert tool.page_paths(Path("g.png"), 1) == [Path("g.png")]
    assert [p.name for p in tool.page_paths(Path("out/g.png"), 3)] == ["g_1.png", "g_2.png", "g_3.png"]


def test_cli_png(tmp_path, monkeypatch):
    out = tmp_path / "grid.png"
    assert run(monkeypatch, EXAMPLE, "-o", out, "--cols", 3, "--size", 200) == 0
    assert out.read_bytes()[:8] == PNG_MAGIC
    img = Image.open(out)
    assert img.size == (600, 400)  # 3 columns x 2 rows of 200 px
    assert img.info["dpi"][0] == pytest.approx(300, abs=1)


def test_cli_figure_width_and_dpi(tmp_path, monkeypatch):
    out = tmp_path / "fig.png"
    assert run(monkeypatch, EXAMPLE, "-o", out, "--cols", 3, "--width-in", 6.5, "--dpi", 300) == 0
    img = Image.open(out)
    assert img.size[0] == 1950  # 6.5 in x 300 dpi
    assert img.info["dpi"][0] == pytest.approx(300, abs=1)


def test_cli_svg(tmp_path, monkeypatch):
    out = tmp_path / "grid.svg"
    assert run(monkeypatch, EXAMPLE, "-o", out, "--legend", "formula") == 0
    text = out.read_text()
    assert "<svg" in text and "</svg>" in text


def test_cli_pagination(tmp_path, monkeypatch):
    out = tmp_path / "pages.png"
    assert run(monkeypatch, EXAMPLE, "-o", out, "--per-page", 4, "--cols", 2) == 0
    assert (tmp_path / "pages_1.png").is_file() and (tmp_path / "pages_2.png").is_file()
    assert not out.exists()


def test_cli_highlight_align_and_options(tmp_path, monkeypatch, capsys):
    out = tmp_path / "h.png"
    assert run(monkeypatch, EXAMPLE, "-o", out, "--highlight", "c1ccccc1", "--align",
               "--bw", "--atom-indices", "--legend", "none") == 0
    text = capsys.readouterr().out
    assert "found in 4 of 6" in text and "oriented 4" in text
    assert out.is_file()


def test_cli_invalid_smiles_skipped_exit_1(tmp_path, monkeypatch, capsys):
    src = tmp_path / "in.smi"
    src.write_text("CCO ok\nnot_a_smiles bad\n")
    out = tmp_path / "g.png"
    assert run(monkeypatch, src, "-o", out) == 1
    assert out.is_file()
    assert "[skip] bad" in capsys.readouterr().err


def test_cli_errors(tmp_path, monkeypatch, capsys):
    bad = tmp_path / "bad.smi"
    bad.write_text("nope1 a\nnope2 b\n")
    assert run(monkeypatch, bad, "-o", tmp_path / "x.png") == 2
    assert run(monkeypatch, EXAMPLE, "-o", tmp_path / "x.jpg") == 2
    assert run(monkeypatch, EXAMPLE, "-o", tmp_path / "x.png", "--highlight", "c1cc(") == 2
    assert run(monkeypatch, EXAMPLE, "-o", tmp_path / "x.png", "--align") == 2
    assert run(monkeypatch, EXAMPLE, "-o", tmp_path / "x.png", "--cols", 0) == 2
    assert run(monkeypatch, EXAMPLE, "-o", tmp_path / "x.png", "--width-in", -2) == 2
    assert "needs --highlight" in capsys.readouterr().err
