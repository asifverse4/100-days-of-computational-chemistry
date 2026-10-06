from pathlib import Path

import pytest
import tool

EXAMPLE = Path(__file__).parent / "example"
WATER = "3\nwater\nO 0.0 0.0 0.117\nH 0.0 0.757 -0.469\nH 0.0 -0.757 -0.469\n"


def frames_equal(a, b, tol=1e-5):
    assert [f.symbols for f in a] == [f.symbols for f in b]
    for fa, fb in zip(a, b):
        for ca, cb in zip(fa.coords, fb.coords):
            assert ca == pytest.approx(cb, abs=tol)


def test_read_xyz_multi_frame_and_blank_lines():
    frames = tool.read_xyz(WATER + "\n" + WATER.replace("water", "water2"))
    assert len(frames) == 2
    assert frames[1].comment == "water2"
    assert frames[0].symbols == ["O", "H", "H"]


def test_read_xyz_normalizes_symbols_and_ignores_extra_columns():
    f = tool.read_xyz("2\n\nCL 0 0 0 0.1\nna 1 0 0 0.2\n")[0]
    assert f.symbols == ["Cl", "Na"]


@pytest.mark.parametrize("text, fragment", [
    ("hello\nworld\n", "atom count"),
    ("3\nc\nO 0 0 0\n", "does not have them"),
    ("1\nc\nO 0 0\n", "symbol x y z"),
    ("1\nc\nO a b c\n", "line 3"),
    ("1\nc\nO1 0 0 0\n", "element symbol"),
    ("", "no structures"),
])
def test_read_xyz_errors(text, fragment):
    with pytest.raises(ValueError, match=fragment):
        tool.read_xyz(text)


def test_formula_hill_order():
    assert tool.formula(tool.read_xyz("4\n\nC 0 0 0\nH 0 0 1\nH 0 1 0\nO 1 0 0\n")[0]) == "CH2O"
    assert tool.formula(tool.read_xyz(WATER)[0]) == "H2O"
    assert tool.formula(tool.read_xyz("3\n\nCl 0 0 0\nNa 1 0 0\nC 2 0 0\n")[0]) == "CClNa"


def test_centroid_and_center_of_mass():
    f = tool.read_xyz(WATER)[0]
    c = tool.centered(f, "centroid")
    assert tool.centroid(c) == pytest.approx([0, 0, 0], abs=1e-12)
    m = tool.centered(f, "mass")
    assert tool.center_of_mass(m) == pytest.approx([0, 0, 0], abs=1e-12)
    # oxygen is heavy, so the center of mass sits closer to O (z = 0.117) than the centroid does
    assert tool.center_of_mass(f)[2] > tool.centroid(f)[2]


def test_center_of_mass_unknown_element():
    f = tool.read_xyz("1\n\nUu 0 0 0\n")[0]
    with pytest.raises(ValueError, match="no atomic mass"):
        tool.center_of_mass(f)


def test_pdb_columns_and_round_trip():
    frames = tool.read_xyz(WATER)
    pdb = tool.write_pdb(frames)
    atom_lines = [ln for ln in pdb.splitlines() if ln.startswith("HETATM")]
    assert len(atom_lines) == 3
    assert all(len(ln) == 78 for ln in atom_lines)
    assert float(atom_lines[1][38:46]) == pytest.approx(0.757)
    assert atom_lines[0][76:78] == " O"
    frames_equal(frames, tool.read_pdb(pdb), tol=1e-3)


def test_pdb_multi_model_and_element_from_name():
    frames = tool.read_xyz(WATER) * 2
    back = tool.read_pdb(tool.write_pdb(frames))
    assert len(back) == 2
    line = "ATOM      1 CL1  MOL A   1       1.000   2.000   3.000  1.00  0.00"  # 2-letter names start in column 13
    assert tool.read_pdb(line + "\n")[0].symbols == ["Cl"]
    line = "ATOM      1  C1  MOL A   1       1.000   2.000   3.000  1.00  0.00"
    assert tool.read_pdb(line + "\n")[0].symbols == ["C"]


def test_coord_round_trip_uses_bohr():
    frames = tool.read_xyz(WATER)
    text = tool.write_coord(frames)
    assert text.startswith("$coord") and text.rstrip().endswith("$end")
    assert text.splitlines()[1].split()[-1] == "o"
    frames_equal(frames, tool.read_coord(text), tol=1e-9)
    # 1 angstrom = 1.8897 bohr
    one = tool.write_coord([tool.Frame(["H"], [[1.0, 0.0, 0.0]])])
    assert float(one.splitlines()[1].split()[0]) == pytest.approx(1.8897261, abs=1e-6)


def test_coord_rejects_multiple_frames_and_bad_files():
    with pytest.raises(ValueError, match="one structure"):
        tool.write_coord(tool.read_xyz(WATER) * 2)
    with pytest.raises(ValueError, match=r"\$coord"):
        tool.read_coord("nothing here")


def test_file_format_detection(tmp_path):
    assert tool.file_format(Path("a.XYZ")) == "xyz"
    assert tool.file_format(Path("coord")) == "coord"
    with pytest.raises(ValueError, match="unsupported"):
        tool.file_format(Path("a.cif"))


def run(monkeypatch, *argv):
    monkeypatch.setattr("sys.argv", ["tool.py", *map(str, argv)])
    return tool.main()


def test_cli_merge_split_round_trip(tmp_path, monkeypatch):
    merged = tmp_path / "all.xyz"
    assert run(monkeypatch, "merge", EXAMPLE / "ethanol.xyz", EXAMPLE / "benzene.xyz", "-o", merged) == 0
    assert [len(f.symbols) for f in tool.read_xyz(merged.read_text())] == [9, 12]
    assert run(monkeypatch, "split", merged, "--outdir", tmp_path / "fr") == 0
    parts = sorted((tmp_path / "fr").glob("*.xyz"))
    assert [p.name for p in parts] == ["all_001.xyz", "all_002.xyz"]
    frames_equal(tool.read_xyz(parts[1].read_text()), tool.read_xyz((EXAMPLE / "benzene.xyz").read_text()))


def test_cli_merge_strict_rejects_different_molecules(tmp_path, monkeypatch, capsys):
    rc = run(monkeypatch, "merge", EXAMPLE / "ethanol.xyz", EXAMPLE / "benzene.xyz",
             "-o", tmp_path / "x.xyz", "--strict")
    assert rc == 2
    assert "--strict" in capsys.readouterr().err


def test_cli_center_default_output_name_and_modes(tmp_path, monkeypatch):
    src = tmp_path / "w.xyz"
    src.write_text(WATER)
    assert run(monkeypatch, "center", src, "--mode", "mass") == 0
    out = tool.read_xyz((tmp_path / "w_centered.xyz").read_text())[0]
    assert tool.center_of_mass(out) == pytest.approx([0, 0, 0], abs=1e-5)


def test_cli_convert_and_info(tmp_path, monkeypatch, capsys):
    pdb = tmp_path / "e.pdb"
    assert run(monkeypatch, "convert", EXAMPLE / "ethanol.xyz", "-o", pdb) == 0
    assert run(monkeypatch, "info", pdb) == 0
    out = capsys.readouterr().out
    assert "atoms: 9" in out and "formula: C2H6O" in out


def test_cli_errors_return_2(tmp_path, monkeypatch, capsys):
    assert run(monkeypatch, "info", tmp_path / "missing.xyz") == 2
    assert "cannot read" in capsys.readouterr().err
    bad = tmp_path / "bad.xyz"
    bad.write_text("nope")
    assert run(monkeypatch, "info", bad) == 2
