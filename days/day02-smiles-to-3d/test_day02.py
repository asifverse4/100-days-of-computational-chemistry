import pytest
import tool


def test_ethanol_geometry_and_forcefield():
    mol, ff, energy = tool.smiles_to_mol("CCO")
    assert mol.GetNumAtoms() == 9  # 3 heavy atoms + 6 hydrogens
    assert ff == "MMFF94"
    assert mol.GetNumConformers() == 1
    assert energy == pytest.approx(tool.smiles_to_mol("CCO")[2])  # same seed, same result


def test_aspirin_atom_count():
    assert tool.smiles_to_mol("CC(=O)Oc1ccccc1C(=O)O")[0].GetNumAtoms() == 21


def test_invalid_smiles_raises():
    with pytest.raises(ValueError):
        tool.smiles_to_mol("not_a_smiles")


def test_write_xyz_format(tmp_path):
    mol, _, _ = tool.smiles_to_mol("CCO")
    out = tmp_path / "ethanol.xyz"
    tool.write_xyz(mol, out, "ethanol test")
    lines = out.read_text().splitlines()
    assert lines[0] == "9" and lines[1] == "ethanol test"
    assert len(lines) == 11
    assert lines[2].split()[0] == "C" and len(lines[2].split()) == 4


def test_batch_cli_writes_files_and_reports_failures(tmp_path, monkeypatch, capsys):
    smi = tmp_path / "in.smi"
    smi.write_text("# header\nCCO ethanol\nbroken bad\n")
    monkeypatch.setattr("sys.argv", ["tool.py", str(smi), "--outdir", str(tmp_path / "out")])
    assert tool.main() == 1
    assert (tmp_path / "out" / "ethanol.xyz").is_file()
    assert "[skip] bad" in capsys.readouterr().err
