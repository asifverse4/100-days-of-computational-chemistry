import pytest
from tool import describe, read_entries


def test_aspirin():
    r = describe("CC(=O)Oc1ccccc1C(=O)O", "aspirin")
    assert r["formula"] == "C9H8O4"
    assert r["mw"] == pytest.approx(180.16, abs=0.01)
    assert r["hbd"] == 1
    assert r["hba"] == 4  # Lipinski definition: N + O atoms
    assert r["lipinski_violations"] == 0


def test_ethanol_logp_not_negative_zero():
    r = describe("CCO", "ethanol")
    assert str(r["logp"]) == "0.0"


def test_caffeine_hba_lipinski():
    assert describe("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "caffeine")["hba"] == 6


def test_violation_counted():
    # tetradecane-like long alkane: logP > 5
    assert describe("C" * 20, "alkane")["lipinski_violations"] >= 1


def test_invalid_smiles():
    with pytest.raises(ValueError):
        describe("not_a_smiles", "bad")


def test_long_smiles_not_treated_as_path():
    assert read_entries("C" * 400) == [("C" * 400, "molecule")]
