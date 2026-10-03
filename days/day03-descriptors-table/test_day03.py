import pytest
from tool import describe


def test_aspirin():
    r = describe("CC(=O)Oc1ccccc1C(=O)O", "aspirin")
    assert r["formula"] == "C9H8O4"
    assert r["mw"] == pytest.approx(180.16, abs=0.01)
    assert r["hbd"] == 1
    assert r["lipinski_violations"] == 0


def test_invalid_smiles():
    with pytest.raises(ValueError):
        describe("not_a_smiles", "bad")
