import csv

import pytest
import tool


def run(monkeypatch, *argv):
    monkeypatch.setattr("sys.argv", ["tool.py", *map(str, argv)])
    return tool.main()


# ---------- data ----------
def test_dataset_is_complete_and_consistent():
    els = tool.load_elements()
    assert [e["atomic_number"] for e in els] == list(range(1, 119))
    assert len({e["symbol"] for e in els}) == 118
    assert all(e["atomic_mass"] > 0 and 1 <= e["period"] <= 7 and e["block"] in "spdf" for e in els)
    ens = [e["electronegativity"] for e in els if e["electronegativity"] is not None]
    assert all(0.5 < v < 4.1 for v in ens)


@pytest.mark.parametrize("symbol, mass, group, period, en", [
    ("H", 1.008, 1, 1, 2.2), ("C", 12.011, 14, 2, 2.55), ("O", 15.999, 16, 2, 3.44),
    ("Fe", 55.845, 8, 4, 1.83), ("Au", 196.966569, 11, 6, 2.4),
])
def test_spot_values(symbol, mass, group, period, en):
    e = tool.find_element(symbol)
    assert e["atomic_mass"] == pytest.approx(mass, abs=1e-6)
    assert (e["group"], e["period"], e["electronegativity"]) == (group, period, en)


def test_blank_values_are_none():
    assert tool.find_element("Og")["electronegativity"] is None
    assert tool.find_element("Ce")["group"] is None  # Ce-Lu and Th-Lr have no group number
    assert tool.find_element("La")["group"] == 3


# ---------- lookup ----------
@pytest.mark.parametrize("query", ["26", "Fe", "fe", "FE", "iron", "Iron", " iron "])
def test_find_element_variants(query):
    assert tool.find_element(query)["symbol"] == "Fe"


@pytest.mark.parametrize("query, symbol", [("aluminium", "Al"), ("Sulphur", "S"), ("caesium", "Cs")])
def test_british_spellings(query, symbol):
    assert tool.find_element(query)["symbol"] == symbol


def test_unknown_element_suggests_close_match():
    with pytest.raises(ValueError, match="Did you mean Iron"):
        tool.find_element("iorn")
    with pytest.raises(ValueError, match="unknown element"):
        tool.find_element("Xx")


def test_select_filters():
    els = tool.load_elements()
    assert [e["symbol"] for e in tool.select(els, group=17)] == ["F", "Cl", "Br", "I", "At", "Ts"]
    assert [e["symbol"] for e in tool.select(els, period=1)] == ["H", "He"]
    assert [e["symbol"] for e in tool.select(els, block="d", period=4)] == [
        "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn"]
    assert {e["category"] for e in tool.select(els, category="noble")} == {"Noble gases"}


# ---------- formulas ----------
@pytest.mark.parametrize("formula, counts", [
    ("H2O", {"H": 2, "O": 1}),
    ("C9H8O4", {"C": 9, "H": 8, "O": 4}),
    ("Ca(OH)2", {"Ca": 1, "O": 2, "H": 2}),
    ("K4[Fe(CN)6]", {"K": 4, "Fe": 1, "C": 6, "N": 6}),
    ("CuSO4.5H2O", {"Cu": 1, "S": 1, "O": 9, "H": 10}),
    ("CuSO4\u00b75H2O", {"Cu": 1, "S": 1, "O": 9, "H": 10}),
    ("CO", {"C": 1, "O": 1}),
    ("Co", {"Co": 1}),
    ("2H2O", {"H": 4, "O": 2}),
])
def test_parse_formula(formula, counts):
    assert tool.parse_formula(formula) == counts


@pytest.mark.parametrize("bad", ["Ca(OH", "H2O)", "h2o", "", "H2O$"])
def test_parse_formula_errors(bad):
    with pytest.raises(ValueError):
        tool.parse_formula(bad)


def test_molar_masses():
    assert tool.molar_mass("H2O")[0] == pytest.approx(18.015, abs=1e-3)
    assert tool.molar_mass("NaCl")[0] == pytest.approx(58.44, abs=0.01)
    assert tool.molar_mass("C9H8O4")[0] == pytest.approx(180.16, abs=0.01)  # aspirin, matches Day 3
    assert tool.molar_mass("CuSO4.5H2O")[0] == pytest.approx(249.68, abs=0.01)


def test_molar_mass_unknown_element():
    with pytest.raises(ValueError, match="unknown element"):
        tool.molar_mass("Xx2")


# ---------- CLI ----------
def test_cli_card(monkeypatch, capsys):
    assert run(monkeypatch, "Fe") == 0
    out = capsys.readouterr().out
    assert out.startswith("Fe  Iron  (Z = 26)")
    assert "55.845 u" in out and "d-block, Transition metals" in out


def test_cli_card_handles_missing_values(monkeypatch, capsys):
    assert run(monkeypatch, "Og") == 0
    assert "n/a" in capsys.readouterr().out


def test_cli_property(monkeypatch, capsys):
    assert run(monkeypatch, "O", "--property", "en") == 0
    assert capsys.readouterr().out.strip() == "3.44"
    assert run(monkeypatch, "Cu", "Zn", "--property", "mass") == 0
    assert capsys.readouterr().out.splitlines() == ["Cu\t63.546", "Zn\t65.38"]
    assert run(monkeypatch, "Fe", "--property", "nonsense") == 2


def test_cli_filter_table_and_no_match(monkeypatch, capsys):
    assert run(monkeypatch, "--group", 17) == 0
    out = capsys.readouterr().out
    assert out.splitlines()[0].startswith("Z") and "Chlorine" in out
    assert run(monkeypatch, "--group", 99) == 1


def test_cli_filter_combined_with_query(monkeypatch, capsys):
    assert run(monkeypatch, "Fe", "Na", "--block", "d") == 0
    out = capsys.readouterr().out
    assert "Iron" in out and "Sodium" not in out


def test_cli_formula(monkeypatch, capsys):
    assert run(monkeypatch, "--formula", "CuSO4.5H2O") == 0
    out = capsys.readouterr().out
    assert out.startswith("CuSO4.5H2O: 249.677 g/mol")
    assert run(monkeypatch, "--formula", "Ca(OH") == 2


def test_cli_csv_output(tmp_path, monkeypatch):
    out = tmp_path / "halogens.csv"
    assert run(monkeypatch, "--category", "halogens", "-o", out) == 0
    rows = list(csv.DictReader(out.open()))
    assert [r["symbol"] for r in rows] == ["F", "Cl", "Br", "I", "At", "Ts"]
    assert "electron_configuration" in rows[0]
    assert rows[-1]["electronegativity"] == ""


def test_cli_needs_something_to_do(monkeypatch):
    with pytest.raises(SystemExit) as exc:
        run(monkeypatch)
    assert exc.value.code == 2


def test_cli_unknown_element(monkeypatch, capsys):
    assert run(monkeypatch, "iorn") == 2
    assert "Did you mean Iron" in capsys.readouterr().err
