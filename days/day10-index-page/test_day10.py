from pathlib import Path

import pytest
import tool

GOOD_README = "# Day {n}: {title}\n\n**What:** Does a thing | with a pipe.\n**Why:** Because.\n"


def make_day(root: Path, n: int, slug: str, title="Thing", readme=None, tool_py=True, tests=True, example=True):
    folder = root / "days" / f"day{n:02d}-{slug}"
    folder.mkdir(parents=True)
    if readme is not False:
        (folder / "README.md").write_text(readme if readme is not None else GOOD_README.format(n=n, title=title))
    if tool_py:
        (folder / "tool.py").write_text("print('x')\n")
    if tests:
        (folder / "test_x.py").write_text("def test_x():\n    pass\n")
    if example:
        (folder / "example").mkdir()
        (folder / "example" / "in.txt").write_text("x")
    return folder


def make_repo(tmp_path, badge=2, link_all=True):
    make_day(tmp_path, 2, "first", "First")
    make_day(tmp_path, 3, "second", "Second")
    links = " ".join(f"[d](days/{p.name}/)" for p in sorted((tmp_path / "days").iterdir())) if link_all else ""
    (tmp_path / "README.md").write_text(f"![streak](https://img.shields.io/badge/streak-day%20{badge}%2F100-blue)\n{links}\n")
    return tmp_path


def run(monkeypatch, *argv):
    monkeypatch.setattr("sys.argv", ["tool.py", *map(str, argv)])
    return tool.main()


# ---------- scanning one day ----------
def test_valid_day(tmp_path):
    day = tool.scan_day(make_day(tmp_path, 2, "first", "First"))
    assert (day.number, day.title, day.tests) == (2, "First", True)
    assert day.what == "Does a thing | with a pipe."
    assert day.path == "days/day02-first/"
    assert day.errors == [] and day.warnings == []


@pytest.mark.parametrize("kwargs, fragment", [
    ({"tool_py": False}, "missing tool.py"),
    ({"readme": False}, "missing README.md"),
    ({"readme": "no heading\n**What:** a\n**Why:** b\n"}, "must start with"),
    ({"readme": "# Day 5: Wrong\n**What:** a\n**Why:** b\n"}, "folder is day02"),
    ({"readme": "# Day 2: T\n**Why:** b\n"}, "'**What:**'"),
    ({"readme": "# Day 2: T\n**What:** a\n"}, "'**Why:**'"),
    ({"readme": "# Day 2: T\n**What:** one sentence on what it does.\n**Why:** b\n"}, "placeholder"),
])
def test_errors(tmp_path, kwargs, fragment):
    day = tool.scan_day(make_day(tmp_path, 2, "first", **kwargs))
    assert any(fragment in e for e in day.errors), day.errors


def test_warnings_for_missing_example_and_tests(tmp_path):
    day = tool.scan_day(make_day(tmp_path, 2, "first", tests=False, example=False))
    assert day.errors == []
    assert sorted(w.split(": ")[1] for w in day.warnings) == ["no example/ files", "no test_*.py"]
    assert day.tests is False


def test_duplicate_day_numbers(tmp_path):
    make_day(tmp_path, 2, "a")
    make_day(tmp_path, 2, "b")
    errors = [e for d in tool.scan_days(tmp_path) for e in d.errors]
    assert any("also used by" in e for e in errors)


def test_scan_days_sorted_and_ignores_other_folders(tmp_path):
    make_day(tmp_path, 10, "z")
    make_day(tmp_path, 2, "a")
    (tmp_path / "days" / "notes").mkdir()
    (tmp_path / "days" / "day-x").mkdir()
    assert [d.number for d in tool.scan_days(tmp_path)] == [2, 10]


# ---------- rendering ----------
def test_clip():
    assert tool.clip("short") == "short"
    assert tool.clip("a|b") == "a\\|b"
    long = "word " * 60
    out = tool.clip(long, 40)
    assert out.endswith("...") and len(out) <= 43 and not out.startswith(" ")


def test_render_index(tmp_path):
    make_repo(tmp_path)
    make_day(tmp_path, 12, "later", "Later", tests=False)
    text = tool.render_index(tool.scan_days(tmp_path))
    assert text.count("## Phase") == 2
    assert "## Phase 1 (Days 1-10): Foundations and template" in text
    assert "## Phase 2 (Days 11-20): Nanomaterial and spectroscopy data" in text
    assert "## Phase 3" not in text
    assert "| 1 | [Repo setup](template/) |" in text
    assert "| 2 | [First](days/day02-first/) | Does a thing \\| with a pipe. | yes |" in text
    assert "| 12 | [Later](days/day12-later/) | Does a thing \\| with a pipe. | - |" in text
    assert "**4 of 100 days done.**" in text
    assert text == tool.render_index(tool.scan_days(tmp_path))  # deterministic


# ---------- main README check ----------
def test_main_readme_ok(tmp_path):
    make_repo(tmp_path, badge=3)
    assert tool.check_main_readme(tmp_path, tool.scan_days(tmp_path)) == []


def test_main_readme_problems(tmp_path):
    make_repo(tmp_path, badge=2, link_all=False)
    errors = tool.check_main_readme(tmp_path, tool.scan_days(tmp_path))
    assert any("does not link days/day03-second/" in e for e in errors)
    assert any("badge says day 2" in e for e in errors)
    (tmp_path / "README.md").write_text("nothing\n")
    assert any("no streak badge" in e for e in tool.check_main_readme(tmp_path, tool.scan_days(tmp_path)))
    (tmp_path / "README.md").unlink()
    assert tool.check_main_readme(tmp_path, []) == ["README.md is missing"]


# ---------- CLI ----------
def test_cli_write_then_check(tmp_path, monkeypatch, capsys):
    make_repo(tmp_path, badge=3)
    assert run(monkeypatch, "--root", tmp_path) == 0
    index = tmp_path / "INDEX.md"
    assert index.is_file()
    assert run(monkeypatch, "--root", tmp_path, "--check") == 0
    assert "checked 2 days: 0 errors" in capsys.readouterr().out


def test_cli_check_detects_stale_and_missing_index(tmp_path, monkeypatch, capsys):
    make_repo(tmp_path, badge=3)
    assert run(monkeypatch, "--root", tmp_path, "--check") == 1  # no INDEX.md yet
    assert not (tmp_path / "INDEX.md").exists()  # --check never writes
    run(monkeypatch, "--root", tmp_path)
    (tmp_path / "days" / "day02-first" / "README.md").write_text(
        "# Day 2: First\n\n**What:** Changed text.\n**Why:** Because.\n")
    assert run(monkeypatch, "--root", tmp_path, "--check") == 1
    assert "out of date" in capsys.readouterr().out


def test_cli_reports_errors_with_exit_1(tmp_path, monkeypatch, capsys):
    make_repo(tmp_path, badge=3)
    make_day(tmp_path, 4, "bad", readme="oops")
    assert run(monkeypatch, "--root", tmp_path) == 1
    assert "[error] day04-bad: README must start with" in capsys.readouterr().out


def test_cli_no_days_folder(tmp_path, monkeypatch, capsys):
    assert run(monkeypatch, "--root", tmp_path) == 2
    assert "no days/ folder" in capsys.readouterr().err


# ---------- this repository ----------
def test_this_repo_passes_its_own_rules():
    days = tool.scan_days(tool.ROOT)
    assert days, "no day folders found"
    assert [e for d in days for e in d.errors] == []
    assert tool.check_main_readme(tool.ROOT, days) == []
    index = tool.ROOT / "INDEX.md"
    assert index.is_file() and index.read_text() == tool.render_index(days), \
        "INDEX.md is out of date: run python days/day10-index-page/tool.py"


# ---------- agent files ----------
SKILL = """---
name: demo-skill
description: Does a thing. Use when asked to do the thing.
license: MIT
metadata:
  version: "1.0"
---
# Demo skill
"""


def add_skill(root: Path, text: str = SKILL, folder: str = ".claude/skills", name: str = "demo-skill") -> Path:
    path = root / folder / name / "SKILL.md"
    path.parent.mkdir(parents=True)
    path.write_text(text)
    return path


def test_parse_frontmatter_plain_folded_and_metadata():
    fields, meta = tool.parse_frontmatter(SKILL)
    assert fields["name"] == "demo-skill" and fields["license"] == "MIT" and meta == {"version": '"1.0"'}
    folded = "---\nname: x\ndescription: >\n  first part\n  second part\n---\n"
    assert tool.parse_frontmatter(folded)[0]["description"] == "first part second part"
    assert tool.parse_frontmatter("no frontmatter here") is None
    assert tool.parse_frontmatter("---\nname: x\n") is None  # never closed


def test_valid_skill_has_no_errors(tmp_path):
    assert tool.check_skill(add_skill(tmp_path), tmp_path) == []


@pytest.mark.parametrize("text, fragment", [
    ("# no frontmatter\n", "missing YAML frontmatter"),
    (SKILL.replace("license: MIT", "license: MIT\nowner: me"), "unknown top-level field 'owner'"),
    (SKILL.replace("name: demo-skill\n", ""), "missing name"),
    (SKILL.replace("name: demo-skill", "name: Demo-Skill"), "lowercase letters"),
    (SKILL.replace("name: demo-skill", "name: demo--skill"), "lowercase letters"),
    (SKILL.replace("name: demo-skill", "name: " + "a" * 65), "lowercase letters"),
    (SKILL.replace("name: demo-skill", "name: other-name"), "must match the folder name"),
    (SKILL.replace("description: Does a thing. Use when asked to do the thing.\n", ""), "missing description"),
    (SKILL.replace("Does a thing.", "x" * 1030), "limit is 1024"),
    (SKILL.replace("metadata:\n  version: \"1.0\"", 'metadata: {version: "1.0"}'), "block mapping"),
    (SKILL.replace('version: "1.0"', "version: 1.0"), "quoted number"),
])
def test_invalid_skills(tmp_path, text, fragment):
    errors = tool.check_skill(add_skill(tmp_path, text), tmp_path)
    assert any(fragment in e for e in errors), errors


def test_check_agent_files(tmp_path):
    errors, warnings = tool.check_agent_files(tmp_path)
    assert errors == [] and any("AGENTS.md is missing" in w for w in warnings)
    (tmp_path / "AGENTS.md").write_text("# Agents\n")
    add_skill(tmp_path)
    add_skill(tmp_path, SKILL.replace("demo-skill", "other-skill"), folder=".agents/skills", name="other-skill")
    assert tool.check_agent_files(tmp_path) == ([], [])
    add_skill(tmp_path, "broken", name="bad-skill")
    errors, _ = tool.check_agent_files(tmp_path)
    assert len(errors) == 1 and ".claude/skills/bad-skill/SKILL.md" in errors[0]


def test_cli_agent_file_problems(tmp_path, monkeypatch, capsys):
    make_repo(tmp_path, badge=3)
    assert run(monkeypatch, "--root", tmp_path) == 0  # a missing AGENTS.md is only a warning
    assert "AGENTS.md is missing" in capsys.readouterr().out
    add_skill(tmp_path, "broken")
    assert run(monkeypatch, "--root", tmp_path, "--check") == 1
    assert "missing YAML frontmatter" in capsys.readouterr().out


def test_this_repo_agent_files_are_valid():
    errors, warnings = tool.check_agent_files(tool.ROOT)
    assert errors == [] and warnings == []
