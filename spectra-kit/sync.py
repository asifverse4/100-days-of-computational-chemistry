"""Keep spectra-kit in step with the daily tools.

    python spectra-kit/sync.py          # rewrite the package files from days/
    python spectra-kit/sync.py --check  # exit 1 if any copy is out of date (CI runs this)

Each tool.py is copied into the package with its module docstring rewritten, so `--help` shows the installed
command name (spectra-uvvis) and the packaged example paths instead of `python tool.py` and `example/`.
Nothing else in the code changes. The example files and the LICENSE are copied byte for byte.
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
KIT = ROOT / "spectra-kit"
PACKAGE = KIT / "src" / "spectra_kit"
sys.path.insert(0, str(KIT / "src"))  # use the registry in this checkout, not an installed copy

from spectra_kit.registry import TOOLS, Tool, command


def render_module(tool: Tool) -> str:
    """The package copy of a daily tool: a generated-file header and an adjusted module docstring."""
    text = (ROOT / tool.source / "tool.py").read_text()
    match = re.match(r'"""(.*?)"""', text, re.DOTALL)
    if not match:
        raise SystemExit(f"{tool.source}/tool.py must start with a module docstring")
    doc = match[1].replace("python tool.py", command(tool)).replace("example/", f"spectra-examples/{tool.key}/")
    header = f"# Generated from {tool.source}/tool.py by spectra-kit/sync.py. Do not edit this copy.\n"
    return header + '"""' + doc + '"""' + text[match.end():]


def expected_files() -> dict[Path, bytes]:
    """Every generated file and its expected content."""
    files: dict[Path, bytes] = {KIT / "LICENSE": (ROOT / "LICENSE").read_bytes()}
    for tool in TOOLS:
        files[PACKAGE / f"{tool.module}.py"] = render_module(tool).encode()
        for name in tool.examples:
            files[PACKAGE / "examples" / tool.key / name] = (ROOT / tool.source / "example" / name).read_bytes()
    return files


def same(a: bytes, b: bytes) -> bool:
    """Equal apart from line endings (a Windows checkout may turn LF into CRLF)."""
    return a.replace(b"\r\n", b"\n") == b.replace(b"\r\n", b"\n")


def stale_files(files: dict[Path, bytes]) -> list[str]:
    problems = []
    for path, content in files.items():
        if not path.is_file():
            problems.append(f"missing: {path.relative_to(ROOT).as_posix()}")
        elif not same(path.read_bytes(), content):
            problems.append(f"out of date: {path.relative_to(ROOT).as_posix()}")
    examples = PACKAGE / "examples"
    if examples.is_dir():
        problems += [f"not in the registry: {p.relative_to(ROOT).as_posix()}"
                     for p in sorted(examples.rglob("*")) if p.is_file() and p not in files]
    return problems


def main() -> int:
    check = "--check" in sys.argv[1:]
    files = expected_files()
    if check:
        problems = stale_files(files)
        for p in problems:
            print(f"[error] {p}")
        if problems:
            print("spectra-kit is out of step with the daily tools. Run: python spectra-kit/sync.py")
            return 1
        print(f"spectra-kit is in step with the daily tools ({len(files)} files)")
        return 0
    for path, content in files.items():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(content.replace(b"\r\n", b"\n"))
    for stale in stale_files(files):
        print(f"[warn] {stale}")
    print(f"Wrote {len(files)} files to {KIT.relative_to(ROOT).as_posix()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
