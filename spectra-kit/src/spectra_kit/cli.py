"""The `spectra` command: one entry point for every tool in the kit.

    spectra list                       the tools
    spectra uvvis spectrum.csv -o p.png     run a tool (same options as the stand-alone spectra-uvvis)
    spectra examples [DIR]             copy the synthetic example files to try the tools on
"""
import importlib
import sys
from importlib import resources
from pathlib import Path

from spectra_kit import __version__
from spectra_kit.registry import BY_KEY, TOOLS, command

DEFAULT_EXAMPLES_DIR = "spectra-examples"

USAGE = """usage: spectra <tool> [options]
       spectra list | examples [DIR] | --version

Run `spectra <tool> --help` for the options of a tool. Tools:
"""


def tool_listing() -> str:
    width = max(len(t.key) for t in TOOLS)
    return "\n".join(f"  {t.key:<{width}}  Day {t.day:<2}  {t.summary}  ({command(t)})" for t in TOOLS)


def example_files() -> list[tuple[str, str]]:
    """(tool key, file name) for every example file that ships with the package."""
    return [(t.key, name) for t in TOOLS for name in t.examples]


def copy_examples(dest: Path, force: bool = False) -> list[Path]:
    """Copy the packaged example files to dest/<tool>/<file>. Refuses to overwrite unless force is set."""
    root = resources.files("spectra_kit") / "examples"
    plan = [(dest / key / name, root / key / name) for key, name in example_files()]
    existing = [target for target, _ in plan if target.exists()]
    if existing and not force:
        raise FileExistsError(f"{len(existing)} file(s) already exist, for example {existing[0]}. Use --force to overwrite")
    written = []
    for target, source in plan:
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(source.read_bytes())
        written.append(target)
    return written


def run_tool(key: str, args: list[str]) -> int:
    """Run one tool's main() with its own command line."""
    module = importlib.import_module(f"spectra_kit.{BY_KEY[key].module}")
    saved = sys.argv
    sys.argv = [f"spectra {key}", *args]
    try:
        return int(module.main())
    finally:
        sys.argv = saved


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in ("-h", "--help"):
        print(USAGE + tool_listing())
        return 0
    command_name, rest = argv[0], argv[1:]
    if command_name in ("--version", "-V"):
        print(f"spectra-kit {__version__}")
        return 0
    if command_name == "list":
        print(tool_listing())
        return 0
    if command_name == "examples":
        force = "--force" in rest
        rest = [a for a in rest if a != "--force"]
        if len(rest) > 1:
            print("error: usage: spectra examples [DIR] [--force]", file=sys.stderr)
            return 2
        try:
            written = copy_examples(Path(rest[0] if rest else DEFAULT_EXAMPLES_DIR), force)
        except FileExistsError as err:
            print(f"error: {err}", file=sys.stderr)
            return 2
        print(f"Copied {len(written)} example files to {rest[0] if rest else DEFAULT_EXAMPLES_DIR}")
        return 0
    if command_name not in BY_KEY:
        print(f"error: unknown tool {command_name!r}. Choose one of: {', '.join(BY_KEY)}", file=sys.stderr)
        return 2
    return run_tool(command_name, rest)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except BrokenPipeError:  # e.g. piping into `head`
        sys.exit(0)
