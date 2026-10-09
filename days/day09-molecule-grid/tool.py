"""Day 9: molecule grid image generator for figures.

Draw a list of SMILES as a labelled grid image (PNG or SVG) for papers, theses and
slides. Optionally highlight and align a substructure, add atom numbers, and size the
image for a given figure width and DPI.

Usage:
    python tool.py example/input.smi -o grid.png
    python tool.py example/input.smi -o grid.svg --cols 3 --legend formula
    python tool.py example/input.smi -o ring.png --highlight "c1ccccc1" --align
    python tool.py example/input.smi -o fig1.png --cols 3 --width-in 6.5 --dpi 300
"""
import argparse
import sys
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Descriptors, Draw, rdDepictor, rdMolDescriptors

LEGENDS = ["name", "smiles", "formula", "mw", "none"]


def read_entries(arg: str) -> list[tuple[str, str]]:
    """Read (smiles, name) pairs from a .smi file, or treat arg as one SMILES."""
    src = Path(arg)
    try:
        is_file = src.is_file()
    except OSError:  # a very long SMILES can exceed the filename length limit
        is_file = False
    if not is_file:
        return [(arg, "molecule")]
    entries = []
    for i, line in enumerate(src.read_text().splitlines(), 1):
        parts = line.split()
        if parts and not parts[0].startswith("#"):
            entries.append((parts[0], parts[1] if len(parts) > 1 else f"mol{i}"))
    return entries


def make_legend(mol, smiles: str, name: str, kind: str) -> str:
    if kind == "name":
        return name
    if kind == "smiles":
        return smiles
    if kind == "formula":
        return rdMolDescriptors.CalcMolFormula(mol)
    if kind == "mw":
        return f"{name}  {Descriptors.MolWt(mol):.1f}"
    return ""


def grid_size(cols: int, size: int, width_in: float | None, dpi: int) -> int:
    """Pixels per molecule: either --size, or derived from a figure width (inches) and DPI."""
    if width_in is None:
        return size
    if width_in <= 0 or dpi <= 0:
        raise ValueError("--width-in and --dpi must be positive")
    return max(50, round(width_in * dpi / cols))


def pattern_from_smarts(smarts: str):
    patt = Chem.MolFromSmarts(smarts)
    if patt is None:
        raise ValueError(f"invalid SMARTS: {smarts!r}")
    return patt


def highlight_lists(mols, patt):
    """Atom and bond highlight lists for the first match of patt in each molecule."""
    atoms, bonds, hits = [], [], 0
    for m in mols:
        match = m.GetSubstructMatch(patt)
        if match:
            hits += 1
        inside = set(match)
        atoms.append(list(match))
        bonds.append([b.GetIdx() for b in m.GetBonds()
                      if b.GetBeginAtomIdx() in inside and b.GetEndAtomIdx() in inside])
    return atoms, bonds, hits


def align_to_pattern(mols, patt) -> int:
    """Redraw molecules so the matched substructure has the same orientation. Returns count aligned."""
    template = Chem.Mol(patt)
    rdDepictor.Compute2DCoords(template)
    done = 0
    for m in mols:
        if m.HasSubstructMatch(patt):
            rdDepictor.GenerateDepictionMatching2DStructure(m, template, acceptFailure=True)
            done += 1
    return done


def build_image(mols, legends, cols, size, *, svg=False, highlight=None, bw=False, indices=False):
    opts = Draw.MolDrawOptions()
    opts.legendFontSize = max(12, size // 12)
    opts.addAtomIndices = indices
    if bw:
        opts.useBWAtomPalette()
    kwargs = {}
    if highlight is not None:
        kwargs["highlightAtomLists"], kwargs["highlightBondLists"] = highlight
    return Draw.MolsToGridImage(mols, molsPerRow=cols, subImgSize=(size, size), legends=legends,
                                useSVG=svg, drawOptions=opts, **kwargs)


def save(image, path: Path, svg: bool, dpi: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if svg:
        path.write_text(image if isinstance(image, str) else image.data)
    else:
        image.save(path, dpi=(dpi, dpi))


def page_paths(out: Path, pages: int) -> list[Path]:
    if pages == 1:
        return [out]
    return [out.with_name(f"{out.stem}_{i}{out.suffix}") for i in range(1, pages + 1)]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", help="SMILES string or .smi file (one SMILES per line, optional name)")
    ap.add_argument("-o", "--output", default="grid.png", help="output .png or .svg (default: grid.png)")
    ap.add_argument("--cols", type=int, default=3, help="molecules per row (default 3)")
    ap.add_argument("--size", type=int, default=300, help="pixels per molecule (default 300)")
    ap.add_argument("--width-in", type=float, help="total figure width in inches (sets the size from --dpi)")
    ap.add_argument("--dpi", type=int, default=300, help="resolution for --width-in and PNG metadata (default 300)")
    ap.add_argument("--legend", choices=LEGENDS, default="name", help="label under each molecule (default: name)")
    ap.add_argument("--per-page", type=int, help="molecules per image; extra pages get _1, _2, ... suffixes")
    ap.add_argument("--highlight", metavar="SMARTS", help="highlight this substructure in every molecule")
    ap.add_argument("--align", action="store_true", help="orient molecules by the --highlight substructure")
    ap.add_argument("--atom-indices", action="store_true", help="show atom numbers")
    ap.add_argument("--bw", action="store_true", help="black and white atoms (journal figures)")
    args = ap.parse_args()

    out = Path(args.output)
    svg = out.suffix.lower() == ".svg"
    if out.suffix.lower() not in (".png", ".svg"):
        print("error: output must end in .png or .svg", file=sys.stderr)
        return 2
    if args.cols < 1 or (args.per_page is not None and args.per_page < 1):
        print("error: --cols and --per-page must be at least 1", file=sys.stderr)
        return 2
    if args.align and not args.highlight:
        print("error: --align needs --highlight SMARTS", file=sys.stderr)
        return 2

    mols, legends, bad = [], [], 0
    for smi, name in read_entries(args.input):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            print(f"[skip] {name}: invalid SMILES {smi!r}", file=sys.stderr)
            bad += 1
            continue
        mols.append(mol)
        legends.append(make_legend(mol, smi, name, args.legend))
    if not mols:
        print("error: no valid molecules to draw", file=sys.stderr)
        return 2

    try:
        patt = pattern_from_smarts(args.highlight) if args.highlight else None
        size = grid_size(min(args.cols, len(mols)), args.size, args.width_in, args.dpi)
    except ValueError as err:
        print(f"error: {err}", file=sys.stderr)
        return 2

    if patt is not None:
        hits = sum(m.HasSubstructMatch(patt) for m in mols)
        print(f"highlight: pattern found in {hits} of {len(mols)} molecules")
        if args.align:
            print(f"align: oriented {align_to_pattern(mols, patt)} molecules")

    per_page = args.per_page or len(mols)
    chunks = [slice(i, i + per_page) for i in range(0, len(mols), per_page)]
    for path, sl in zip(page_paths(out, len(chunks)), chunks):
        page = mols[sl]
        highlight = highlight_lists(page, patt)[:2] if patt is not None else None
        image = build_image(page, legends[sl], min(args.cols, len(page)), size,
                            svg=svg, highlight=highlight, bw=args.bw, indices=args.atom_indices)
        save(image, path, svg, args.dpi)
        print(f"[ok] {len(page)} molecules -> {path}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
