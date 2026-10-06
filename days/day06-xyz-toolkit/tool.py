"""Day 6: XYZ file toolkit.

Small, dependency-free helpers for the structure files you juggle before and
after QM runs: inspect, split, merge, center, align and convert
(xyz, pdb, Turbomole coord).

Usage:
    python tool.py info example/ethanol.xyz
    python tool.py split trajectory.xyz --outdir frames
    python tool.py merge example/ethanol.xyz example/benzene.xyz -o all.xyz
    python tool.py center example/ethanol.xyz --mode mass -o ethanol_com.xyz
    python tool.py align example/ethanol.xyz shifted.xyz -o aligned.xyz --rmsd
    python tool.py convert example/ethanol.xyz -o ethanol.pdb
"""
import argparse
import math
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

BOHR_TO_ANG = 0.529177210903  # CODATA 2018

# Standard atomic weights (g/mol), enough for most organic, organometallic and QD work
MASSES = {
    "H": 1.008, "He": 4.0026, "Li": 6.94, "Be": 9.0122, "B": 10.81, "C": 12.011,
    "N": 14.007, "O": 15.999, "F": 18.998, "Ne": 20.180, "Na": 22.990, "Mg": 24.305,
    "Al": 26.982, "Si": 28.085, "P": 30.974, "S": 32.06, "Cl": 35.45, "Ar": 39.948,
    "K": 39.098, "Ca": 40.078, "Sc": 44.956, "Ti": 47.867, "V": 50.942, "Cr": 51.996,
    "Mn": 54.938, "Fe": 55.845, "Co": 58.933, "Ni": 58.693, "Cu": 63.546, "Zn": 65.38,
    "Ga": 69.723, "Ge": 72.630, "As": 74.922, "Se": 78.971, "Br": 79.904, "Kr": 83.798,
    "Rb": 85.468, "Sr": 87.62, "Y": 88.906, "Zr": 91.224, "Mo": 95.95, "Ru": 101.07,
    "Rh": 102.91, "Pd": 106.42, "Ag": 107.87, "Cd": 112.41, "Sn": 118.71, "Sb": 121.76,
    "Te": 127.60, "I": 126.90, "Xe": 131.29, "Pt": 195.08, "Au": 196.97, "Hg": 200.59,
    "Pb": 207.2,
}


@dataclass
class Frame:
    symbols: list[str]
    coords: list[list[float]]  # angstrom
    comment: str = field(default="")


def norm_symbol(token: str) -> str:
    if not re.fullmatch(r"[A-Za-z]{1,2}", token):
        raise ValueError(f"not an element symbol: {token!r}")
    return token[0].upper() + token[1:].lower()


# ---------- readers ----------
def read_xyz(text: str) -> list[Frame]:
    lines, i, frames = text.splitlines(), 0, []
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        try:
            n = int(lines[i].split()[0])
        except ValueError:
            raise ValueError(f"line {i + 1}: expected an atom count, got {lines[i].strip()!r}") from None
        if n <= 0 or i + 2 + n > len(lines):
            raise ValueError(f"line {i + 1}: frame declares {n} atoms but the file does not have them")
        symbols, coords = [], []
        for j in range(n):
            ln = i + 3 + j
            parts = lines[i + 2 + j].split()
            if len(parts) < 4:
                raise ValueError(f"line {ln}: expected 'symbol x y z'")
            try:
                symbols.append(norm_symbol(parts[0]))
                coords.append([float(v) for v in parts[1:4]])
            except ValueError as err:
                raise ValueError(f"line {ln}: {err}") from None
        frames.append(Frame(symbols, coords, lines[i + 1].strip()))
        i += 2 + n
    if not frames:
        raise ValueError("no structures found")
    return frames


def _pdb_element(line: str) -> str:
    if len(line) >= 78 and line[76:78].strip():
        return norm_symbol(line[76:78].strip())
    letters = re.sub(r"[^A-Za-z]", "", line[12:16])
    if not letters:
        raise ValueError("cannot determine element")
    return norm_symbol(letters[:1] if line[12:13] == " " else letters[:2])


def read_pdb(text: str) -> list[Frame]:
    frames, cur = [], Frame([], [])
    for k, line in enumerate(text.splitlines(), 1):
        rec = line[:6].strip()
        if rec in ("ATOM", "HETATM"):
            try:
                cur.symbols.append(_pdb_element(line))
                cur.coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            except ValueError as err:
                raise ValueError(f"line {k}: {err}") from None
        elif rec == "ENDMDL" and cur.symbols:
            frames.append(cur)
            cur = Frame([], [])
    if cur.symbols:
        frames.append(cur)
    if not frames:
        raise ValueError("no ATOM/HETATM records found")
    return frames


def read_coord(text: str) -> list[Frame]:
    """Turbomole/xtb 'coord' file: coordinates in bohr, lowercase symbols."""
    frame, inside = Frame([], []), False
    for k, line in enumerate(text.splitlines(), 1):
        s = line.strip()
        if s.startswith("$coord"):
            inside = True
        elif s.startswith("$"):
            inside = False
        elif inside and s:
            p = s.split()
            try:
                frame.coords.append([float(v) * BOHR_TO_ANG for v in p[:3]])
                frame.symbols.append(norm_symbol(p[3]))
            except (ValueError, IndexError):
                raise ValueError(f"line {k}: expected 'x y z symbol'") from None
    if not frame.symbols:
        raise ValueError("no $coord block found")
    return [frame]


# ---------- writers ----------
def write_xyz(frames: list[Frame]) -> str:
    out = []
    for f in frames:
        out += [str(len(f.symbols)), f.comment]
        out += [f"{s:<2} {x:12.6f} {y:12.6f} {z:12.6f}" for s, (x, y, z) in zip(f.symbols, f.coords)]
    return "\n".join(out) + "\n"


def write_pdb(frames: list[Frame]) -> str:
    out = []
    for m, f in enumerate(frames, 1):
        if len(frames) > 1:
            out.append(f"MODEL     {m:>4}")
        for i, (s, (x, y, z)) in enumerate(zip(f.symbols, f.coords), 1):
            nm = f"{s}{i}"[:4]
            name = (" " + nm).ljust(4)[:4] if len(s) == 1 else nm.ljust(4)
            out.append(f"HETATM{i % 100000:5d} {name} MOL A{1:4d}    "
                       f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {s:>2}")
        out.append("ENDMDL" if len(frames) > 1 else "TER")
    out.append("END")
    return "\n".join(out) + "\n"


def write_coord(frames: list[Frame]) -> str:
    if len(frames) != 1:
        raise ValueError("the coord format holds one structure; split the file first")
    rows = [f"{x / BOHR_TO_ANG:20.14f}{y / BOHR_TO_ANG:20.14f}{z / BOHR_TO_ANG:20.14f}      {s.lower()}"
            for s, (x, y, z) in zip(frames[0].symbols, frames[0].coords)]
    return "\n".join(["$coord", *rows, "$end"]) + "\n"


READERS = {"xyz": read_xyz, "pdb": read_pdb, "coord": read_coord}
WRITERS = {"xyz": write_xyz, "pdb": write_pdb, "coord": write_coord}


def file_format(path: Path) -> str:
    ext = path.suffix.lower().lstrip(".")
    if ext in READERS:
        return ext
    if path.name.lower() == "coord":
        return "coord"
    raise ValueError(f"unsupported file type: {path.name} (use .xyz, .pdb or .coord)")


def read_frames(path: Path) -> list[Frame]:
    try:
        return READERS[file_format(path)](path.read_text())
    except ValueError as err:
        raise ValueError(f"{path}: {err}") from None
    except OSError as err:
        raise ValueError(f"cannot read {path}: {err.strerror}") from None


def write_frames(frames: list[Frame], path: Path) -> None:
    text = WRITERS[file_format(path)](frames)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


# ---------- geometry ----------
def formula(frame: Frame) -> str:
    counts: dict[str, int] = {}
    for s in frame.symbols:
        counts[s] = counts.get(s, 0) + 1
    if "C" in counts:
        order = ["C"] + (["H"] if "H" in counts else []) + sorted(e for e in counts if e not in ("C", "H"))
    else:
        order = sorted(counts)
    return "".join(f"{e}{counts[e] if counts[e] > 1 else ''}" for e in order)


def centroid(frame: Frame) -> list[float]:
    n = len(frame.coords)
    return [sum(c[k] for c in frame.coords) / n for k in range(3)]


def center_of_mass(frame: Frame) -> list[float]:
    unknown = sorted({s for s in frame.symbols if s not in MASSES})
    if unknown:
        raise ValueError(f"no atomic mass for: {', '.join(unknown)}")
    masses = [MASSES[s] for s in frame.symbols]
    total = sum(masses)
    return [sum(m * c[k] for m, c in zip(masses, frame.coords)) / total for k in range(3)]


def centered(frame: Frame, mode: str = "centroid") -> Frame:
    origin = center_of_mass(frame) if mode == "mass" else centroid(frame)
    coords = [[c[k] - origin[k] for k in range(3)] for c in frame.coords]
    return Frame(list(frame.symbols), coords, frame.comment)


def bounding_box(frame: Frame) -> list[float]:
    return [max(c[k] for c in frame.coords) - min(c[k] for c in frame.coords) for k in range(3)]


def _dot(a: list[float], b: list[float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a: list[float], b: list[float]) -> list[float]:
    return [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]


def _norm(v: list[float]) -> float:
    return math.sqrt(_dot(v, v))


def _normalized(v: list[float], eps: float = 1e-12) -> list[float] | None:
    n = _norm(v)
    if n <= eps:
        return None
    return [x / n for x in v]


def _transpose(m: list[list[float]]) -> list[list[float]]:
    return [[m[j][i] for j in range(3)] for i in range(3)]


def _mat_mul(a: list[list[float]], b: list[list[float]]) -> list[list[float]]:
    return [[sum(a[i][k] * b[k][j] for k in range(3)) for j in range(3)] for i in range(3)]


def _vec_mat_mul(v: list[float], m: list[list[float]]) -> list[float]:
    return [sum(v[k] * m[k][j] for k in range(3)) for j in range(3)]


def _det3(m: list[list[float]]) -> float:
    return (
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    )


def _identity3() -> list[list[float]]:
    return [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]


def _jacobi_eigen_symmetric3(a: list[list[float]], max_iter: int = 50, eps: float = 1e-12) -> tuple[list[float], list[list[float]]]:
    """Eigenvalues and eigenvectors (columns) for a symmetric 3x3 matrix."""
    m = [row[:] for row in a]
    v = _identity3()
    for _ in range(max_iter):
        p, q = 0, 1
        off = abs(m[0][1])
        for i, j in ((0, 2), (1, 2)):
            cur = abs(m[i][j])
            if cur > off:
                p, q, off = i, j, cur
        if off < eps:
            break
        app, aqq, apq = m[p][p], m[q][q], m[p][q]
        phi = 0.5 * math.atan2(2.0 * apq, aqq - app)
        c, s = math.cos(phi), math.sin(phi)
        for k in range(3):
            mkp, mkq = m[k][p], m[k][q]
            m[k][p], m[k][q] = c * mkp - s * mkq, s * mkp + c * mkq
        for k in range(3):
            mpk, mqk = m[p][k], m[q][k]
            m[p][k], m[q][k] = c * mpk - s * mqk, s * mpk + c * mqk
        m[p][q] = 0.0
        m[q][p] = 0.0
        for k in range(3):
            vkp, vkq = v[k][p], v[k][q]
            v[k][p], v[k][q] = c * vkp - s * vkq, s * vkp + c * vkq
    eigvals = [m[i][i] for i in range(3)]
    order = sorted(range(3), key=lambda i: eigvals[i], reverse=True)
    vals = [eigvals[i] for i in order]
    vecs = [[v[row][i] for i in order] for row in range(3)]
    return vals, vecs


def _orthonormal_columns(cols: list[list[float]]) -> list[list[float]]:
    basis: list[list[float]] = []
    candidates = cols + [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    for cand in candidates:
        w = cand[:]
        for b in basis:
            proj = _dot(w, b)
            w = [w[i] - proj * b[i] for i in range(3)]
        nw = _normalized(w)
        if nw is not None:
            basis.append(nw)
        if len(basis) == 3:
            break
    return basis


def _kabsch_rotation(mobile_centered: list[list[float]], reference_centered: list[list[float]]) -> list[list[float]]:
    cov = [[0.0, 0.0, 0.0] for _ in range(3)]
    for m, r in zip(mobile_centered, reference_centered):
        for i in range(3):
            for j in range(3):
                cov[i][j] += m[i] * r[j]
    ct_c = _mat_mul(_transpose(cov), cov)
    eigvals, v = _jacobi_eigen_symmetric3(ct_c)
    sigmas = [math.sqrt(max(ev, 0.0)) for ev in eigvals]
    v_cols = [[v[row][i] for row in range(3)] for i in range(3)]
    u_cols: list[list[float]] = []
    for sigma, v_col in zip(sigmas, v_cols):
        if sigma <= 1e-12:
            continue
        u_col = _normalized([sum(cov[i][k] * v_col[k] for k in range(3)) / sigma for i in range(3)])
        if u_col is not None:
            u_cols.append(u_col)
    u_cols = _orthonormal_columns(u_cols)
    if len(u_cols) != 3:
        raise ValueError("cannot determine a stable rotation for these coordinates")
    v_cols = _orthonormal_columns(v_cols)
    u = [[u_cols[col][row] for col in range(3)] for row in range(3)]
    vv = [[v_cols[col][row] for col in range(3)] for row in range(3)]
    rot = _mat_mul(u, _transpose(vv))
    if _det3(rot) < 0.0:
        u_cols[2] = [-x for x in u_cols[2]]
        u = [[u_cols[col][row] for col in range(3)] for row in range(3)]
        rot = _mat_mul(u, _transpose(vv))
    if _det3(rot) < 0.0:
        raise ValueError("cannot determine a proper rotation for these coordinates")
    return rot


def rmsd(frame_a: Frame, frame_b: Frame) -> float:
    if len(frame_a.coords) != len(frame_b.coords):
        raise ValueError("RMSD requires the same atom count")
    if not frame_a.coords:
        raise ValueError("RMSD requires at least one atom")
    total = 0.0
    for ca, cb in zip(frame_a.coords, frame_b.coords):
        total += sum((ca[k] - cb[k]) ** 2 for k in range(3))
    return math.sqrt(total / len(frame_a.coords))


def aligned_to(reference: Frame, mobile: Frame) -> Frame:
    if len(reference.symbols) != len(mobile.symbols):
        raise ValueError("frames differ in atom count")
    if reference.symbols != mobile.symbols:
        raise ValueError("frames differ in element order; align matching atoms in the same order")
    if not reference.coords:
        raise ValueError("cannot align an empty structure")
    for c in reference.coords + mobile.coords:
        if not all(math.isfinite(v) for v in c):
            raise ValueError("coordinates must be finite numbers")
    cref = centroid(reference)
    cmob = centroid(mobile)
    ref_centered = [[c[k] - cref[k] for k in range(3)] for c in reference.coords]
    mob_centered = [[c[k] - cmob[k] for k in range(3)] for c in mobile.coords]
    rot = _kabsch_rotation(mob_centered, ref_centered)
    coords = []
    for c in mobile.coords:
        shifted = [c[k] - cmob[k] for k in range(3)]
        turned = _vec_mat_mul(shifted, rot)
        coords.append([turned[k] + cref[k] for k in range(3)])
    return Frame(list(mobile.symbols), coords, mobile.comment)


# ---------- commands ----------
def cmd_info(args) -> None:
    frames = read_frames(Path(args.input))
    print(f"file: {args.input}\nframes: {len(frames)}")
    for n, f in enumerate(frames[:3], 1):
        cx, cy, cz = centroid(f)
        print(f"[frame {n}] atoms: {len(f.symbols)}  formula: {formula(f)}")
        print(f"  centroid (A):       {cx:10.4f} {cy:10.4f} {cz:10.4f}")
        try:
            mx, my, mz = center_of_mass(f)
            print(f"  center of mass (A): {mx:10.4f} {my:10.4f} {mz:10.4f}")
        except ValueError as err:
            print(f"  center of mass (A): n/a ({err})")
        dx, dy, dz = bounding_box(f)
        print(f"  bounding box (A):   {dx:10.4f} {dy:10.4f} {dz:10.4f}")
        if f.comment:
            print(f"  comment: {f.comment}")
    if len(frames) > 3:
        print(f"... and {len(frames) - 3} more frames")


def cmd_split(args) -> None:
    src = Path(args.input)
    frames = read_frames(src)
    width = max(3, len(str(len(frames))))
    prefix = args.prefix or src.stem
    for n, f in enumerate(frames, 1):
        out = Path(args.outdir) / f"{prefix}_{n:0{width}d}.{args.format}"
        write_frames([f], out)
        print(f"[ok] frame {n}: {len(f.symbols)} atoms -> {out}")


def cmd_merge(args) -> None:
    frames: list[Frame] = []
    for p in args.inputs:
        frames += read_frames(Path(p))
    if args.strict and any(f.symbols != frames[0].symbols for f in frames):
        raise ValueError("--strict: frames differ in atom count or element order")
    write_frames(frames, Path(args.output))
    print(f"[ok] merged {len(frames)} structures from {len(args.inputs)} files -> {args.output}")


def cmd_center(args) -> None:
    src = Path(args.input)
    frames = [centered(f, args.mode) for f in read_frames(src)]
    out = Path(args.output) if args.output else src.with_name(f"{src.stem}_centered{src.suffix}")
    write_frames(frames, out)
    what = "center of mass" if args.mode == "mass" else "centroid"
    print(f"[ok] moved the {what} of {len(frames)} structure(s) to the origin -> {out}")


def cmd_convert(args) -> None:
    frames = read_frames(Path(args.input))
    write_frames(frames, Path(args.output))
    print(f"[ok] {args.input} ({file_format(Path(args.input))}) -> {args.output} ({file_format(Path(args.output))})")


def _pick_frame(frames: list[Frame], which: str, idx1: int) -> Frame:
    if idx1 <= 0:
        raise ValueError(f"{which} frame index must be >= 1")
    if idx1 > len(frames):
        raise ValueError(f"{which} file has {len(frames)} frame(s), cannot use frame {idx1}")
    return frames[idx1 - 1]


def cmd_align(args) -> None:
    ref_path = Path(args.reference)
    mob_path = Path(args.mobile)
    ref_frame = _pick_frame(read_frames(ref_path), "reference", args.reference_frame)
    mob_frame = _pick_frame(read_frames(mob_path), "mobile", args.mobile_frame)
    aligned = aligned_to(ref_frame, mob_frame)
    out = Path(args.output) if args.output else mob_path.with_name(f"{mob_path.stem}_aligned{mob_path.suffix}")
    write_frames([aligned], out)
    before = rmsd(ref_frame, mob_frame)
    after = rmsd(ref_frame, aligned)
    print(f"[ok] aligned {mob_path} frame {args.mobile_frame} to {ref_path} frame {args.reference_frame} -> {out}")
    if args.rmsd:
        print(f"  RMSD before alignment (A): {before:.6f}")
        print(f"  RMSD after alignment (A):  {after:.6f}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="command", required=True)

    p = sub.add_parser("info", help="atoms, formula, centroid, center of mass, size")
    p.add_argument("input")
    p.set_defaults(func=cmd_info)

    p = sub.add_parser("split", help="split a multi-structure file into one file per structure")
    p.add_argument("input")
    p.add_argument("--outdir", default="frames")
    p.add_argument("--prefix", help="file name prefix (default: input file name)")
    p.add_argument("--format", choices=sorted(WRITERS), default="xyz")
    p.set_defaults(func=cmd_split)

    p = sub.add_parser("merge", help="merge files into one multi-structure file")
    p.add_argument("inputs", nargs="+")
    p.add_argument("-o", "--output", required=True)
    p.add_argument("--strict", action="store_true", help="require identical atoms in every frame (trajectories)")
    p.set_defaults(func=cmd_merge)

    p = sub.add_parser("center", help="move each structure to the origin")
    p.add_argument("input")
    p.add_argument("-o", "--output", help="default: <input>_centered.<ext>")
    p.add_argument("--mode", choices=["centroid", "mass"], default="centroid")
    p.set_defaults(func=cmd_center)

    p = sub.add_parser("convert", help="convert between .xyz, .pdb and Turbomole coord")
    p.add_argument("input")
    p.add_argument("-o", "--output", required=True)
    p.set_defaults(func=cmd_convert)

    p = sub.add_parser("align", help="align one structure onto a reference with a Kabsch fit")
    p.add_argument("reference")
    p.add_argument("mobile")
    p.add_argument("-o", "--output", help="default: <mobile>_aligned.<ext>")
    p.add_argument("--reference-frame", type=int, default=1, help="1-based frame index in reference file")
    p.add_argument("--mobile-frame", type=int, default=1, help="1-based frame index in mobile file")
    p.add_argument("--rmsd", action="store_true", help="print RMSD before and after alignment")
    p.set_defaults(func=cmd_align)

    args = ap.parse_args()
    try:
        args.func(args)
    except ValueError as err:
        print(f"error: {err}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
