"""Scaffold a new day: python tools/new_day.py 2 smiles-to-3d"""
import re
import shutil
import sys
from pathlib import Path

root = Path(__file__).resolve().parent.parent
n, slug = int(sys.argv[1]), sys.argv[2]
dest = root / "days" / f"day{n:02d}-{slug}"
if dest.exists():
    sys.exit(f"{dest} already exists")
shutil.copytree(root / "template" / "dayNN-name", dest)

title = slug.replace("-", " ").title()
for f in dest.rglob("*"):
    if f.suffix in {".md", ".py"}:
        f.write_text(f.read_text().replace("Day NN", f"Day {n}").replace("<Tool name>", title))

readme = root / "README.md"
text = readme.read_text()
row = f"| {n} | {title} | [`days/{dest.name}/`](days/{dest.name}/) |\n"
text = text.replace("\n## Contributing", row + "\n## Contributing", 1) if row not in text else text
text = re.sub(r"streak-day%20\d+%2F100", f"streak-day%20{n}%2F100", text)
readme.write_text(text)
print(f"Created {dest}")
