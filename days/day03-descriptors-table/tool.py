"""Day 3: Descriptors Table. One-line description."""
import argparse


def main(path: str) -> None:
    print(f"Processing {path}")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("path", help="input file")
    main(p.parse_args().path)
