"""Command-line interface for smiles2iupac."""

import argparse
import sys

from . import smiles_to_iupac


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="smiles2iupac",
        description="Convert a SMILES string to an IUPAC 2013 preferred name.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
examples:
  smiles2iupac "CC"
  smiles2iupac "c1ccccc1"
  smiles2iupac "CC(=O)O"
  smiles2iupac "C[C@@H](O)CC(=O)O"

  # read multiple SMILES from stdin (one per line):
  printf 'CC\\nCCC\\nCCCC' | smiles2iupac -

  # pipe from a file:
  smiles2iupac - < smiles.txt
""",
    )
    parser.add_argument(
        "smiles",
        nargs="?",
        metavar="SMILES",
        help='SMILES string to convert. Use "-" to read one SMILES per line from stdin.',
    )
    parser.add_argument(
        "--version",
        action="version",
        version="smiles2iupac 0.1.0",
    )

    args = parser.parse_args()

    if args.smiles is None:
        parser.print_help(sys.stderr)
        sys.exit(1)

    if args.smiles == "-":
        _process_stdin()
    else:
        _convert_one(args.smiles)


def _convert_one(smiles: str) -> None:
    try:
        result = smiles_to_iupac(smiles)
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        sys.exit(1)
    print(result)


def _process_stdin() -> None:
    exit_code = 0
    for line in sys.stdin:
        smiles = line.rstrip("\n")
        if not smiles:
            continue
        try:
            print(smiles_to_iupac(smiles))
        except Exception as exc:
            print(f"error: {exc}", file=sys.stderr)
            exit_code = 1
    sys.exit(exit_code)


if __name__ == "__main__":
    main()
