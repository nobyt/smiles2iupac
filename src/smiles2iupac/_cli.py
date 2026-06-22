"""Command-line interface for smiles2iupac."""

import argparse
import sys

from . import smiles_to_iupac
from .pubchem import is_valid_iupac, get_iupac_by_smiles_using_inchikey


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
    parser.add_argument(
        "--validate",
        action="store_true",
        help="Check the generated IUPAC name against PubChem and fail if not matched",
    )

    args = parser.parse_args()

    if args.smiles is None:
        parser.print_help(sys.stderr)
        sys.exit(1)

    if args.smiles == "-":
        _process_stdin(args.validate)
    else:
        _convert_one(args.smiles, validate=args.validate)


def _convert_one(smiles: str, validate: bool = False) -> None:
    try:
        result = smiles_to_iupac(smiles)
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        sys.exit(1)
    print(result)
    if validate:
        try:
            pub_name = get_iupac_by_smiles_using_inchikey(smiles)
            ok = False if pub_name is None else (" ".join(pub_name.strip().lower().split()) == " ".join(result.strip().lower().split()))
        except Exception:
            ok = False
        if not ok:
            print("error: generated IUPAC name was not matched by PubChem (via InChIKey)", file=sys.stderr)
            sys.exit(2)


def _process_stdin(validate: bool = False) -> None:
    exit_code = 0
    for line in sys.stdin:
        smiles = line.rstrip("\n")
        if not smiles:
            continue
        try:
            name = smiles_to_iupac(smiles)
            print(name)
            if validate:
                try:
                    pub_name = get_iupac_by_smiles_using_inchikey(smiles)
                    ok = False if pub_name is None else (" ".join(pub_name.strip().lower().split()) == " ".join(name.strip().lower().split()))
                except Exception:
                    ok = False
                if not ok:
                    print(f"error: generated IUPAC name for '{smiles}' was not matched by PubChem (via InChIKey)", file=sys.stderr)
                    exit_code = 2
        except Exception as exc:
            print(f"error: {exc}", file=sys.stderr)
            exit_code = 1
    sys.exit(exit_code)


if __name__ == "__main__":
    main()
