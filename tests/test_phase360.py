"""Phase 360: IUPAC 2013 alternating enclosing marks (P-16.3.4).

When a substituent name itself contains parentheses, the outer enclosure must
use brackets [ ] instead of parentheses ( ), per IUPAC 2013 P-16.3.4:
  Level 1: ( )
  Level 2: [ ]
  Level 3: { }

Examples:
  "dimethyl sulfide"  -- level-1 only, keeps ()
  "N-[(2E)-but-2-en-1-yl]urea"  -- inner (2E) forces outer to []
  "2-[4-(2-methylpropyl)phenyl]propanoic acid"  -- inner () forces outer []

Simple stereo descriptors like "(E)-" or "(1R,2S)-" at the start of the
full name are at level-0 (top), so they keep () and are never upgraded.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Level-1 only: () retained
    ("CS(=O)C",              "dimethyl sulfoxide"),
    ("CSC",                  "dimethyl sulfide"),
    ("OC(C)c1ccccc1",        "1-phenylethanol"),
    # Top-level stereo descriptor: always keeps ()
    ("OC(=O)/C=C/CC",        "(2E)-pent-2-enoic acid"),
    # Level-2: [] when substituent contains ()
    ("C/C=C/CNC(=O)O",       "N-[(2E)-but-2-en-1-yl]carbamic acid"),
    (r"C/C=C\CNC(=O)O",      "N-[(2Z)-but-2-en-1-yl]carbamic acid"),
    ("C/C=C/CNC(=O)N",       "N-[(2E)-but-2-en-1-yl]urea"),
    ("C/C=C/CNc1ccccc1",     "N-[(2E)-but-2-en-1-yl]aniline"),
    ("C/C=C/CNO",            "N-[(2E)-but-2-en-1-yl]hydroxylamine"),
    # Nested substituted phenyl: inner () forces outer []
    ("CC(C)Cc1ccc(cc1)C(C)C(=O)O",  "2-[4-(2-methylpropyl)phenyl]propanoic acid"),
    # Regression: top-level (E)- not converted
    ("C/C=C/C",              "(2E)-but-2-ene"),
    # Regression: simple (methylsulfonyl) not converted
    ("CS(=O)(=O)C",          "dimethyl sulfone"),
    ("CC(O)S(=O)C",          "1-(methylsulfinyl)ethanol"),
])
def test_phase360_alternating_enclosing_marks(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
