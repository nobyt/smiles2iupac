"""Phase 291: no elision of terminal 'e' before '-dione'/'-trione' (IUPAC 2013 P-31.1.3.4).

The terminal 'e' of a parent name is elided only before suffixes beginning with
a vowel (a, e, i, o, u).  '-dione' and '-trione' begin with 'd' (consonant),
so the 'e' is retained:
  pyrrolidine-2,5-dione  (not pyrrolidin-2,5-dione)
  piperidine-2,6-dione   (not piperidin-2,6-dione)
  oxetane-2,4-dione      (not oxetan-2,4-dione)

Contrast: '-one' begins with 'o' → elision IS correct:
  pyrrolidin-2-one, piperidin-2-one
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ── dione: terminal 'e' retained ─────────────────────────────────────
    ("O=C1NC(=O)CC1",    "pyrrolidine-2,5-dione"),
    ("O=C1CCC(=O)N1",    "pyrrolidine-2,5-dione"),
    ("O=C1NC(=O)CCC1",   "piperidine-2,6-dione"),
    ("O=C1CCCC(=O)N1",   "piperidine-2,6-dione"),
    ("O=C1CC(=O)O1",     "oxetane-2,4-dione"),

    # ── -one: elision IS correct ('o' is a vowel) ────────────────────────
    ("O=C1CCCN1",        "pyrrolidin-2-one"),
    ("O=C1CCCCN1",       "piperidin-2-one"),

    # ── regressions: all-carbon diones unchanged ─────────────────────────
    ("O=C1CCCC1=O",      "cyclopentane-1,2-dione"),
    ("O=C1CC(=O)CCC1",   "cyclohexane-1,3-dione"),
])
def test_phase291_dione_elision(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
