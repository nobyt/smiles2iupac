"""Phase 277: enyne numbering — double bonds get lower locants than triple bonds
when combined locant sets are equal (IUPAC 2013 P-31.1.6.3).

When a chain has both C=C and C≡C and both numbering directions give the same
combined locant set, the direction giving the double bond the lower locant is chosen.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ── symmetric locant set → double bond wins ───────────────────────────
    # same molecule as C=CC#C; combined set {1,3} from both ends
    ("C#CC=C",    "but-1-en-3-yne"),
    # same molecule: locant set {1,4} from both ends
    ("C#CCC=C",   "pent-1-en-4-yne"),
    # locant set {1,5} from both ends
    ("C#CCCC=C",  "hex-1-en-5-yne"),

    # ── asymmetric locant set → overall minimum wins (not the db rule) ────
    # C=CCCC#C: from left db=1,tb=5 {1,5}; from right tb=1,db=5 {1,5} → tie → db lower → db=1
    ("C=CCCC#C",  "hex-1-en-5-yne"),
    # C#CC=CC: from C# end: tb=1,db=3 {1,3}; from C end: db=2,tb=4 {2,4} → {1,3} wins
    ("C#CC=CC",   "pent-3-en-1-yne"),

    # ── regressions: pure alkene/alkyne unaffected ────────────────────────
    ("C=CC=C",    "buta-1,3-diene"),
    ("C#CC#C",    "buta-1,3-diyne"),
    ("C=CC#C",    "but-1-en-3-yne"),
])
def test_phase277_enyne_numbering(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
