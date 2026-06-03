"""Phase 218: cyclic dione/trione terminal 'e' NOT elided (IUPAC 2013 P-31.1.3.4).

The terminal 'e' is elided only before suffixes starting with a vowel (a, e, i, o, u).
'-dione' and '-trione' begin with 'd' (consonant) → no elision.
'-one' begins with 'o' (vowel) → elision IS correct.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # imide series (N flanked by two C=O) — no elision before -dione
    ("O=C1NC(=O)CC1",   "pyrrolidine-2,5-dione"),
    ("O=C1NC(=O)CCC1",  "piperidine-2,6-dione"),
    ("O=C1NC(=O)CCCC1", "azepane-2,7-dione"),
    # cyclic anhydride series (O flanked by two C=O) — no elision
    ("O=C1CC(=O)O1",    "oxetane-2,4-dione"),
    ("O=C1CCC(=O)O1",   "oxolane-2,5-dione"),
    ("O=C1CCCC(=O)O1",  "oxane-2,6-dione"),
    # other heterocyclic diones — no elision
    ("O=C1NC(=O)CN1",   "imidazolidine-2,4-dione"),
    ("O=C1CCC(=O)N1",   "pyrrolidine-2,5-dione"),
    # trione — no elision
    ("O=C1NC(=O)NC(=O)N1", "1,3,5-triazinane-2,4,6-trione"),
    # regression: single C=O lactam — elision IS correct ('-one' starts with 'o')
    ("O=C1CCCN1",  "pyrrolidin-2-one"),
    ("O=C1CCCCN1", "piperidin-2-one"),
    ("O=C1CCCCCN1", "azepan-2-one"),
    # regression: single C=O lactone — elision IS correct
    ("O=C1CCCO1",  "oxolan-2-one"),
    ("O=C1CCCCO1", "oxan-2-one"),
    # regression: acyclic dione unchanged
    ("CC(=O)CC(=O)C", "pentane-2,4-dione"),
    # regression: carbocyclic dione unchanged
    ("O=C1CCCCC(=O)C1", "cycloheptane-1,3-dione"),
])
def test_phase218_dione_elision(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
