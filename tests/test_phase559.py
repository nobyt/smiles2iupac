"""Phase 559: Substituted acridine naming.
Acridine has C2v symmetry (C1≡C8, C2≡C7, C3≡C6, C4≡C5); preferred IUPAC
name uses the lower locant. C9 is the unique central-ring carbon gamma to N10.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # unsubstituted
    ("c1ccc2nc3ccccc3cc2c1",    "acridine"),
    # C9 (unique central-ring position adjacent to N)
    ("Cc1c2ccccc2nc2ccccc12",   "9-methylacridine"),
    # outer-ring positions (lower locant wins due to C2v symmetry)
    # Phase 854: corrected 1<->4, 2<->3 (locant map ring-swap fix)
    ("Cc1cccc2cc3ccccc3nc12",   "4-methylacridine"),
    ("Cc1ccc2cc3ccccc3nc2c1",   "3-methylacridine"),
    ("Cc1ccc2nc3ccccc3cc2c1",   "2-methylacridine"),
    ("Cc1cccc2nc3ccccc3cc12",   "1-methylacridine"),
    # hydroxy at C9 converts to lactam (Phase 761: gamma tautomer)
    ("Oc1c2ccccc2nc2ccccc12",   "acridin-9(10H)-one"),
    ("Nc1c2ccccc2nc2ccccc12",   "acridin-9-amine"),
])
def test_phase559_acridine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
