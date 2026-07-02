"""Phase 760: 1,2,4-benzotriazin-3-ol → 1,2,4-benzotriazin-3(2H)-one (IUPAC 2013).

Extends tautomeric rules to benzo-fused 1,2,4-triazine:
- 1,2,4-benzotriazin-3-ol → 1,2,4-benzotriazin-3(2H)-one
- 1,2,4-benzotriazin-3-thiol → 1,2,4-benzotriazin-3(2H)-thione
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("Oc1nnc2ccccc2n1",                "1,2,4-benzotriazin-3(2H)-one"),
    ("Sc1nnc2ccccc2n1",                "1,2,4-benzotriazin-3(2H)-thione"),
    # Regression: parent ring unaffected
    ("c1nnc2ccccc2n1",                 "1,2,4-benzotriazine"),
    # Regression: 1,2,3-triazine (Phase 757) unchanged
    ("Oc1ccnnn1",                      "1,2,3-triazin-4(3H)-one"),
    # Regression: 1,2,4-triazine (Phase 756) unchanged
    ("Oc1cncnn1",                      "1,2,4-triazin-6(1H)-one"),
])
def test_phase760_benzotriazine_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
