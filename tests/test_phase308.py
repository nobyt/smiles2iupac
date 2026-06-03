"""Phase 308: diselenoamide naming (IUPAC 2013 P-65.1.2.4).

Two selenoamide groups on a chain → "ethanediselenoamide" etc.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("NC(=[Se])C(=[Se])N",   "ethanediselenoamide"),
    ("NC(=[Se])CC(=[Se])N",  "propanediselenoamide"),
    # regression: mono selenoamide unchanged
    ("CC(=[Se])N",           "ethaneselenoamide"),
])
def test_phase308_diselenoamide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
