"""Phase 843: drop indicated H from names when N is substituted (not NH).

IUPAC 2013 P-14.7: indicated H notation (nH- prefix or (nH) locant) is used
only for H atoms at sp3 positions in an otherwise maximally unsaturated ring.
When N instead carries a substituent (e.g. methyl, phenyl), no H is present
at that position, so the indicated H notation must be omitted.

Two patterns fixed:
1. isoindole-1,3(2H)-dione (phthalimide): drop (2H) when N2 is substituted.
2. 1H-pyrrole-2,5-dione (maleimide): drop 1H- when N1 is substituted.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # isoindole-1,3(2H)-dione: N2-H → keep (2H)
    ("O=C1NC(=O)c2ccccc21",            "isoindole-1,3(2H)-dione"),
    # N2-methyl → drop (2H)
    ("CN1C(=O)c2ccccc2C1=O",           "2-methylisoindole-1,3-dione"),
    # N2-ethyl → drop (2H)
    ("CCN1C(=O)c2ccccc2C1=O",          "2-ethylisoindole-1,3-dione"),
    # benzo-ring substituents: N2 still has H → keep (2H)
    ("O=C1NC(=O)c2c(C)cccc21",         "4-methylisoindole-1,3(2H)-dione"),
    ("O=C1NC(=O)c2cc(C)ccc21",         "5-methylisoindole-1,3(2H)-dione"),
    # 1H-pyrrole-2,5-dione (maleimide): N1-H → keep 1H-
    ("O=C1C=CC(=O)N1",                 "1H-pyrrole-2,5-dione"),
    # N1-phenyl → drop 1H-
    ("O=C1C=CC(=O)N1c1ccccc1",         "1-phenylpyrrole-2,5-dione"),
])
def test_phase843(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
