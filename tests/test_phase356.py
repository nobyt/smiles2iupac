"""Phase 356: Gem-disubstituted ring naming (IUPAC 2013).

When two or more substituents share the same ring locant (gem pattern),
IUPAC names require explicit locants such as "1,1-dimethylcyclohexane".
Previously the locant was dropped producing wrong names like "(dimethyl)cyclohexane".
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1,1-disubstituted cyclohexane
    ("CC1(C)CCCCC1",           "1,1-dimethylcyclohexane"),
    ("CC1(CC)CCCCC1",          "1-ethyl-1-methylcyclohexane"),
    ("ClC1(Cl)CCCCC1",         "1,1-dichlorocyclohexane"),
    # Smaller rings
    ("CC1(C)CCCC1",            "1,1-dimethylcyclopentane"),
    ("CC1(C)CCC1",             "1,1-dimethylcyclobutane"),
    ("CC1(C)CC1",              "1,1-dimethylcyclopropane"),
    # Regressions: mono-substituted still has no locant
    ("CC1CCCCC1",              "methylcyclohexane"),
    ("ClC1CCCCC1",             "chlorocyclohexane"),
    # Multi-position stays correct
    ("CC1CCCCC1C",             "1,2-dimethylcyclohexane"),
    ("CC1CCC(C)CC1",           "1,4-dimethylcyclohexane"),
])
def test_phase356_gem_disubstituted_rings(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
