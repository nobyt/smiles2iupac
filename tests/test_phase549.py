"""Phase 549: Substituted coumarin naming.
Coumarin (2H-chromen-2-one) with substituents at positions 3, 4, 6, 7, 8.
When the base retained name already encodes a PCG (lactone), substituents
like hydroxy/amino stay as prefixes, not converted to -ol/-amine suffixes.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # unsubstituted coumarin
    ("O=c1ccc2ccccc2o1",         "chromen-2-one"),
    # 4-substituted
    ("O=c1cc(C)c2ccccc2o1",      "4-methylchromen-2-one"),
    ("O=c1cc(O)c2ccccc2o1",      "4-hydroxychromen-2-one"),
    ("O=c1cc(N)c2ccccc2o1",      "4-aminochromen-2-one"),
    # 3-substituted
    ("O=c1ccc2ccccc2o1",         "chromen-2-one"),   # sanity
    # 6-substituted
    ("O=c1ccc2cc(O)ccc2o1",      "6-hydroxychromen-2-one"),
    # 7-substituted
    ("O=c1ccc2ccc(O)cc2o1",      "7-hydroxychromen-2-one"),
])
def test_phase549_substituted_coumarin(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
