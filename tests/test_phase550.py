"""Phase 550: Substituted isocoumarin naming.
Isocoumarin (1H-2-benzopyran-1-one / 1H-isochromen-1-one) with substituents
at positions 3, 4, 5-8. Hydroxy/amino stay as prefixes (PCG already encoded).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # unsubstituted
    ("O=c1occc2ccccc12",      "isocoumarin"),
    # 3-substituted
    ("O=c1oc(C)cc2ccccc12",   "3-methylisocoumarin"),
    ("O=c1oc(O)cc2ccccc12",   "3-hydroxyisocoumarin"),
    # 4-substituted
    ("O=c1occ(C)c2ccccc12",   "4-methylisocoumarin"),
    ("O=c1occ(O)c2ccccc12",   "4-hydroxyisocoumarin"),
    # 3,4-disubstituted
    ("O=c1oc(C)c(C)c2ccccc12", "3,4-dimethylisocoumarin"),
])
def test_phase550_substituted_isocoumarin(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
