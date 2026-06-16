"""Phase 548: 1,3,5-trithiane and substituted fused heterocycle naming.
1,3,5-trithiane: 6-membered ring with 3 S atoms.
Substituted fused heterocycles (e.g. 2-phenylchromone = flavone): sub-cluster
ring matching allows correct parent identification.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1,3,5-trithiane
    ("C1SCSCS1",   "1,3,5-trithiane"),
    # substituted chromone (flavone)
    ("O=c1cc(-c2ccccc2)oc2ccccc12",  "2-phenylchromone"),
    ("O=C1C=C(c2ccccc2)Oc2ccccc21",  "2-phenylchromone"),
    # substituted flavanone
    ("O=C1CC(c2ccccc2)Oc2ccccc21",   "2-phenylchroman-4-one"),
    # methylchromone regression
    ("Cc1cc(=O)c2ccccc2o1",          "2-methylchromone"),
])
def test_phase548_trithiane_and_substituted_fused(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
