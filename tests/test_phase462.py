"""Phase 462: benzo[c]thiophene, benzo[b/c]selenophene, benzoselenadiazoles
(IUPAC 2013 P-31.1.3 retained names).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("c1ccc2cscc2c1",    "benzo[c]thiophene"),
    ("c1ccc2[se]ccc2c1", "benzo[b]selenophene"),
    ("c1ccc2c[se]cc2c1", "benzo[c]selenophene"),
    ("c1ccc2[se]nnc2c1", "1,2,3-benzoselenadiazole"),
    ("c1ccc2n[se]nc2c1", "2,1,3-benzoselenadiazole"),
])
def test_phase462_benzo_selenophene_selenadiazole(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
