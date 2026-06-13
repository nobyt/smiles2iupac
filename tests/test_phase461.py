"""Phase 461: benzo-fused 5-membered thiadiazole/oxadiazole/isothiazole
(IUPAC 2013 P-31.1.3 retained names).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("c1ccc2snnc2c1",   "1,2,3-benzothiadiazole"),
    ("c1ccc2nsnc2c1",   "2,1,3-benzothiadiazole"),
    ("c1ccc2nonc2c1",   "2,1,3-benzoxadiazole"),
    ("c1ccc2nscc2c1",   "1,2-benzisothiazole"),
])
def test_phase461_benzo_thiadiazole_isothiazole(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
