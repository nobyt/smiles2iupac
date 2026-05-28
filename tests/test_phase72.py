"""Phase 72: 酸ハライド / エステルの鎖上置換基"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 酸ハライドの鎖上置換基
    ("ClCC(=O)Cl", "2-chloroethanoyl chloride"),
    ("BrCC(=O)Cl", "2-bromoethanoyl chloride"),
    ("ClCC(=O)Br", "2-chloroethanoyl bromide"),
    ("OCCC(=O)Cl", "3-hydroxypropanoyl chloride"),
    ("NCCC(=O)Cl", "3-aminopropanoyl chloride"),
    ("ClCCC(=O)Cl", "3-chloropropanoyl chloride"),
    # エステルの鎖上置換基
    ("ClCC(=O)OC", "methyl 2-chloroacetate"),
    ("BrCC(=O)OCC", "ethyl 2-bromoacetate"),
    ("OCCC(=O)OC", "methyl 3-hydroxypropanoate"),
    ("NCCC(=O)OCC", "ethyl 3-aminopropanoate"),
])
def test_phase72_acid_halide_ester_substituents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
