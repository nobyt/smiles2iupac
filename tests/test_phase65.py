"""Phase 65: ハロホルメート拡張 (F, Cl, Br, I → carbonofluoridate/chloridate/bromate/iodate)"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # クロロホルメート (chloroformate → carbonochloridate)
    ("ClC(=O)OC", "methyl carbonochloridate"),
    ("ClC(=O)OCC", "ethyl carbonochloridate"),
    ("ClC(=O)OCCC", "propyl carbonochloridate"),
    ("ClC(=O)OC(C)C", "propan-2-yl carbonochloridate"),
    # フルオロホルメート (carbonofluoridate)
    ("FC(=O)OC", "methyl carbonofluoridate"),
    ("FC(=O)OCC", "ethyl carbonofluoridate"),
    # ブロモホルメート (carbonobromate)
    ("BrC(=O)OCC", "ethyl carbonobromate"),
    ("BrC(=O)OCCC", "propyl carbonobromate"),
    # ヨードホルメート (carbonoiodate)
    ("IC(=O)OCCC", "propyl carbonoiodate"),
    ("IC(=O)OCC", "ethyl carbonoiodate"),
])
def test_phase65_haloformate(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
