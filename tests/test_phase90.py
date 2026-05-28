"""Phase 90: カルバメート + アルケン鎖 (ethenyl carbamate 等)"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ビニル/アルケニル O 側
    ("C=COC(=O)N", "ethenyl carbamate"),
    ("C=CCOC(=O)N", "prop-2-en-1-yl carbamate"),
    # 回帰: 飽和
    ("NC(=O)OC", "methyl carbamate"),
    ("NC(=O)OCC", "ethyl carbamate"),
    ("CNC(=O)OCC", "ethyl N-methylcarbamate"),
    ("CN(C)C(=O)OC", "methyl N,N-dimethylcarbamate"),
])
def test_phase90_carbamate_with_ene(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
