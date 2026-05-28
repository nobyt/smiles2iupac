"""Phase 107: ベンゾイルハライド命名 + N-置換ベンズアミド (IUPAC P-65.1.2.1, P-66.4)"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ベンゾイルハライド (aryl acid halide)
    ("ClC(=O)c1ccccc1", "benzoyl chloride"),
    ("BrC(=O)c1ccccc1", "benzoyl bromide"),
    ("FC(=O)c1ccccc1", "benzoyl fluoride"),
    ("IC(=O)c1ccccc1", "benzoyl iodide"),
    # N-置換ベンズアミド
    ("c1ccc(C(=O)NC)cc1", "N-methylbenzamide"),
    ("c1ccc(C(=O)N(C)C)cc1", "N,N-dimethylbenzamide"),
    ("c1ccc(C(=O)NCC)cc1", "N-ethylbenzamide"),
    # 回帰: 一級ベンズアミド
    ("c1ccc(C(=O)N)cc1", "benzamide"),
    # 回帰: 非芳香族酸ハライド変わらず
    ("CC(=O)Cl", "acetyl chloride"),
    ("CCC(=O)Cl", "propanoyl chloride"),
    # 回帰: アリールケトン (フェニルが置換基側)
    ("O=C(C)c1ccccc1", "acetophenone"),
])
def test_phase107_benzoyl_halide_and_n_aryl_amide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
