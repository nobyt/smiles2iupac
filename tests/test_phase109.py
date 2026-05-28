"""Phase 109: N-置換基括弧表記 + 複合アリールアミド命名 (IUPAC P-14.5.7, P-66.4)"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ヘテロアリール N-置換 → 括弧必要
    ("CC(=O)Nc1ccncc1", "N-(pyridin-4-yl)acetamide"),
    ("CC(=O)Nc1cccs1", "N-(thiophen-2-yl)acetamide"),
    # 置換フェニル N-置換 → 括弧必要
    ("CC(=O)Nc1ccc(Cl)cc1", "N-(4-chlorophenyl)acetamide"),
    ("CC(=O)Nc1ccc(C)cc1", "N-(4-methylphenyl)acetamide"),
    # ヘテロアリール on ベンズアミド
    ("c1ccc(C(=O)Nc2ccncc2)cc1", "N-(pyridin-4-yl)benzamide"),
    # 非置換フェニル → 括弧不要
    ("CC(=O)Nc1ccccc1", "N-phenylacetamide"),
    # 単純アルキル → 括弧不要
    ("CC(=O)NC", "N-methylacetamide"),
    ("CC(=O)NCC", "N-ethylacetamide"),
    # 回帰: pyridine-3-carboxylic acid (pyridine 上の外環 COOH)
    ("c1ccc(cn1)C(=O)O", "pyridine-3-carboxylic acid"),
    # 回帰: N-phenylacetamide (通常の N-aryl amide)
    ("c1ccc(NC(=O)C)cc1", "N-phenylacetamide"),
])
def test_phase109_n_substituent_parentheses(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
