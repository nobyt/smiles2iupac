"""Phase 110: アリールスルホニル命名 — benzene 直結 S(=O)(=O)/S(=O) (IUPAC P-65.3)"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ベンゼンスルホン酸
    ("c1ccc(S(=O)(=O)O)cc1", "benzenesulfonic acid"),
    # ベンゼンスルホニルクロライド
    ("c1ccc(S(=O)(=O)Cl)cc1", "benzenesulfonyl chloride"),
    # ベンゼンスルホンアミド
    ("c1ccc(S(=O)(=O)N)cc1", "benzenesulfonamide"),
    # N-置換ベンゼンスルホンアミド
    ("c1ccc(S(=O)(=O)NC)cc1", "N-methylbenzenesulfonamide"),
    # ベンゼンスルフィン酸
    ("c1ccc(S(=O)O)cc1", "benzenesulfinic acid"),
    # 4-メチルベンゼンスルホン酸 (tosylic acid)
    ("Cc1ccc(S(=O)(=O)O)cc1", "4-methylbenzenesulfonic acid"),
    # ベンゼンスルホン酸エステル
    ("c1ccc(S(=O)(=O)OCC)cc1", "ethyl benzenesulfonate"),
    # 回帰: アルキルスルホン酸は変わらず
    ("CS(=O)(=O)O", "methanesulfonic acid"),
    ("CS(=O)(=O)Cl", "methanesulfonyl chloride"),
    ("CS(=O)(=O)N", "methanesulfonamide"),
    ("CS(=O)O", "methanesulfinic acid"),
    ("CS(=O)(=O)OC", "methyl methanesulfonate"),
])
def test_phase110_aryl_sulfonyl(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
