"""Phase 124: 分岐置換基 propan-2-yl 形式 (IUPAC 2013 P-31.1.3.4)"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # isopropyl → propan-2-yl
    ("CC(C)c1ccccc1", "(propan-2-yl)benzene"),
    # sec-butyl → butan-2-yl
    ("CCC(C)c1ccccc1", "(butan-2-yl)benzene"),
    ("CC(CC)c1ccccc1", "(butan-2-yl)benzene"),
    # propan-2-yl attached to chain
    ("CC(C)CC", "2-methylbutane"),
    ("CC(C)CCC", "2-methylpentane"),
    # 2x propan-2-yl (bis notation)
    ("CCCC(C(C)C)(C(C)C)CCC", "4,4-bis(propan-2-yl)heptane"),
    # 回帰: isobutyl (primary attachment) = 2-methylpropyl, unchanged
    ("CC(C)CC(=O)O", "3-methylbutanoic acid"),
    # 回帰: straight-chain alkyl unchanged
    ("CCCc1ccccc1", "propylbenzene"),
    ("CCc1ccccc1", "ethylbenzene"),
])
def test_phase124_branched_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
