"""Phase 122: ジニトリル命名 (IUPAC 2013 P-66.5.1)"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ethanedinitrile (oxalonitrile)
    ("N#CC#N", "ethanedinitrile"),
    # propanedinitrile (malononitrile)
    ("N#CCC#N", "propanedinitrile"),
    # butanedinitrile (succinonitrile)
    ("N#CCCC#N", "butanedinitrile"),
    # pentanedinitrile (glutaronitrile)
    ("N#CCCCC#N", "pentanedinitrile"),
    # hexanedinitrile (adiponitrile)
    ("N#CCCCCC#N", "hexanedinitrile"),
    # 回帰: 単一ニトリル
    ("N#CC", "acetonitrile"),
    ("N#CCCC", "butanenitrile"),
])
def test_phase122_dinitrile(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
