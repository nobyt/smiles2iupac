"""Phase 114: 2炭素鎖単一置換基ロカント省略 (IUPAC 2013 P-14.5)"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ハロゲン置換エタン (locant 省略)
    ("CCCl", "chloroethane"),
    ("CCF", "fluoroethane"),
    ("CCBr", "bromoethane"),
    ("CCI", "iodoethane"),
    # ニトロエタン
    ("CC[N+](=O)[O-]", "nitroethane"),
    # エーテル (methoxy/ethoxy on ethane)
    ("CCOC", "methoxyethane"),
    ("CCOCC", "ethoxyethane"),
    # 1,1-二置換: ロカント保持
    ("CC(Cl)Cl", "1,1-dichloroethane"),
    ("CC(F)F", "1,1-difluoroethane"),
    # 1,2-二置換: ロカント保持
    ("ClCCCl", "1,2-dichloroethane"),
    # 3炭素以上: ロカント保持
    ("CCC[N+](=O)[O-]", "1-nitropropane"),
    ("CCCCl", "1-chloropropane"),
])
def test_phase114_ethane_locant_omission(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
