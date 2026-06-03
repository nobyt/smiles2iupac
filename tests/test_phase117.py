"""Phase 117: IUPAC 2013 保留名 (ethanol, acetic acid, acetate esters)"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ethanol (IUPAC 2013 PIN: P-63.1.1)
    ("CCO", "ethanol"),
    # acetic acid (IUPAC 2013 PIN: P-65.1.1)
    ("CC(=O)O", "acetic acid"),
    # acetate esters
    ("COC(=O)C", "methyl acetate"),
    ("CCOC(=O)C", "ethyl acetate"),
    ("CC(=O)Oc1ccccc1", "phenyl acetate"),
    # 置換エタノール (ethanol 保留名ベース)
    ("NCCO", "2-aminoethanol"),
    ("CCO", "ethanol"),
    # 置換酢酸 (acetic acid 保留名ベース)
    ("c1ccc(CC(=O)O)cc1", "phenylacetic acid"),
    # 回帰: 3炭素以上はそのまま
    ("CCCO", "propan-1-ol"),
    ("CCC(=O)O", "propanoic acid"),
    ("CCCOC(=O)C", "propyl acetate"),
])
def test_phase117_retained_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
