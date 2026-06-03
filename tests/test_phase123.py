"""Phase 123: ジアゾ化合物命名 (IUPAC 2013 P-68.4.1)"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # diazomethane
    ("C=[N+]=[N-]", "diazomethane"),
    ("[N-]=[N+]=C", "diazomethane"),
    # diazoethane
    ("CC=[N+]=[N-]", "diazoethane"),
    ("[N-]=[N+]=CC", "diazoethane"),
    # diazopropane
    ("CCC=[N+]=[N-]", "diazopropane"),
    # diazobutane
    ("CCCC=[N+]=[N-]", "diazobutane"),
    # 回帰: azide (R-N=N=N) は影響なし
    ("CN=[N+]=[N-]", "azidomethane"),
    # 回帰: azo compound は影響なし
    ("c1ccc(N=Nc2ccccc2)cc1", "azobenzene"),
])
def test_phase123_diazo(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
