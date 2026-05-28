"""Phase 73: カルボジイミド (R-N=C=N-R → N,N'-di{alkyl}carbodiimide)"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 対称カルボジイミド (Phase 116: N,N'- 表記)
    ("CN=C=NC", "N,N'-dimethylcarbodiimide"),
    ("CCN=C=NCC", "N,N'-diethylcarbodiimide"),
    ("CCCN=C=NCCC", "N,N'-dipropylcarbodiimide"),
    # 非対称カルボジイミド
    ("CN=C=NCC", "N-ethyl-N'-methylcarbodiimide"),
])
def test_phase73_carbodiimide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
