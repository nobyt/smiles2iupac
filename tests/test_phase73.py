"""Phase 73: カルボジイミド (R-N=C=N-R → N,N'-di{alkyl}methanediimine)"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 対称カルボジイミド (Phase 116: N,N'- 表記)
    ("CN=C=NC", "N,N'-dimethylmethanediimine"),
    ("CCN=C=NCC", "N,N'-diethylmethanediimine"),
    ("CCCN=C=NCCC", "N,N'-dipropylmethanediimine"),
    # 非対称カルボジイミド
    ("CN=C=NCC", "N-ethyl-N'-methylmethanediimine"),
])
def test_phase73_methanediimine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
