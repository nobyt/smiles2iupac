"""Phase 116: カルボジイミドのN,N'- 表記 (IUPAC P-66.5)"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("CN=C=NC", "N,N'-dimethylmethanediimine"),
    ("CCN=C=NCC", "N,N'-diethylmethanediimine"),
    ("CCCN=C=NCCC", "N,N'-dipropylmethanediimine"),
    ("CN=C=NCC", "N-ethyl-N'-methylmethanediimine"),
])
def test_phase116_methanediimine_nn_prime(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
