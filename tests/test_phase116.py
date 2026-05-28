"""Phase 116: カルボジイミドのN,N'- 表記 (IUPAC P-66.5)"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("CN=C=NC", "N,N'-dimethylcarbodiimide"),
    ("CCN=C=NCC", "N,N'-diethylcarbodiimide"),
    ("CCCN=C=NCCC", "N,N'-dipropylcarbodiimide"),
    ("CN=C=NCC", "N-ethyl-N'-methylcarbodiimide"),
])
def test_phase116_carbodiimide_nn_prime(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
