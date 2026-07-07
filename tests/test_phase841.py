"""Phase 841: imidazol-2(1H)-one parent and 4-methyl derivative."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("O=C1NC=CN1",    "imidazol-2(1H)-one"),
    ("CC1=CNC(=O)N1", "4-methylimidazol-2(1H)-one"),
])
def test_phase841(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
