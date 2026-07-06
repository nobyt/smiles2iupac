"""Phase 840: pyrazino[2,3-b]pyrazine — parent and methyl derivatives."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("c1cnc2nccnc2n1",    "pyrazino[2,3-b]pyrazine"),
    ("Cc1cnc2nccnc2n1",   "2-methylpyrazino[2,3-b]pyrazine"),
])
def test_phase840(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
