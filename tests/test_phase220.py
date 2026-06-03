"""Phase 220: 1,2-oxazetidine (4-membered ring with adjacent O and N).

IUPAC 2013: 4-membered O-N heterocycle where O and N are at positions 1 and 2.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1,2-oxazetidine: 4-membered ring, O-N adjacent
    ("C1CNO1",  "1,2-oxazetidine"),
    # regression: 1,3-oxazetidine (4-membered, O and N non-adjacent) still works
    ("O1CNC1",  "1,3-oxazetidine"),
])
def test_phase220_oxazetidine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
