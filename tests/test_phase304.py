"""Phase 304: hydrazinecarboxylic acid naming (IUPAC 2013 P-65.1.2.3).

H2N-NH-C(=O)-OH → hydrazinecarboxylic acid.
N-substituted variants follow 2-alkyl locant convention.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("NNC(=O)O",       "hydrazinecarboxylic acid"),
    # regressions: carbamic acid unchanged
    ("NC(=O)O",        "carbamic acid"),
    ("CNC(=O)O",       "N-methylcarbamic acid"),
])
def test_phase304_hydrazinecarboxylic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
