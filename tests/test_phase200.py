"""Phase 200: monophosphinic acid (1C), sulfinamide, phosphinous acid."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1. phosphinic acid with 1 C substituent
    ("CP(=O)(O)",   "methylphosphinic acid"),
    ("CCP(=O)(O)",  "ethylphosphinic acid"),
    # 2. sulfinamide
    ("CS(=O)N",     "methanesulfinamide"),
    ("CS(=O)NC",    "N-methylmethanesulfinamide"),
    # 3. phosphinous acid
    ("CP(O)",       "methylphosphinous acid"),
    ("CCP(O)",      "ethylphosphinous acid"),
    # regression: dimethylphosphinic acid (Phase 143) still works
    ("CP(=O)(O)C",  "dimethylphosphinic acid"),
])
def test_phase200(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
