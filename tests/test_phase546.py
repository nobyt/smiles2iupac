"""Phase 546: Zwitterionic sulfoxide [S+][O-] notation (IUPAC 2013 P-65.3).
C[S+](C)[O-] is equivalent to CS(=O)C and should give the same name.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("C[S+](C)[O-]",             "dimethyl sulfoxide"),
    ("CC[S+](CC)[O-]",           "diethyl sulfoxide"),
    ("C[S+](CC)[O-]",            "ethyl methyl sulfoxide"),
    ("c1ccc([S+](c2ccccc2)[O-])cc1", "diphenyl sulfoxide"),
])
def test_phase546_zwitterionic_sulfoxide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
