"""Phase 514: sulfanylmethyl substituent naming (thiol-bearing alkyl branches)
(IUPAC 2013 P-63.6: sulfanyl substituent prefixes for -SH on chain carbons).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # -CH2-SH branch on longer chain (suffix correctly wrapped in parens)
    ("CCC(CS)C(=O)O",    "2-(sulfanylmethyl)butanoic acid"),
    ("OCC(CS)CO",         "2-(sulfanylmethyl)propane-1,3-diol"),
    ("NCC(CS)C(=O)O",    "3-amino-2-(sulfanylmethyl)propanoic acid"),
    # Equal-length chain: prefer more substituents (2-methyl-3-sulfanyl wins)
    ("CC(CS)C(=O)O",     "2-methyl-3-sulfanylpropanoic acid"),
    # Thioether on branch → branched substituent name
    ("CCC(SC)C(=O)O",    "2-(methylsulfanyl)butanoic acid"),
    # Regression: plain thiols on main chain unaffected
    ("SCCC",             "propane-1-thiol"),
    ("CC(S)CC",          "butane-2-thiol"),
])
def test_phase514_sulfanylmethyl(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
