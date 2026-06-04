"""Phase 383: Terminal nitrile N named as 'cyano' substituent (IUPAC 2013 P-65.2.1.1).

When -C≡N is at the end of a chain and COOH is the principal group,
the terminal N was returned as '(N)' instead of 'cyano'.

Fix: _name_nitrogen_substituent now detects N triple-bonded to C (terminal)
and returns 'cyano'.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # nitrile as substituent (COOH principal)
    ("N#CCC(=O)O",    "3-cyanopropanoic acid"),
    ("N#CCCCC(=O)O",  "5-cyanopentanoic acid"),
    ("OC(=O)CC#N",    "3-cyanopropanoic acid"),
    # cyano on branched chain
    ("N#CC(C)C(=O)O", "3-cyano-2-methylpropanoic acid"),
    # nitrile alone still works
    ("N#CCCC",        "butanenitrile"),
    ("N#CC",          "acetonitrile"),
])
def test_phase383_cyano_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
