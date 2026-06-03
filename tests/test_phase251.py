"""Phase 251: O-thioester and S-dithioate ester naming (IUPAC 2013 P-65.1.2.2).

  CC(=S)OR → O-alkyl alkanethioate  (O-ester of thioic O-acid)
  CC(=S)SR → S-alkyl alkanedithioate (S-ester of dithioic acid)
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # O-thioesters
    ("CC(=S)OC",      "O-methyl ethanethioate"),
    ("CC(=S)OCC",     "O-ethyl ethanethioate"),
    ("CCC(=S)OC",     "O-methyl propanethioate"),
    # S-dithioate esters
    ("CC(=S)SC",      "S-methyl ethanedithioate"),
    ("CC(=S)SCC",     "S-ethyl ethanedithioate"),
    ("CCC(=S)SC",     "S-methyl propanedithioate"),
    # regression: S-thioester (C=O) unchanged
    ("CC(=O)SC",      "S-methyl ethanethioate"),
    ("CCC(=O)SC",     "S-methyl propanethioate"),
    # regression: thioic acids unchanged
    ("CC(=S)O",       "ethanethioic O-acid"),
    ("CC(=O)S",       "ethanethioic S-acid"),
    ("CC(=S)S",       "ethanedithioic acid"),
])
def test_phase251_thioester_variants(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
