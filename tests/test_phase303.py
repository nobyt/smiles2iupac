"""Phase 303: N-substituted sulfamic acid naming (IUPAC 2013 P-65.3.1.2).

RHN-SO3H → N-Rsulfamic acid; R2N-SO3H → N,N-diRsulfamic acid.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("CNS(=O)(=O)O",    "N-methylsulfamic acid"),
    ("CCNS(=O)(=O)O",   "N-ethylsulfamic acid"),
    ("CN(C)S(=O)(=O)O", "N,N-dimethylsulfamic acid"),
    # regression: unsubstituted unchanged
    ("NS(=O)(=O)O",     "sulfamic acid"),
])
def test_phase303_n_substituted_sulfamic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
