"""Phase 547: Zwitterionic aliphatic azide [N-][N+]#N notation (IUPAC 2013 P-68.3.1).
R-[N-]-[N+]#N is equivalent to R-N=[N+]=[N-] and should give the same name.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("C[N-][N+]#N",   "azidomethane"),
    ("CC[N-][N+]#N",  "azidoethane"),
    ("CCC[N-][N+]#N", "azidopropane"),
    ("CCCC[N-][N+]#N", "azidobutane"),
])
def test_phase547_zwitterionic_azide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
