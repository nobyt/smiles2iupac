"""Phase 871: phosphonimidic acids (R-P(=NH)(OH)2).

The imido (=NH replacing =O) analogue of phosphonic acid was mis-named:
CP(=N)(O)O -> "methylphosphonous acid" (dropping the =NH and demoting P(V) to
P(III)).

Naming (IUPAC 2013 P-67.1.4), parallel to phosphonamidic/phosphonohalidic: the
C-on-P substituent takes a P- locant, imido-N substituents take N- locants:

  CP(=N)(O)O   -> P-methylphosphonimidic acid
  CP(=NC)(O)O  -> N,P-dimethylphosphonimidic acid

The carbon-free parent HN=P(OH)3 is out of scope (namer needs a carbon).
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("CP(=N)(O)O",    "P-methylphosphonimidic acid"),
    ("CCP(=N)(O)O",   "P-ethylphosphonimidic acid"),
    ("CP(=NC)(O)O",   "N,P-dimethylphosphonimidic acid"),
    ("CP(=NCC)(O)O",  "N-ethyl-P-methylphosphonimidic acid"),
    # regression: neighbouring phosphorus acid classes unchanged
    ("CP(=O)(O)O",    "methylphosphonic acid"),
    ("CP(O)O",        "methylphosphonous acid"),
    ("CP(=S)(O)O",    "methylphosphonothioic acid"),
    ("CP(=O)(N)O",    "P-methylphosphonamidic acid"),
    ("CP(=O)(Cl)O",   "P-methylphosphonochloridic acid"),
    ("CNP(=O)(O)O",   "N-methylphosphoramidic acid"),
    ("CP(C)(=O)O",    "dimethylphosphinic acid"),
])
def test_phase871_phosphonimidic_acids(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
