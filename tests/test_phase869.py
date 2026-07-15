"""Phase 869: phosphono acid halides (R-P(=O)(X)(OH)).

The halide analogue of the Phase 868 phosphonamidic acids (OH replaced by a
halogen instead of NH2). R-P(=O)(Cl)(OH) was confidently mis-named
"methylphosphinic acid" -- dropping the halogen and mis-classifying the acid.

Naming (IUPAC 2013 P-67.1.4), parallel to phosphonamidic: the C-on-P
substituent takes a P- locant, halogen -> chlorid/bromid/fluorid/iodid:

  CP(=O)(Cl)O -> P-methylphosphonochloridic acid
  CP(=O)(F)O  -> P-methylphosphonofluoridic acid

Only members with a single OH remaining are acids; the dichloride
CP(=O)(Cl)Cl (no OH) is a different class and stays out of scope.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("CP(=O)(Cl)O",   "P-methylphosphonochloridic acid"),
    ("CP(=O)(Br)O",   "P-methylphosphonobromidic acid"),
    ("CP(=O)(F)O",    "P-methylphosphonofluoridic acid"),
    ("CP(=O)(I)O",    "P-methylphosphonoiodidic acid"),
    ("CCP(=O)(Cl)O",  "P-ethylphosphonochloridic acid"),
    # regression: neighbouring phosphorus classes unchanged
    ("CP(=O)(O)O",    "methylphosphonic acid"),
    ("CP(C)(=O)O",    "dimethylphosphinic acid"),
    ("CP(=O)(N)O",    "P-methylphosphonamidic acid"),
    ("CNP(=O)(O)O",   "N-methylphosphoramidic acid"),
    ("CP(=O)(Cl)OC",  "methyl methylphosphonate"),
])
def test_phase869_phosphono_halide_acids(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
