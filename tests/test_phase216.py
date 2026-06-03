"""Phase 216: amino-substituted alkyl chains as ring substituents (IUPAC 2013)

  NCCc1ccccc1             → 2-phenylethylamine  (already works as main chain)
  NCCc1ccc(O)c(O)c1       → 4-(2-aminoethyl)benzene-1,2-diol  (was wrong)
  NCCc1ccc(O)cc1           → 4-(2-aminoethyl)phenol

N-containing substituents on alkyl chains need amino prefix naming.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Aminoethyl substituents on benzene rings
    ("NCCc1ccc(O)c(O)c1",  "4-(2-aminoethyl)benzene-1,2-diol"),
    ("NCCc1ccc(O)cc1",     "4-(2-aminoethyl)phenol"),
    ("NCc1ccccc1",         "phenylmethanamine"),  # benzylamine - main chain case
    # regression: straight chain amine still works
    ("NCCc1ccccc1",        "2-phenylethanamine"),
    # regression: phenol still works
    ("Oc1ccccc1",          "phenol"),
])
def test_phase216_amino_substituted_chain(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
