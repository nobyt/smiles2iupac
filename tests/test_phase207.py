"""Phase 207: cyclic ketone oxime naming (IUPAC 2013 P-68.3.1)

  ON=C1CCCC1    → cyclopentanone oxime
  ON=C1CCCCC1   → cyclohexanone oxime

Ring ketones (cycloalkanones) form oximes via C=N-OH on the ring carbon.
The name is the cycloalkanone name + ' oxime' (retained method).
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Simple cyclic ketone oximes
    ("ON=C1CCCC1",   "cyclopentanone oxime"),
    ("ON=C1CCCCC1",  "cyclohexanone oxime"),
    ("ON=C1CCC1",    "cyclobutanone oxime"),
    ("ON=C1CC1",     "cyclopropanone oxime"),
    # regression: acyclic ketoxime unaffected
    ("ON=C(C)C",     "propan-2-one oxime"),
    ("ON=C(CC)CC",   "pentan-3-one oxime"),
    # regression: acyclic aldoxime unaffected
    ("ON=CC",        "ethanal oxime"),
])
def test_phase207_cyclic_ketoxime(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
