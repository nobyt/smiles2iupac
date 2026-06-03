"""Phase 206: benzaldoxime and ring-attached aldoxime naming (IUPAC 2013 P-68.3.1)

  ON=Cc1ccccc1         → benzaldoxime  (retained: benzaldehyde oxime)
  ON=Cc1ccc(O)cc1      → 4-hydroxybenzaldoxime

Arene-attached aldoximes use benzaldoxime as retained parent name.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Benzaldoxime retained name
    ("ON=Cc1ccccc1",      "benzaldoxime"),
    # With ring substituents
    ("ON=Cc1ccc(O)cc1",   "4-hydroxybenzaldoxime"),
    ("ON=Cc1ccc(Cl)cc1",  "4-chlorobenzaldoxime"),
    # regression: non-aromatic aldoximes still work
    ("ON=CC",             "ethanal oxime"),
    ("ON=CCC",            "propanal oxime"),
    # regression: ketoxime unaffected
    ("ON=C(C)C",          "propan-2-one oxime"),
])
def test_phase206_benzaldoxime(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
