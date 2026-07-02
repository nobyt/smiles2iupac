"""Phase 753: 1,3,5-triazin-2-ol → lactam; pyrimidine-2,4-dithiol → dithione (IUPAC 2013).

Two additions:
- 2-hydroxy-1,3,5-triazine → 1,3,5-triazin-2(1H)-one (α-position, Phase 744 pattern)
- pyrimidine-2,4-dithiol → pyrimidine-2,4(1H,3H)-dithione (parallel to Phase 752 diol)
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1,3,5-triazin-2-ol → lactam
    ("Oc1ncncn1",                  "1,3,5-triazin-2(1H)-one"),
    # pyrimidine-2,4-dithiol → dithione
    ("Sc1ccnc(S)n1",               "pyrimidine-2,4(1H,3H)-dithione"),
    # Regression: parent rings unaffected
    ("c1ncncn1",                   "1,3,5-triazine"),
    ("c1ccncn1",                   "pyrimidine"),
    # Regression: mono cases unchanged
    ("Oc1ncccn1",                  "1H-pyrimidin-2-one"),
    ("Sc1ccccn1",                  "pyridin-2(1H)-thione"),
    # Regression: diol unchanged
    ("Oc1ccnc(O)n1",               "pyrimidine-2,4(1H,3H)-dione"),
])
def test_phase753_triazinone_and_dithione(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
