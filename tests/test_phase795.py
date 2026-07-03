"""Phase 795: 1H-[1,2,3]triazolo[4,5-e]pyrazine, pyrido[2,3-e]pyrimidine α-ol/thiol → tautomers.

- 1H-[1,2,3]triazolo[4,5-e]pyrazine C5 → 5(4H)-one/thione
- pyrido[2,3-e]pyrimidine C2 → 2(1H)-one/thione; C4 → 4(3H)-one/thione; C6 → 6(5H)-one/thione
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1H-[1,2,3]triazolo[4,5-e]pyrazine C5-OH/SH
    ("Oc1cnc2[nH]nnc2n1",  "1H-[1,2,3]triazolo[4,5-e]pyrazin-5(4H)-one"),
    ("Sc1cnc2[nH]nnc2n1",  "1H-[1,2,3]triazolo[4,5-e]pyrazin-5(4H)-thione"),
    # pyrido[2,3-e]pyrimidine C2-OH/SH
    ("Oc1ncc2ncccc2n1",    "pyrido[2,3-e]pyrimidin-2(1H)-one"),
    ("Sc1ncc2ncccc2n1",    "pyrido[2,3-e]pyrimidin-2(1H)-thione"),
    # pyrido[2,3-e]pyrimidine C4-OH/SH
    ("Oc1ncnc2cccnc12",    "pyrido[2,3-e]pyrimidin-4(3H)-one"),
    ("Sc1ncnc2cccnc12",    "pyrido[2,3-e]pyrimidin-4(3H)-thione"),
    # pyrido[2,3-e]pyrimidine C6-OH/SH
    ("Oc1ccc2ncncc2n1",    "pyrido[2,3-e]pyrimidin-6(5H)-one"),
    ("Sc1ccc2ncncc2n1",    "pyrido[2,3-e]pyrimidin-6(5H)-thione"),
    # Regressions: parent rings unchanged
    ("c1cnc2[nH]nnc2n1",   "1H-[1,2,3]triazolo[4,5-e]pyrazine"),
    ("c1ccc2ncncc2n1",     "pyrido[2,3-e]pyrimidine"),
])
def test_phase795_triazolo_pyrido_pyrimidine_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
