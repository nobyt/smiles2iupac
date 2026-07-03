"""Phase 787: 7H-purine, 1H-imidazo[4,5-e]pyrazine, thieno/furo[2,3-c]pyridazine α-ol/thiol → tautomers.

- 7H-purine C8 → 7H-purin-8(9H)-one/thione
- 1H-imidazo[4,5-e]pyrazine C2 → 2(3H)-one/thione; C5 → 5(4H)-one/thione
- thieno[2,3-c]pyridazine C3 → 3(2H)-one/thione
- furo[2,3-c]pyridazine C3 → 3(2H)-one/thione
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 7H-purine C8-OH/SH
    ("Oc1nc2ncncc2[nH]1",   "7H-purin-8(9H)-one"),
    ("Sc1nc2ncncc2[nH]1",   "7H-purin-8(9H)-thione"),
    # 1H-imidazo[4,5-e]pyrazine C2-OH/SH
    ("Oc1nc2nccnc2[nH]1",   "1H-imidazo[4,5-e]pyrazin-2(3H)-one"),
    ("Sc1nc2nccnc2[nH]1",   "1H-imidazo[4,5-e]pyrazin-2(3H)-thione"),
    # 1H-imidazo[4,5-e]pyrazine C5-OH/SH
    ("Oc1cnc2nc[nH]c2n1",   "1H-imidazo[4,5-e]pyrazin-5(4H)-one"),
    ("Sc1cnc2nc[nH]c2n1",   "1H-imidazo[4,5-e]pyrazin-5(4H)-thione"),
    # thieno[2,3-c]pyridazine C3-OH/SH
    ("Oc1cc2ccsc2nn1",       "thieno[2,3-c]pyridazin-3(2H)-one"),
    ("Sc1cc2ccsc2nn1",       "thieno[2,3-c]pyridazin-3(2H)-thione"),
    # furo[2,3-c]pyridazine C3-OH/SH
    ("Oc1cc2ccoc2nn1",       "furo[2,3-c]pyridazin-3(2H)-one"),
    ("Sc1cc2ccoc2nn1",       "furo[2,3-c]pyridazin-3(2H)-thione"),
    # Regressions: parent rings unchanged
    ("c1nc2ncncc2[nH]1",    "7H-purine"),
    ("c1nc2nccnc2[nH]1",    "1H-imidazo[4,5-e]pyrazine"),
    ("c1csc2nnccc12",        "thieno[2,3-c]pyridazine"),
    ("c1coc2nnccc12",        "furo[2,3-c]pyridazine"),
])
def test_phase787_purine_imidazo_pyridazine_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
