"""Phase 781: pyrazolo[1,5-a]pyrimidine and 1H-imidazo[4,5-b]pyridine α-ol/thiol → tautomers.

pyrazolo[1,5-a]pyrimidine alpha positions:
- C2 (alpha to N1) → 2(1H)-one/thione
- C5 (alpha to N4 junction) → 5(4H)-one/thione
- C7 (alpha to N4 junction) → 7(4H)-one/thione

1H-imidazo[4,5-b]pyridine alpha positions:
- C2 (alpha to N3, with permanent N1-H) → 2(3H)-one/thione
- C5 (alpha to N4 junction) → 5(4H)-one/thione
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # pyrazolo[1,5-a]pyrimidine C2 (alpha to N1)
    ("Oc1cc2ncccn2n1",       "pyrazolo[1,5-a]pyrimidin-2(1H)-one"),
    ("Sc1cc2ncccn2n1",       "pyrazolo[1,5-a]pyrimidin-2(1H)-thione"),
    # pyrazolo[1,5-a]pyrimidine C5 (alpha to junction N4)
    ("Oc1ccn2nccc2n1",       "pyrazolo[1,5-a]pyrimidin-5(4H)-one"),
    ("Sc1ccn2nccc2n1",       "pyrazolo[1,5-a]pyrimidin-5(4H)-thione"),
    # pyrazolo[1,5-a]pyrimidine C7 (alpha to junction N4)
    ("Oc1ccnc2ccnn12",       "pyrazolo[1,5-a]pyrimidin-7(4H)-one"),
    ("Sc1ccnc2ccnn12",       "pyrazolo[1,5-a]pyrimidin-7(4H)-thione"),
    # 1H-imidazo[4,5-b]pyridine C2 (alpha to N3, H on N3)
    ("Oc1nc2cccnc2[nH]1",   "1H-imidazo[4,5-b]pyridin-2(3H)-one"),
    ("Sc1nc2cccnc2[nH]1",   "1H-imidazo[4,5-b]pyridin-2(3H)-thione"),
    # 1H-imidazo[4,5-b]pyridine C5 (alpha to junction N4)
    ("Oc1ccc2nc[nH]c2n1",   "1H-imidazo[4,5-b]pyridin-5(4H)-one"),
    ("Sc1ccc2nc[nH]c2n1",   "1H-imidazo[4,5-b]pyridin-5(4H)-thione"),
    # Regressions: parent rings unchanged
    ("c1cnc2ccnn2c1",        "pyrazolo[1,5-a]pyrimidine"),
    ("c1cnc2[nH]cnc2c1",    "1H-imidazo[4,5-b]pyridine"),
])
def test_phase781_pyrazolo_imidazo_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
