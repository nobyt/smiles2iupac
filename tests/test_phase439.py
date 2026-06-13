"""Phase 439: oxazolo/thiazolo-pyridine and -pyrazine bicyclic retained names
(IUPAC 2013 P-31.1.3).

Names corrected in Phase 470: [4,5-b]/[5,4-b] were swapped; pyrimidine entries
were actually pyrazine fusions (ring6 N's in 1,4 positions, not 1,3).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # oxazolo[5,4-b]pyridine — O adjacent to lower junction (C5 of oxazole at C2-pyr)
    ("c1cnc2ocnc2c1",              "oxazolo[5,4-b]pyridine"),
    # oxazolo[4,5-b]pyridine — N adjacent to lower junction (C4 of oxazole at C2-pyr)
    ("c1cnc2ncoc2c1",              "oxazolo[4,5-b]pyridine"),
    # thiazolo[5,4-b]pyridine — S adjacent to lower junction (C5 of thiazole at C2-pyr)
    ("c1cnc2scnc2c1",              "thiazolo[5,4-b]pyridine"),
    # thiazolo[4,5-b]pyridine — N adjacent to lower junction (C4 of thiazole at C2-pyr)
    ("c1cnc2ncsc2c1",              "thiazolo[4,5-b]pyridine"),
    # ring6 is pyrazine (N at 1,4 positions), not pyrimidine
    ("c1cnc2ocnc2n1",              "oxazolo[4,5-e]pyrazine"),
    ("c1cnc2scnc2n1",              "thiazolo[4,5-e]pyrazine"),
    # regression: 3H-imidazo[4,5-b]pyridine unchanged
    ("c1cnc2nc[nH]c2c1",           "3H-imidazo[4,5-b]pyridine"),
    # regression: furo[2,3-b]pyridine unchanged (Phase 436)
    ("c1cnc2occc2c1",              "furo[2,3-b]pyridine"),
])
def test_phase439_oxazolo_thiazolo(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
