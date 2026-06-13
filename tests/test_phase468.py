"""Phase 468: thiazolo[4,5-c], oxazolo[4,5-c], isoxazolo-c pyridine isomers,
2H-[1,2,3]triazolo b/c pyridine, thiazolo-c-pyridazine isomers
(IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # thiazolo[4,5-c]pyridine (S at 1, N at 3; between=[N,C,S] → lower=4,upper=5 → [4,5-c])
    ("c1cc2scnc2cn1",    "thiazolo[4,5-c]pyridine"),
    # 2H-[1,2,3]triazolo[4,5-c]pyridine
    ("c1cc2n[nH]nc2cn1", "2H-[1,2,3]triazolo[4,5-c]pyridine"),
    # 1H-[1,2,3]triazolo[5,4-c]pyridine
    ("c1cc2nn[nH]c2cn1", "1H-[1,2,3]triazolo[5,4-c]pyridine"),
    # 2H-[1,2,3]triazolo[4,5-b]pyridine
    ("c1cnc2n[nH]nc2c1", "2H-[1,2,3]triazolo[4,5-b]pyridine"),
    # thiazolo-c-pyridazine isomers
    ("c1cc2scnc2nn1",    "thiazolo[4,5-c]pyridazine"),
    ("c1cc2ncsc2nn1",    "thiazolo[5,4-c]pyridazine"),
    # oxazolo[4,5-c]pyridine
    ("c1cc2ocnc2cn1",    "oxazolo[4,5-c]pyridine"),
    # isoxazolo-c-pyridine isomers
    ("c1cc2cnoc2cn1",    "isoxazolo[5,4-c]pyridine"),
    ("c1cc2oncc2cn1",    "isoxazolo[4,5-c]pyridine"),
    ("c1cc2nocc2cn1",    "isoxazolo[4,3-c]pyridine"),
])
def test_phase468_thiazolo_oxazolo_isoxazolo_c(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
