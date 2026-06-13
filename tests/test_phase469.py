"""Phase 469: oxazolo/isoxazolo-c-pyridazine, 1H-imidazo[4,5-c]pyridazine,
1H/2H-[1,2,3]triazolo-c/b-pyridazine isomers
(IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # oxazolo-c-pyridazine
    ("c1cc2ncoc2nn1",    "oxazolo[5,4-c]pyridazine"),
    ("c1cc2ocnc2nn1",    "oxazolo[4,5-c]pyridazine"),
    # isoxazolo-c-pyridazine
    ("c1cc2cnoc2nn1",    "isoxazolo[5,4-c]pyridazine"),
    ("c1cc2oncc2nn1",    "isoxazolo[4,5-c]pyridazine"),
    ("c1cc2nocc2nn1",    "isoxazolo[4,3-c]pyridazine"),
    # imidazo-c-pyridazine
    ("c1cc2[nH]cnc2nn1", "1H-imidazo[4,5-c]pyridazine"),
    # [1,2,3]triazolo-c-pyridazine isomers
    ("c1cc2[nH]nnc2nn1", "1H-[1,2,3]triazolo[4,5-c]pyridazine"),
    ("c1cc2n[nH]nc2nn1", "2H-[1,2,3]triazolo[4,5-c]pyridazine"),
    ("c1cc2nn[nH]c2nn1", "1H-[1,2,3]triazolo[5,4-c]pyridazine"),
    # corrected Phase 474: ring6 is pyrazine (c1cnc2 template), not pyridazine
    ("c1cnc2[nH]nnc2n1", "1H-[1,2,3]triazolo[4,5-e]pyrazine"),
    ("c1cnc2n[nH]nc2n1", "2H-[1,2,3]triazolo[4,5-e]pyrazine"),
])
def test_phase469_oxazolo_isoxazolo_triazolo_pyridazine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
