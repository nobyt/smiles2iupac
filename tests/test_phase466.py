"""Phase 466: thieno/furo/pyrrolo-[x,y-c]pyridine isomers,
thieno/pyrrolo-[x,y-c]pyridazine, and fix thieno[3,2-c]pyridine
(IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # corrected from Phase 442 (was thieno[2,3-c]pyridine)
    ("c1cc2sccc2cn1",     "thieno[3,2-c]pyridine"),
    # new c-bond pyridine isomers
    ("c1cc2ccsc2cn1",     "thieno[2,3-c]pyridine"),
    ("c1cc2cscc2cn1",     "thieno[3,4-c]pyridine"),
    ("c1cc2ccoc2cn1",     "furo[2,3-c]pyridine"),
    ("c1cc2cocc2cn1",     "furo[3,4-c]pyridine"),
    ("c1cc2cc[nH]c2cn1",  "1H-pyrrolo[2,3-c]pyridine"),
    # c-bond pyridazine fusions
    ("c1cc2sccc2nn1",     "thieno[3,2-c]pyridazine"),
    ("c1cc2nccc-2[nH]n1",  "1H-pyrrolo[3,2-c]pyridazine"),
])
def test_phase466_thienofuro_c_pyridine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
