"""Phase 495: 1H-pyrazolo and 1H-pyrrolo fused with pyridine at b- and c-bonds
(IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("c1cnc2c[nH]nc2c1",  "1H-pyrazolo[4,3-b]pyridine"),
    ("c1c[nH]c2cncc-2c1",  "1H-pyrrolo[3,4-b]pyridine"),
    ("c1cc2cn[nH]c2cn1",  "1H-pyrazolo[5,4-c]pyridine"),
    ("C1=NCc2ccncc21",  "1H-pyrrolo[3,4-c]pyridine"),
])
def test_phase495_pyrazolo_pyrrolo_pyridine_bc(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
