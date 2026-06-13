"""Phase 480: [1,2,4]triazolo[1,5-a]pyridine and related missing fused names
(IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # [1,2,4]triazolo[1,5-a]pyridine: N1-C5 junction, a-bond of pyridine
    ("c1ccn2ncnc2c1",  "[1,2,4]triazolo[1,5-a]pyridine"),
])
def test_phase480_triazolo_a_pyridine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
