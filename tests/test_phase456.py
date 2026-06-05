"""Phase 456: pyrido[4,3-g]quinoline, pyrido[3,2-g]quinazoline, pyrimido[4,5-b]quinoline,
plus regressions for Phase 445/453/454 corrections (IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # new linear pyridine-benzene-pyridine isomers
    ("c1cnc2cc3cnccc3cc2c1",      "pyrido[4,3-g]quinoline"),
    # pyridine fused to quinazoline at bond g (N-adjacent end at C7)
    ("c1cnc2cc3ncncc3cc2c1",      "pyrido[3,2-g]quinazoline"),
    # pyrimidine fused to quinoline at bond b
    ("c1ccc2nc3ncncc3cc2c1",      "pyrimido[4,5-b]quinoline"),
    # regressions: Phase 445 fix
    ("c1cnc2cc3ccncc3cc2c1",      "pyrido[3,4-g]quinoline"),
    # regressions: Phase 453 fix
    ("c1cnc2cc3ncccc3cc2c1",      "pyrido[3,2-g]quinoline"),
    # regressions: Phase 454 fix
    ("c1cnc2cc3cccnc3cc2c1",      "pyrido[2,3-g]quinoline"),
])
def test_phase456_pyridoquinolines(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
