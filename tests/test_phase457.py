"""Phase 457: pyrido[2,3-g]quinazoline, pyrimido[5,4-b/g] series,
pyrido[2,3-g][1,5]naphthyridine (IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # pyrido-quinazoline: N of pyridine adjacent to C6 (lower end of bond g)
    ("c1cnc2cc3cncnc3cc2c1",      "pyrido[2,3-g]quinazoline"),
    # pyrimido[5,4-g][1,5]naphthyridine: pyrimidine N3 at C7 (higher) of [1,5]naphthyridine
    ("c1cnc2cc3ncncc3nc2c1",      "pyrimido[5,4-g][1,5]naphthyridine"),
    # pyrimido[5,4-b]quinoline: pyrimidine N3 at C3 (higher) of quinoline bond b
    ("c1ccc2nc3cncnc3cc2c1",      "pyrimido[5,4-b]quinoline"),
    # pyrimido[5,4-g][1,8]naphthyridine: pyrimidine N3 at C7 (higher) of [1,8]naphthyridine
    ("c1cnc2nc3ncncc3cc2c1",      "pyrimido[5,4-g][1,8]naphthyridine"),
    # pyrido[2,3-g][1,5]naphthyridine: [1,5]naphthyridine base, pyridine N at C6 (lower)
    ("c1cnc2nc3cccnc3cc2c1",      "pyrido[2,3-g][1,5]naphthyridine"),
])
def test_phase457_pyrimido_fused(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
