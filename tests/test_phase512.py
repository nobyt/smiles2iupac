"""Phase 512: partially saturated 1H-pyrrole naming (2,3-dihydro-1H-pyrrole etc.)
(IUPAC 2013 P-31.1.6: dihydro/tetrahydro naming for 5-membered N-heterocycles).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 5-membered N ring, one double bond (pyrroline)
    ("C1=CCCN1",  "2,3-dihydro-1H-pyrrole"),
    ("C1=CCNC1",  "2,5-dihydro-1H-pyrrole"),
    # 5-membered O/S rings still correct
    ("C1=CCCO1",  "2,3-dihydrofuran"),
    ("C1=CCCS1",  "2,3-dihydrothiophene"),
    # 6-membered N rings unaffected
    ("C1=CCCCN1", "1,2,3,4-tetrahydropyridine"),
    ("C1CCCCN1",  "piperidine"),
])
def test_phase512_dihydropyrrole(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
