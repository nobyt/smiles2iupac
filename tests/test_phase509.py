"""Phase 509: partially saturated 2-heteroatom 6-membered rings
(IUPAC 2013 P-31.1.6 dihydro/tetrahydro naming for pyrimidine, pyridazine,
1,3-oxazine, 1,3-thiazine, and related partially saturated heterocycles).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Pyrimidine (N at 1,3) partially saturated forms
    ("C1=NCCCN1",  "1,4,5,6-tetrahydropyrimidine"),   # 1 C=N double bond
    ("C1=CCNC=N1", "1,6-dihydropyrimidine"),            # 2 double bonds
    # Pyridazine (N at 1,2) partially saturated
    ("C1=NNCCC1",  "1,4,5,6-tetrahydropyridazine"),
    # N+O 1,3-oxazine partially saturated
    ("C1=NCCCO1",  "5,6-dihydro-4H-1,3-oxazine"),
    # N+S 1,3-thiazine partially saturated
    ("C1=NCCCS1",  "5,6-dihydro-4H-1,3-thiazine"),
    # Fully saturated forms still correct
    ("C1NCCCN1",   "hexahydropyrimidine"),
    ("C1NCCCO1",   "1,3-oxazinane"),
    ("C1NCCCS1",   "1,3-thiazinane"),
    ("C1COCCN1",   "morpholine"),
])
def test_phase509_partial_unsat_6mem_dihetero(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
