"""Phase 508: partially saturated 2-heteroatom 5-membered rings
(IUPAC 2013 P-31.1.6 dihydro naming for thiazoline, oxazoline, etc.)
and pyrazolidine retained name fix.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Partially saturated (4,5-dihydro) forms
    ("C1=NCCS1",  "4,5-dihydro-1,3-thiazole"),
    ("C1=NCCO1",  "4,5-dihydro-1,3-oxazole"),
    ("C1=NSCC1",  "4,5-dihydroisothiazole"),
    ("C1=NOCC1",  "4,5-dihydroisoxazole"),
    ("C1=NCCN1",  "4,5-dihydro-1H-imidazole"),
    ("C1=NNCC1",  "4,5-dihydro-1H-pyrazole"),
    # Fully saturated forms still correct
    ("C1CSCN1",   "thiazolidine"),
    ("C1COCN1",   "oxazolidine"),
    ("C1CNSC1",   "isothiazolidine"),
    ("C1CNOC1",   "isoxazolidine"),
    ("C1CNCN1",   "imidazolidine"),
    ("C1CCNN1",   "pyrazolidine"),
])
def test_phase508_partial_unsat_dihetero(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
