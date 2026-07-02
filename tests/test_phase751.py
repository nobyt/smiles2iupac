"""Phase 751: hydroxy takes precedence over amino as principal characteristic group
in N-heterocycles (IUPAC PCG hierarchy: -ol > -amine).

When both OH (at α-position) and NH2 are present, the OH generates the -one suffix
and the amino group becomes a substituent prefix — e.g. cytosine.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Cytosine: 4-aminopyrimidin-2-ol → 4-amino-1H-pyrimidin-2-one
    ("Nc1ccnc(O)n1",            "4-amino-1H-pyrimidin-2-one"),
    # Regression: amino-only (no OH present) unchanged
    ("Nc1ccncc1",               "pyridin-4-amine"),
    ("Nc1ccccn1",               "pyridin-2-amine"),
    # Regression: hydroxy-only still correct
    ("Oc1ccncn1",               "pyrimidin-4(3H)-one"),
    ("Oc1ncccn1",               "1H-pyrimidin-2-one"),
    ("Oc1ccccn1",               "1H-pyridin-2-one"),
    # Regression: beta-OH stays as -ol even with amino group
    ("Nc1cccnc1",               "pyridin-3-amine"),
    # Phase 101 regression: thiol tautomer updated
    ("Sc1ccncc1",               "pyridin-4(1H)-thione"),
])
def test_phase751_hydroxy_over_amino_pcg(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
