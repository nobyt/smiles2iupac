"""Phase 510: polyol tetraol/pentaol/hexaol suffix support
(IUPAC 2013 P-63.1.2: alkanetetraol, alkanepentaol, alkanehexaol naming).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("OCC(O)C(O)CO",         "butane-1,2,3,4-tetraol"),
    ("OCC(O)C(O)C(O)CO",     "pentane-1,2,3,4,5-pentaol"),
    ("OCC(O)C(O)C(O)C(O)CO", "hexane-1,2,3,4,5,6-hexaol"),
    # Lower polyols still correct
    ("OCCO",    "ethane-1,2-diol"),
    ("OCC(O)CO", "propane-1,2,3-triol"),
])
def test_phase510_higher_polyol(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
