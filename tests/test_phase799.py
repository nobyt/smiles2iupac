"""Phase 799: imidazo[1,2-a]pyrimidine α-ol/thiol → tautomers.

- C3 → 3(2H)-one/thione (adjacent to N2)
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # imidazo[1,2-a]pyrimidine C3-OH/SH
    ("Oc1cnc2ncccn12",   "imidazo[1,2-a]pyrimidin-3(2H)-one"),
    ("Sc1cnc2ncccn12",   "imidazo[1,2-a]pyrimidin-3(2H)-thione"),
    # Regression: parent ring unchanged
    ("c1cnc2ncccn12",    "imidazo[1,2-a]pyrimidine"),
])
def test_phase799_imidazo_pyrimidine_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
