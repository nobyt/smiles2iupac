"""Phase 499: isothiazolo/isoxazolo/oxazolo/thiazolo fused with pyridazine at d-bond
(IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("c1nncc2sncc12",  "isothiazolo[5,4-d]pyridazine"),
    ("c1nncc2nscc12",  "isothiazolo[3,4-d]pyridazine"),
    ("c1nncc2oncc12",  "isoxazolo[5,4-d]pyridazine"),
    ("c1nncc2nocc12",  "isoxazolo[3,4-d]pyridazine"),
    ("c1nc2cnncc2o1",  "oxazolo[4,5-d]pyridazine"),
    ("c1nc2cnncc2s1",  "thiazolo[4,5-d]pyridazine"),
])
def test_phase499_d_bond_pyridazine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
