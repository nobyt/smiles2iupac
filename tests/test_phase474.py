"""Phase 474: isoxazolo/pyrazolo/[1,2,3]triazolo[x,y-e]pyrazine bicyclic retained names
(IUPAC 2013 P-31.1.3 fusion nomenclature).

Pyrazine C2 symmetry: [x,y] and [y,x] are the same compound; prefer lower locants.
[1,2,3]triazolo entries also fix Phase 469 wrong "b]pyridazine" labels.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # isoxazolo-e-pyrazine (2 unique compounds due to pyrazine C2 symmetry)
    ("c1cnc2nocc2n1",    "isoxazolo[3,4-e]pyrazine"),
    ("c1cnc2oncc2n1",    "isoxazolo[4,5-e]pyrazine"),
    # 1H-pyrazolo-e-pyrazine ([4,5] and [3,4] unique; C2 makes [5,4]=[4,5], [4,3]=[3,4])
    ("c1cnc2[nH]ncc2n1", "1H-pyrazolo[4,5-e]pyrazine"),
    ("c1cnc2n[nH]cc2n1", "1H-pyrazolo[3,4-e]pyrazine"),
    # [1,2,3]triazolo-e-pyrazine (corrected from wrong "b]pyridazine" in Phase 469)
    ("c1cnc2[nH]nnc2n1", "1H-[1,2,3]triazolo[4,5-e]pyrazine"),
    ("c1cnc2n[nH]nc2n1", "2H-[1,2,3]triazolo[4,5-e]pyrazine"),
])
def test_phase474_pyrazine_e_bond_5membered(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
