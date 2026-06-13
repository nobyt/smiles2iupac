"""Phase 470: fix Phase 439 b-bond naming errors (oxazolo/thiazolo [4,5] vs [5,4]
swapped; pyrimidine entries were pyrazine); add real oxazolo/thiazolo[x,y-d]pyrimidine,
isoxazolo[x,y-b]pyridine, and 1H-pyrazolo[3,4-b]pyridine
(IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # --- Phase 439 corrections ---
    # b-bond pyridine: [4,5] vs [5,4] was swapped; O-adjacent lower junc → [5,4]
    ("c1cnc2ocnc2c1",    "oxazolo[5,4-b]pyridine"),
    ("c1cnc2ncoc2c1",    "oxazolo[4,5-b]pyridine"),
    ("c1cnc2scnc2c1",    "thiazolo[5,4-b]pyridine"),
    ("c1cnc2ncsc2c1",    "thiazolo[4,5-b]pyridine"),
    # Phase 439 pyrimidine entries were actually pyrazine (e-bond, N at 1,4)
    ("c1cnc2ocnc2n1",    "oxazolo[4,5-e]pyrazine"),
    ("c1cnc2scnc2n1",    "thiazolo[4,5-e]pyrazine"),
    # --- New pyrimidine d-bond fusions (N at 1,3 in ring6) ---
    ("c1ncc2ncoc2n1",    "oxazolo[4,5-d]pyrimidine"),
    ("c1ncc2ocnc2n1",    "oxazolo[5,4-d]pyrimidine"),
    ("c1ncc2ncsc2n1",    "thiazolo[4,5-d]pyrimidine"),
    ("c1ncc2scnc2n1",    "thiazolo[5,4-d]pyrimidine"),
    # --- isoxazolo-b-pyridine (O=1,N=2 isoxazole) ---
    ("c1cnc2nocc2c1",    "isoxazolo[3,4-b]pyridine"),
    ("c1cnc2cnoc2c1",    "isoxazolo[4,5-b]pyridine"),
    ("c1cnc2oncc2c1",    "isoxazolo[5,4-b]pyridine"),
    # --- 1H-pyrazolo[3,4-b]pyridine ---
    ("c1cnc2n[nH]cc2c1", "1H-pyrazolo[3,4-b]pyridine"),
])
def test_phase470_b_bond_fixes_and_pyrimidine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
