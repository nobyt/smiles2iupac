"""Phase 473: fix Phase 435/436/437 pyrazine-mislabeled-as-pyrimidine entries;
add remaining pyrazine e-bond fusions (IUPAC 2013 P-31.1.3).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # --- Phase 435/436/437 corrections (ring6 is pyrazine, not pyrimidine) ---
    ("c1c[nH]c2ccnc-2n1",   "1H-pyrrolo[2,3-e]pyrazine"),
    ("c1cnc2[nH]ncc2n1",   "1H-pyrazolo[4,5-e]pyrazine"),
    ("c1cnc2sccc2n1",      "thieno[2,3-e]pyrazine"),
    ("c1cnc2occc2n1",      "furo[2,3-e]pyrazine"),
    ("c1cnc2cocc2n1",      "furo[3,4-e]pyrazine"),
    ("c1cnc2[nH]cnc2n1",   "1H-imidazo[4,5-e]pyrazine"),
    ("c1cnc2cscc2n1",      "thieno[3,4-e]pyrazine"),
    # --- New pyrazine e-bond fusions ---
    # (pyrazine C2 symmetry: [2,3]/[3,2] same compound; [4,5]/[5,4] same for oxazolo/thiazolo)
    ("c1c[nH]c2cncc-2n1",   "1H-pyrrolo[3,4-e]pyrazine"),
    # imidazo[4,5-e] and [5,4-e] are the SAME compound (pyrazine C2); preferred name [4,5]
    # --- Regressions: already-correct pyrazine entries (Phase 439/470) ---
    ("c1cnc2ocnc2n1",      "oxazolo[4,5-e]pyrazine"),
    ("c1cnc2scnc2n1",      "thiazolo[4,5-e]pyrazine"),
])
def test_phase473_pyrazine_e_bond_fusions(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
