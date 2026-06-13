"""Phase 472: isoxazolo/pyrazolo/[1,2,3]triazolo[x,y-d]pyrimidine
(IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # isoxazolo-d-pyrimidine (O=1, N=2 isoxazole)
    ("c1ncc2nocc2n1",    "isoxazolo[4,3-d]pyrimidine"),
    ("c1ncc2oncc2n1",    "isoxazolo[4,5-d]pyrimidine"),
    ("c1ncc2conc2n1",    "isoxazolo[3,4-d]pyrimidine"),
    ("c1ncc2cnoc2n1",    "isoxazolo[5,4-d]pyrimidine"),
    # 1H-pyrazolo-d-pyrimidine (N=1,N=2 pyrazole, NH at 1)
    ("c1ncc2[nH]ncc2n1", "1H-pyrazolo[4,5-d]pyrimidine"),
    ("c1ncc2n[nH]cc2n1", "1H-pyrazolo[4,3-d]pyrimidine"),
    ("c1ncc2cn[nH]c2n1", "1H-pyrazolo[5,4-d]pyrimidine"),
    ("c1ncc2c[nH]nc2n1", "1H-pyrazolo[3,4-d]pyrimidine"),
    # [1,2,3]triazolo-d-pyrimidine
    ("c1ncc2[nH]nnc2n1", "1H-[1,2,3]triazolo[4,5-d]pyrimidine"),
    ("c1ncc2n[nH]nc2n1", "2H-[1,2,3]triazolo[4,5-d]pyrimidine"),
    ("c1ncc2nn[nH]c2n1", "1H-[1,2,3]triazolo[5,4-d]pyrimidine"),
])
def test_phase472_d_pyrimidine_fusions(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
