"""Phase 130: 追加縮合ヘテロ芳香族保留名 (IUPAC 2013 P-31.1.3)

benzoxazole, benzothiazole, quinazoline, quinoxaline,
cinnoline, phthalazine, phenazine, pteridine
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1,3-benzoxazole (benzene + oxazole)
    ("c1ccc2ocnc2c1", "1,3-benzoxazole"),
    # 1,3-benzothiazole (benzene + thiazole)
    ("c1ccc2scnc2c1", "1,3-benzothiazole"),
    # quinazoline (benzo[d]pyrimidine, N at 1,3)
    ("c1ccc2ncncc2c1", "quinazoline"),
    # quinoxaline (benzo[g]pyrazine, N at 1,4)
    ("c1ccc2nccnc2c1", "quinoxaline"),
    # cinnoline (benzo[c]pyridazine, N at 1,2)
    ("c1ccc2nnccc2c1", "cinnoline"),
    # phthalazine (benzo[d]pyridazine)
    ("c1ccc2cnncc2c1", "phthalazine"),
    # phenazine (dibenzo[b,e]pyrazine)
    ("c1ccc2nc3ccccc3nc2c1", "phenazine"),
    # pteridine (pyrimido[4,5-d]pyrimidine)
    ("c1cnc2nccnc2n1", "pteridine"),
    # 回帰: previously added fused heterocycles unchanged
    ("c1ccc2[nH]cnc2c1", "1H-benzimidazole"),
    ("c1ccc2ncccc2c1", "quinoline"),
    ("c1ccc2occc2c1", "benzofuran"),
    ("c1ncc2[nH]cnc2n1", "7H-purine"),
])
def test_phase130_fused_heteroaromatics(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
