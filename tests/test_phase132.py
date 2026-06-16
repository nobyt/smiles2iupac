"""Phase 132: ナフチリジン全異性体・追加縮合ヘテロ芳香族保留名 (IUPAC 2013 P-31.1.3)

naphthyridine isomers (all 6), indolizine, indazole, benzotriazole,
benzisoxazole, benzoxadiazole
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ナフチリジン全6異性体 (两N in naphthalene, correct canonical SMILES)
    ("c1cnc2cccnc2c1", "1,5-naphthyridine"),
    ("c1cnc2ccncc2c1", "1,6-naphthyridine"),
    ("c1cnc2cnccc2c1", "1,7-naphthyridine"),
    ("c1cnc2ncccc2c1", "1,8-naphthyridine"),
    ("c1cc2cnccc2cn1", "2,6-naphthyridine"),
    ("c1cc2ccncc2cn1", "2,7-naphthyridine"),
    # indolizine (pyrrolo[1,2-a]pyridine, N at ring junction)
    ("c1ccn2cccc2c1", "indolizine"),
    # 1H-indazole (benzo[c]pyrazole)
    ("c1ccc2[nH]ncc2c1", "1H-indazole"),
    # 1H-benzotriazole (benzo[d][1,2,3]triazole)
    ("c1ccc2[nH]nnc2c1", "1H-benzotriazole"),
    # 1,2-benzisoxazole (benzo[d]isoxazole)
    ("c1ccc2oncc2c1", "1,2-benzisoxazole"),
    # 2,1,3-benzoxadiazole (benzofurazan)
    ("c1ccc2nonc2c1", "2,1,3-benzoxadiazole"),
    # 回帰: Phase 131 tricyclics unchanged
    ("c1ccc2c(c1)[nH]c1ccccc12", "9H-carbazole"),
    ("c1ccc2c(c1)sc1ccccc12", "dibenzothiophene"),
    ("c1ccc2c(c1)oc1ccccc12", "dibenzofuran"),
    # 回帰: Phase 130 fused compounds unchanged
    ("c1ccc2ocnc2c1", "1,3-benzoxazole"),
    ("c1ccc2scnc2c1", "1,3-benzothiazole"),
    ("c1ccc2nc3ccccc3nc2c1", "phenazine"),
])
def test_phase132_naphthyridines_and_more(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
