"""Phase 142: 追加ヘテロ芳香族 保留名 (IUPAC 2013)

selenophene, imidazo[1,2-a]pyridine, thieno[2,3-b]pyridine,
pyrazolo[1,5-a]pyrimidine, 1H-pyrrolo[2,3-b]pyridine (7-azaindole),
1H-pyrrolo[3,2-b]pyridine (4-azaindole), thieno[3,2-b]thiophene
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # selenophene
    ("c1cc[se]c1",         "selenophene"),
    # imidazo-fused
    ("c1ccn2ccnc2c1",      "imidazo[1,2-a]pyridine"),
    ("c1cnc2ccnn2c1",      "pyrazolo[1,5-a]pyrimidine"),
    # thieno-fused (corrected names per OPSIN)
    ("c1cnc2ccsc2c1",      "thieno[3,2-b]pyridine"),
    ("c1cc2ccsc2s1",       "thieno[2,3-b]thiophene"),
    # pyrazolo-fused
    ("c1cnc2cn[nH]c2c1",   "1H-pyrazolo[4,5-b]pyridine"),
    # pyrrolo-fused (azaindoles)
    ("c1cnc2[nH]ccc2c1",   "1H-pyrrolo[2,3-b]pyridine"),
    ("c1cnc2cc[nH]c2c1",   "1H-pyrrolo[3,2-b]pyridine"),
    # 回帰: previously working heteroaromatics unchanged
    ("c1ccsc1",   "thiophene"),
    ("c1ccoc1",   "furan"),
    ("c1ccc2ncccc2c1", "quinoline"),
    ("c1ccc2[nH]ccc2c1", "1H-indole"),
    ("c1ccc2[nH]cnc2c1", "1H-benzimidazole"),
    ("c1ccc2ocnc2c1",    "1,3-benzoxazole"),
    ("c1ccc2scnc2c1",    "1,3-benzothiazole"),
])
def test_phase142_heteroaromatics(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
