"""Phase 105: ハロゲン・ヘテロ原子置換フェニル基の命名 (IUPAC P-31.1.3.4)"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ハロゲン置換フェニル on ヘテロアリール
    ("c1ccc(Cl)c(-c2ccncc2)c1", "4-(2-chlorophenyl)pyridine"),
    ("Clc1ccc(-c2ccncc2)cc1", "4-(4-chlorophenyl)pyridine"),
    ("Fc1ccc(-c2ccncc2)cc1", "4-(4-fluorophenyl)pyridine"),
    ("Brc1ccc(-c2ccncc2)cc1", "4-(4-bromophenyl)pyridine"),
    # OH・NH₂置換フェニル on ヘテロアリール
    ("Oc1ccc(-c2ccncc2)cc1", "4-(4-hydroxyphenyl)pyridine"),
    ("Nc1ccc(-c2ccncc2)cc1", "4-(4-aminophenyl)pyridine"),
    # ハロゲン on フェニル (メタ位 → 3-chlorophenyl)
    ("c1ccc(Cl)cc1-c1ccncc1", "4-(3-chlorophenyl)pyridine"),
    # 二置換フェニル on ヘテロアリール
    ("Clc1cc(-c2ccncc2)ccc1F", "4-(3-chloro-4-fluorophenyl)pyridine"),
    # 回帰: 無置換フェニルは "phenyl" のまま
    ("c1ccc(-c2ccncc2)cc1", "4-phenylpyridine"),
    # 回帰: 炭素置換フェニルは従来通り
    ("Cc1ccc(-c2ccncc2)cc1", "4-(4-methylphenyl)pyridine"),
    # 回帰: 単純ハロゲン置換ベンゼン
    ("Clc1ccccc1", "chlorobenzene"),
    ("Fc1ccccc1", "fluorobenzene"),
])
def test_phase105_halogenated_phenyl_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
