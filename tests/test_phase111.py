"""Phase 111: ジアリールケトン括弧表記 (IUPAC P-14.5.7, P-66.6.1)"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 異なるアリール基 → 括弧必要な側のみ括弧、または両方括弧
    ("O=C(c1ccccc1)c1ccc(Cl)cc1", "(4-chlorophenyl)(phenyl)methanone"),
    ("O=C(c1ccccc1)c1ccncc1", "phenyl(pyridin-4-yl)methanone"),
    ("O=C(c1ccccc1)c1ccc(C)cc1", "(4-methylphenyl)(phenyl)methanone"),
    ("O=C(c1ccccc1)c1ccc(OC)cc1", "(4-methoxyphenyl)(phenyl)methanone"),
    # 両方同じ → diphenyl
    ("O=C(c1ccccc1)c1ccccc1", "diphenylmethanone"),
    # 回帰: アルキルフェニルケトン (ロカントあり)
    ("O=C(C)c1ccccc1", "acetophenone"),
    ("O=C(CC)c1ccccc1", "1-phenylpropan-1-one"),
])
def test_phase111_diaryl_ketone_parentheses(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
