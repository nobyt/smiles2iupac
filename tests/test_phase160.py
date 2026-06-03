"""Phase 160: アリール基を末端に持つアルキル置換基の正確な命名

鎖状炭素が芳香族環に接続する場合の命名修正:
  c1ccc(CCc2ccccc2)cc1 → (2-phenylethyl)benzene  (旧: octylbenzene)
  c1ccc(Cc2ccccc2)cc1  → (phenylmethyl)benzene   (旧: heptylbenzene)
  c1ccc(C=Cc2ccccc2)cc1 → (2-phenylethenyl)benzene (旧: octylbenzene)

また nitrosobenzene の修正 (芳香環の N=O → nitrosobenzene) も含む。
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # アリール末端アルキル鎖 (Phase 160 主要修正)
    ("c1ccc(CCc2ccccc2)cc1",         "(2-phenylethyl)benzene"),
    ("c1ccc(Cc2ccccc2)cc1",          "(phenylmethyl)benzene"),
    ("c1ccc(C=Cc2ccccc2)cc1",        "(2-phenylethenyl)benzene"),
    # ニトロソベンゼン修正
    ("c1ccc(N=O)cc1",                "nitrosobenzene"),
    # 回帰: 通常の置換基は変わらない
    ("c1ccc(CCCC)cc1",               "butylbenzene"),
    ("c1ccc(C(C)C)cc1",              "(propan-2-yl)benzene"),
    ("OCc1ccccc1",                   "phenylmethanol"),
    ("OC(=O)Cc1ccccc1",              "phenylacetic acid"),
    ("c1ccc([N+](=O)[O-])cc1",       "nitrobenzene"),
    ("Cc1ccc([N+](=O)[O-])cc1",      "1-methyl-4-nitrobenzene"),
])
def test_phase160_aryl_terminated_chain(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
