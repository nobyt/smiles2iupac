"""Phase 148: 縮合ヘテロ芳香族への外環官能基 (carboxylic acid / carbaldehyde 等)

fused heteroaromatic retained-name rings (indole, quinoline, benzofuran, etc.)
with exocyclic principal functional groups named as suffix rather than prefix.
  e.g. OC(=O)c1cc2ccccc2[nH]1 → 1H-indole-2-carboxylic acid
  (not "2-carboxyindole" which was the old incorrect output)
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # カルボン酸
    ("OC(=O)c1cc2ccccc2[nH]1",  "1H-indole-2-carboxylic acid"),
    ("OC(=O)c1ccc2ccccc2n1",    "quinoline-2-carboxylic acid"),
    ("OC(=O)c1cc2ccccc2o1",     "benzofuran-2-carboxylic acid"),
    ("OC(=O)c1cc2ccccc2s1",     "benzo[b]thiophene-2-carboxylic acid"),
    # カルボアルデヒド
    ("O=Cc1cc2ccccc2[nH]1",     "1H-indole-2-carbaldehyde"),
    ("O=Cc1ccc2ccccc2n1",       "quinoline-2-carbaldehyde"),
    # カルボキサミド
    ("NC(=O)c1ccc2ccccc2n1",    "quinoline-2-carboxamide"),
    # カルボニトリル
    ("N#Cc1ccc2ccccc2n1",       "quinoline-2-carbonitrile"),
    # 回帰: ナフタレン外環カルボン酸は既存パスで動作継続
    ("OC(=O)c1cccc2ccccc12",    "naphthalene-1-carboxylic acid"),
    ("OC(=O)c1ccc2ccccc2c1",    "naphthalene-2-carboxylic acid"),
    # 回帰: 単環ヘテロ芳香族カルボン酸
    ("OC(=O)c1cccnc1",          "pyridine-3-carboxylic acid"),
    ("c1ccoc1C(=O)O",           "furan-2-carboxylic acid"),
    ("c1cc[nH]c1C(=O)O",        "1H-pyrrole-2-carboxylic acid"),
])
def test_phase148_fused_hetero_exocyclic_fg(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
