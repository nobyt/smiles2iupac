"""Phase 198: カルボニルジハライド命名 (IUPAC 2013 P-65.1.2.5)

  FC(F)=O   → carbonyl difluoride
  ClC(=O)Cl → carbonyl dichloride
  FC(Cl)=O  → carbonyl chloride fluoride

C=O に 2 個のハロゲンが付く場合は "carbonyl di{halide}" 形式。
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 対称カルボニルジハライド
    ("FC(F)=O",    "carbonyl difluoride"),
    ("ClC(=O)Cl",  "carbonyl dichloride"),
    ("BrC(=O)Br",  "carbonyl dibromide"),
    ("IC(=O)I",    "carbonyl diiodide"),
    # 非対称カルボニルハライド (アルファベット順)
    ("FC(Cl)=O",   "carbonyl chloride fluoride"),
    ("FC(Br)=O",   "carbonyl bromide fluoride"),
    ("ClC(=O)Br",  "carbonyl bromide chloride"),
    # 回帰: 通常のアシルハライドは変わらない
    ("CC(=O)Cl",   "acetyl chloride"),
    ("CC(=O)F",    "acetyl fluoride"),
    ("CCC(=O)Cl",  "propanoyl chloride"),
    ("ClCC(=O)Cl", "2-chloroethanoyl chloride"),
])
def test_phase198_carbonyl_dihalide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
