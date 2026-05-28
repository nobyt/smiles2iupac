"""Phase 125: アミン N-オキシド命名 (IUPAC 2013 P-68.1.1)"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # trimethylamine N-oxide (N,N-dimethylmethanamine N-oxide が PIN)
    ("C[N+](C)(C)[O-]", "N,N-dimethylmethanamine N-oxide"),
    # N-methyl-N-ethylamine N-oxide (ethane is principal chain)
    ("C[N+](CC)[O-]", "N-methylethanamine N-oxide"),
    # triethylamine N-oxide
    ("CC[N+](CC)(CC)[O-]", "N,N-diethylethanamine N-oxide"),
    # 回帰: 通常アミン (N-oxide なし)
    ("CN(C)C", "N,N-dimethylmethanamine"),
    ("CCN(C)C", "N,N-dimethylethanamine"),
    # 回帰: nitroso (N=O 二重結合、N-oxide と異なる → 影響なし)
    ("CN=O", "nitrosomethane"),
])
def test_phase125_amine_n_oxide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
