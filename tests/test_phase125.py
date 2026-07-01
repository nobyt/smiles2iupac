"""Phase 125: アミン N-オキシド命名 (IUPAC 2013 P-62.4.1: suffix is 'oxide' not 'N-oxide')"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # trimethylamine oxide (PIN: N,N-dimethylmethanamine oxide)
    ("C[N+](C)(C)[O-]", "N,N-dimethylmethanamine oxide"),
    # N-methyl-N-ethylamine oxide (ethane is principal chain)
    ("C[NH+](CC)[O-]", "N-methylethanamine oxide"),
    # triethylamine oxide
    ("CC[N+](CC)(CC)[O-]", "N,N-diethylethanamine oxide"),
    # 回帰: 通常アミン (N-oxide なし)
    ("CN(C)C", "N,N-dimethylmethanamine"),
    ("CCN(C)C", "N,N-dimethylethanamine"),
    # 回帰: nitroso (N=O 二重結合、N-oxide と異なる → 影響なし)
    ("CN=O", "nitrosomethane"),
])
def test_phase125_amine_n_oxide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
