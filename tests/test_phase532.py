"""Phase 532: N-substituted amide on heteroaromatic rings (IUPAC 2013 P-66.6.1).

name_heterocycle assembled the exocyclic amide suffix ("carboxamide") but
never collected N-substituents from the amide N, silently dropping them.
Same gap existed in _try_fused_hetero_retained for fused rings.
"""
import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Monosubstituted N-methyl
    ("CNC(=O)c1ccccn1",    "N-methylpyridine-2-carboxamide"),
    ("CNC(=O)c1cccnc1",    "N-methylpyridine-3-carboxamide"),
    ("CNC(=O)c1ccncc1",    "N-methylpyridine-4-carboxamide"),
    ("CNC(=O)c1cccs1",     "N-methylthiophene-2-carboxamide"),
    ("CNC(=O)c1ccco1",     "N-methylfuran-2-carboxamide"),
    ("CNC(=O)c1ccc[nH]1",  "N-methyl-1H-pyrrole-2-carboxamide"),
    # N-ethyl and N,N-dimethyl
    ("CCNC(=O)c1ccccn1",   "N-ethylpyridine-2-carboxamide"),
    ("CN(C)C(=O)c1ccccn1", "N,N-dimethylpyridine-2-carboxamide"),
    ("CN(C)C(=O)c1cccs1",  "N,N-dimethylthiophene-2-carboxamide"),
    # Regression: unsubstituted heteroaromatic amides unchanged
    ("NC(=O)c1ccccn1",     "pyridine-2-carboxamide"),
    ("NC(=O)c1ccc[nH]1",   "1H-pyrrole-2-carboxamide"),
    # Regression: benzene and aliphatic unchanged
    ("CNC(=O)c1ccccc1",    "N-methylbenzamide"),
    ("CN(C)C(=O)c1ccccc1", "N,N-dimethylbenzamide"),
    ("CNC(=O)C",           "N-methylacetamide"),
])
def test_phase532_n_substituted_heteroaromatic_amide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
