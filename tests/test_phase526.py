"""Phase 526: Heteroaromatic carbothioamide naming (IUPAC 2013 P-65.1.2.4).

Pyridine, thiophene, furan, and pyrrole rings bearing -C(=S)NH2 groups
were routed to the heterocycle name_heterocycle path, producing wrong names
like "2-aminomethylpyridine". Now correctly named as "<ring>-N-carbothioamide".
"""
import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Pyridine carbothioamides
    ("S=C(N)c1ccccn1",    "pyridine-2-carbothioamide"),
    ("S=C(N)c1cccnc1",    "pyridine-3-carbothioamide"),
    ("S=C(N)c1ccncc1",    "pyridine-4-carbothioamide"),
    # Other heteroaromatic carbothioamides
    ("S=C(N)c1cccs1",     "thiophene-2-carbothioamide"),
    ("S=C(N)c1ccco1",     "furan-2-carbothioamide"),
    ("S=C(N)c1ccc[nH]1",  "1H-pyrrole-2-carbothioamide"),
    # N-substituted
    ("S=C(NC)c1ccccn1",   "N-methylpyridine-2-carbothioamide"),
    ("S=C(N(C)C)c1ccccn1", "N,N-dimethylpyridine-2-carbothioamide"),
    # Regression: benzene uses retained name
    ("S=C(N)c1ccccc1",    "benzothioamide"),
    # Regression: aliphatic thioamides unchanged
    ("S=C(N)C",           "ethanethioamide"),
    ("S=C(NC)C",          "N-methylethanethioamide"),
    # Regression: cycloalkane thioamide unchanged
    ("NC(=S)C1CCCCC1",    "cyclohexanecarbothioamide"),
])
def test_phase526_heteroaromatic_thioamide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
