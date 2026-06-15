"""Phase 539: N-hydroxy amide (hydroxamic acid) on heteroaromatic rings.

The amide N-substituent fix in phase 532 collected only C-neighbors of N.
N-OH (hydroxamic acid) was silently dropped for heteroaromatic rings.
Now "hydroxy" is included in the prefix list alongside C-substituents.
"""
import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # N-hydroxy heteroaromatic amides (hydroxamic acids)
    ("ONC(=O)c1ccccn1",   "N-hydroxypyridine-2-carboxamide"),
    ("ONC(=O)c1cccnc1",   "N-hydroxypyridine-3-carboxamide"),
    ("ONC(=O)c1ccncc1",   "N-hydroxypyridine-4-carboxamide"),
    ("ONC(=O)c1cccs1",    "N-hydroxythiophene-2-carboxamide"),
    ("ONC(=O)c1ccco1",    "N-hydroxyfuran-2-carboxamide"),
    ("ONC(=O)c1ccc[nH]1", "N-hydroxy-1H-pyrrole-2-carboxamide"),
    # Regression: unsubstituted heteroaromatic amide unchanged
    ("NC(=O)c1ccccn1",    "pyridine-2-carboxamide"),
    # Regression: benzene and aliphatic unchanged
    ("ONC(=O)c1ccccc1",   "N-hydroxybenzamide"),
    ("ONC(=O)C",          "N-hydroxyacetamide"),
])
def test_phase539_n_hydroxy_heteroaromatic_amide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
