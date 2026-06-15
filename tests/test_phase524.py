"""Phase 524: Heteroaromatic sulfonic acid / sulfonamide / sulfonyl halide naming.

Pyridine, thiophene, furan, and pyrrole rings bearing -S(=O)(=O)OH,
-S(=O)(=O)NH2, or -S(=O)(=O)Cl groups were being treated as aliphatic
chains. Now correctly named using the retained heteroaromatic ring name
with an explicit locant (IUPAC 2013 P-65.3.1, P-65.3.2).
"""
import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Pyridine sulfonic acids
    ("OS(=O)(=O)c1ccccn1",   "pyridine-2-sulfonic acid"),
    ("OS(=O)(=O)c1cccnc1",   "pyridine-3-sulfonic acid"),
    ("OS(=O)(=O)c1ccncc1",   "pyridine-4-sulfonic acid"),
    # Heteroaromatic 5-membered sulfonic acids
    ("OS(=O)(=O)c1cccs1",    "thiophene-2-sulfonic acid"),
    ("OS(=O)(=O)c1ccco1",    "furan-2-sulfonic acid"),
    ("OS(=O)(=O)c1ccc[nH]1", "1H-pyrrole-2-sulfonic acid"),
    # Heteroaromatic sulfonamides
    ("NS(=O)(=O)c1ccccn1",   "pyridine-2-sulfonamide"),
    ("NS(=O)(=O)c1cccs1",    "thiophene-2-sulfonamide"),
    # Heteroaromatic sulfonyl chlorides
    ("ClS(=O)(=O)c1ccccn1",  "pyridine-2-sulfonyl chloride"),
    ("ClS(=O)(=O)c1cccs1",   "thiophene-2-sulfonyl chloride"),
    # Regression: benzene derivatives unchanged
    ("OS(=O)(=O)c1ccccc1",   "benzenesulfonic acid"),
    ("NS(=O)(=O)c1ccccc1",   "benzenesulfonamide"),
    ("ClS(=O)(=O)c1ccccc1",  "benzenesulfonyl chloride"),
    # Regression: aliphatic sulfonic acids unchanged
    ("OS(=O)(=O)C",    "methanesulfonic acid"),
    ("OS(=O)(=O)CC",   "ethanesulfonic acid"),
])
def test_phase524_heteroaromatic_sulfonyl(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
