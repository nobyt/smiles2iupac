"""Phase 525: Heteroaromatic carbonyl halide naming (IUPAC 2013 P-65.5.1).

Pyridine, thiophene, furan, and pyrrole rings bearing -C(=O)X groups
were named as "(heteroaryl)methanoyl halide" (chain-based). Now correctly
named as "<ring>-N-carbonyl halide" following the carboxylic acid locant
pattern.
"""
import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Pyridine carbonyl chlorides
    ("O=C(Cl)c1ccccn1",   "pyridine-2-carbonyl chloride"),
    ("O=C(Cl)c1cccnc1",   "pyridine-3-carbonyl chloride"),
    ("O=C(Cl)c1ccncc1",   "pyridine-4-carbonyl chloride"),
    # Other heteroaromatic carbonyl chlorides
    ("O=C(Cl)c1cccs1",    "thiophene-2-carbonyl chloride"),
    ("O=C(Cl)c1ccco1",    "furan-2-carbonyl chloride"),
    ("O=C(Cl)c1ccc[nH]1", "1H-pyrrole-2-carbonyl chloride"),
    # Other halides
    ("O=C(F)c1ccccn1",    "pyridine-2-carbonyl fluoride"),
    ("O=C(Br)c1ccccn1",   "pyridine-2-carbonyl bromide"),
    ("O=C(I)c1ccccn1",    "pyridine-2-carbonyl iodide"),
    # Regression: benzene uses retained name
    ("O=C(Cl)c1ccccc1",   "benzoyl chloride"),
    ("O=C(F)c1ccccc1",    "benzoyl fluoride"),
    # Regression: aliphatic acyl halides unchanged
    ("O=C(Cl)C",          "acetyl chloride"),
    ("O=C(Cl)CC",         "propanoyl chloride"),
    ("O=C(Cl)CCC",        "butanoyl chloride"),
])
def test_phase525_heteroaromatic_acyl_halide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
