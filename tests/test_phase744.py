"""Phase 744: 2- and 4-hydroxypyridine → 1H-pyridin-2/4-one (IUPAC 2013 PINs).

IUPAC 2013 preferred tautomers: α- and γ-hydroxypyridines and
2-hydroxypyrimidine prefer the lactam (keto) form over the enol form.
β-hydroxypyridine (position 3) is not in this class and stays as -ol.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # α-hydroxypyridine → lactam (1H-pyridin-2-one)
    ("Oc1ccccn1",    "1H-pyridin-2-one"),
    # γ-hydroxypyridine → lactam (1H-pyridin-4-one)
    ("Oc1ccncc1",    "1H-pyridin-4-one"),
    # 2-hydroxypyrimidine → lactam
    ("Oc1ncccn1",    "1H-pyrimidin-2-one"),
    # Substituted: substituents become prefix of the lactam base
    ("Cc1ccnc(O)c1", "4-methyl-1H-pyridin-2-one"),
    # Regression: β-hydroxypyridine stays as -ol
    ("Oc1cccnc1",    "pyridin-3-ol"),
    # Regression: furanols unaffected
    ("Oc1ccco1",     "furan-2-ol"),
    # Regression: 4-aminopyridine stays as -amine
    ("Nc1ccncc1",    "pyridin-4-amine"),
])
def test_phase744_hydroxypyridine_tautomer(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
