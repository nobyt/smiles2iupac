"""Phase 541: Aldehyde oxime (carbaldehyde oxime) on heteroaromatic rings.

name_heterocycle's _EXOCYCLIC_SUFFIX only covered carboxylic_acid/aldehyde/
amide/nitrile; aldoxime fell through to substituent naming and was
misidentified as "aminomethyl". Adding "aldoxime": "carbaldehyde oxime" routes
it through the same locant-aware exocyclic path as the other carbonyl suffixes.
"""
import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Aldehyde oxime on heteroaromatic rings
    ("ON=Cc1ccccn1",    "pyridine-2-carbaldehyde oxime"),
    ("ON=Cc1cccnc1",    "pyridine-3-carbaldehyde oxime"),
    ("ON=Cc1ccncc1",    "pyridine-4-carbaldehyde oxime"),
    ("ON=Cc1cccs1",     "thiophene-2-carbaldehyde oxime"),
    ("ON=Cc1ccco1",     "furan-2-carbaldehyde oxime"),
    ("ON=Cc1ccc[nH]1",  "1H-pyrrole-2-carbaldehyde oxime"),
    # Regression: benzene and aliphatic unchanged
    ("ON=Cc1ccccc1",    "benzaldoxime"),
    ("ON=CC",           "ethanal oxime"),
    # Regression: ketoxime on heteroaromatic unchanged
    ("CC(=NO)c1ccccn1", "1-(pyridin-2-yl)ethan-1-one oxime"),
])
def test_phase541_aldoxime_heteroaromatic(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
