"""Phase 540: Nitroso on heteroaromatic rings (IUPAC 2013 P-62.3.1).

_name_nitroso had a benzene-only early return that fell through to
"nitrosobenzene" for ALL aromatic rings, including pyridine/thiophene/furan.
Now uses _aryl_sulfonyl_prefix for heteroaromatic detection; benzene unchanged.
"""
import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Nitroso on heteroaromatic rings
    ("O=Nc1ccccn1",    "2-nitrosopyridine"),
    ("O=Nc1cccnc1",    "3-nitrosopyridine"),
    ("O=Nc1ccncc1",    "4-nitrosopyridine"),
    ("O=Nc1cccs1",     "2-nitrosothiophene"),
    ("O=Nc1ccco1",     "2-nitrosofuran"),
    ("O=Nc1ccc[nH]1",  "2-nitroso-1H-pyrrole"),
    # Regression: benzene and aliphatic unchanged
    ("O=Nc1ccccc1",    "nitrosobenzene"),
    ("O=NC",           "nitrosomethane"),
    ("O=NCC",          "nitrosoethane"),
])
def test_phase540_nitroso_heteroaromatic(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
