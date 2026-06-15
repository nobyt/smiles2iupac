"""Phase 538: Heteroaromatic carboxylic acid anhydride naming (IUPAC 2013 P-65.7).

_name_anhydride's inner _aryl_acid_name helper only handled benzene
("benzoic"); heteroaromatic rings fell through to a 1-carbon chain chain
producing "formic anhydride". Now uses _aryl_sulfonyl_prefix to generate
"pyridine-2-carboxylic" etc. for the acid stem.
"""
import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Symmetric heteroaromatic anhydrides
    ("O=C(OC(=O)c1ccccn1)c1ccccn1",  "pyridine-2-carboxylic anhydride"),
    ("O=C(OC(=O)c1cccnc1)c1cccnc1",  "pyridine-3-carboxylic anhydride"),
    ("O=C(OC(=O)c1ccncc1)c1ccncc1",  "pyridine-4-carboxylic anhydride"),
    ("O=C(OC(=O)c1cccs1)c1cccs1",    "thiophene-2-carboxylic anhydride"),
    ("O=C(OC(=O)c1ccco1)c1ccco1",    "furan-2-carboxylic anhydride"),
    # Mixed anhydrides
    ("CC(=O)OC(=O)c1ccccn1",         "acetic pyridine-2-carboxylic anhydride"),
    # Regression: benzene and aliphatic unchanged
    ("O=C(OC(=O)c1ccccc1)c1ccccc1",  "benzoic anhydride"),
    ("CC(=O)OC(=O)C",                 "acetic anhydride"),
])
def test_phase538_heteroaromatic_carboxylic_anhydride(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
