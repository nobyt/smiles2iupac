"""Phase 528: Heteroaromatic carboximidic acid naming (IUPAC 2013 P-65.1.2.4).

Rings bearing -C(=NH)OH were returning "methanimidic acid" because
_name_imidic_acid used _collect_acid_chain which stops at the aromatic ring,
yielding a 1-carbon chain. Now correctly named as "<ring>-N-carboximidic acid".
Also fixed benzene: "benzenecarboximidic acid" (spec.benzene_name was None).
"""
import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Heteroaromatic carboximidic acids
    ("OC(=N)c1ccccn1",   "pyridine-2-carboximidic acid"),
    ("OC(=N)c1cccnc1",   "pyridine-3-carboximidic acid"),
    ("OC(=N)c1ccncc1",   "pyridine-4-carboximidic acid"),
    ("OC(=N)c1cccs1",    "thiophene-2-carboximidic acid"),
    ("OC(=N)c1ccco1",    "furan-2-carboximidic acid"),
    ("OC(=N)c1ccc[nH]1", "1H-pyrrole-2-carboximidic acid"),
    # N-substituted
    ("CN=C(O)c1ccccn1",  "N-methylpyridine-2-carboximidic acid"),
    ("CN=C(O)c1ccccc1",  "N-methylbenzenecarboximidic acid"),
    # Regression: benzene now produces retained-style name
    ("OC(=N)c1ccccc1",   "benzenecarboximidic acid"),
    # Regression: aliphatic imidic acids unchanged
    ("OC(=N)C",          "ethanimidic acid"),
    ("OC(=N)CC",         "propanimidic acid"),
])
def test_phase528_heteroaromatic_carboximidic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
