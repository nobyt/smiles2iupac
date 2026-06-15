"""Phase 527: Heteroaromatic carboselenoamide / carbohydrazide / carboximidamide naming.

Three functional groups that were misrouted to the heterocycle path for
heteroaromatic ring attachment, producing wrong names:
- selenoamide: same return-None guard as thioamide → ring path used
- hydrazide: break after benzene check → chain naming used
- amidine: only handled in _name_acyclic (not reached for cyclic molecules)

All now produce "<ring>-N-carbo{suffix}" names per IUPAC 2013.
"""
import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Heteroaromatic carboselenoamides
    ("[Se]=C(N)c1ccccn1",  "pyridine-2-carboselenoamide"),
    ("[Se]=C(N)c1cccnc1",  "pyridine-3-carboselenoamide"),
    ("[Se]=C(N)c1cccs1",   "thiophene-2-carboselenoamide"),
    ("[Se]=C(N)c1ccco1",   "furan-2-carboselenoamide"),
    # Heteroaromatic carbohydrazides
    ("NNC(=O)c1ccccn1",    "pyridine-2-carbohydrazide"),
    ("NNC(=O)c1cccnc1",    "pyridine-3-carbohydrazide"),
    ("NNC(=O)c1cccs1",     "thiophene-2-carbohydrazide"),
    ("NNC(=O)c1ccco1",     "furan-2-carbohydrazide"),
    # Heteroaromatic carboximidamides
    ("NC(=N)c1ccccn1",     "pyridine-2-carboximidamide"),
    ("NC(=N)c1cccnc1",     "pyridine-3-carboximidamide"),
    ("NC(=N)c1ccncc1",     "pyridine-4-carboximidamide"),
    ("NC(=N)c1cccs1",      "thiophene-2-carboximidamide"),
    ("NC(=N)c1ccco1",      "furan-2-carboximidamide"),
    # Regression: benzene uses retained names
    ("[Se]=C(N)c1ccccc1",  "benzeneselenoamide"),
    ("O=C(NN)c1ccccc1",    "benzohydrazide"),
    ("NC(=N)c1ccccc1",     "benzenecarboximidamide"),
    # Regression: aliphatic forms unchanged
    ("[Se]=C(N)C",         "ethaneselenoamide"),
    ("NNC(=O)C",           "ethanohydrazide"),
    ("NC(=N)CC",           "propanimidamide"),
    ("CNC(=N)C",           "N-methylethanimidamide"),
])
def test_phase527_heteroaromatic_carbo_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
