"""Phase 534: Heteroaromatic carboperoxoic acid naming (IUPAC 2013 P-65.1.1).

_name_peroxyacid had a benzene-only early return; heteroaromatic rings fell
through to _collect_acid_chain which yielded a 1-carbon chain, giving
"methaneperoxoic acid". Now correctly named as "<ring>-carboperoxoic acid".
"""
import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Heteroaromatic carboperoxoic acids
    ("OOC(=O)c1ccccn1",   "pyridine-2-carboperoxoic acid"),
    ("OOC(=O)c1cccnc1",   "pyridine-3-carboperoxoic acid"),
    ("OOC(=O)c1ccncc1",   "pyridine-4-carboperoxoic acid"),
    ("OOC(=O)c1cccs1",    "thiophene-2-carboperoxoic acid"),
    ("OOC(=O)c1ccco1",    "furan-2-carboperoxoic acid"),
    ("OOC(=O)c1ccc[nH]1", "1H-pyrrole-2-carboperoxoic acid"),
    # Regression: benzene uses retained-style name
    ("OOC(=O)c1ccccc1",   "benzeneperoxoic acid"),
    # Regression: aliphatic peroxyacids unchanged
    ("OOC(=O)C",          "ethaneperoxoic acid"),
    ("OOC(=O)CC",         "propaneperoxoic acid"),
])
def test_phase534_heteroaromatic_carboperoxoic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
