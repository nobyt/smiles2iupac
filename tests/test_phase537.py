"""Phase 537: Heteroaromatic carboperoxoate ester naming (IUPAC 2013 P-65.1.1).

_name_peroxyester had no ring detection; _collect_acid_chain returned a
1-carbon chain for ring-adjacent carbonyl C, giving wrong names. Now uses
_aryl_sulfonyl_prefix for both benzene and heteroaromatic.
"""
import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Heteroaromatic carboperoxoate esters
    ("CCOOC(=O)c1ccccn1",  "ethyl pyridine-2-carboperoxoate"),
    ("CCOOC(=O)c1cccnc1",  "ethyl pyridine-3-carboperoxoate"),
    ("CCOOC(=O)c1ccncc1",  "ethyl pyridine-4-carboperoxoate"),
    ("CCOOC(=O)c1cccs1",   "ethyl thiophene-2-carboperoxoate"),
    ("CCOOC(=O)c1ccco1",   "ethyl furan-2-carboperoxoate"),
    # Benzene now also fixed
    ("CCOOC(=O)c1ccccc1",  "ethyl benzenecarboperoxoate"),
    # Regression: aliphatic peroxyesters unchanged
    ("CCOOC(=O)C",         "ethyl ethaneperoxoate"),
    ("CCOOC(=O)CC",        "ethyl propaneperoxoate"),
])
def test_phase537_heteroaromatic_carboperoxoate_ester(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
