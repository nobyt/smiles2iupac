"""Phase 535: Heteroaromatic O-thioate ester naming (IUPAC 2013 P-65.1.6).

_name_o_thioester had no ring detection; _collect_acid_chain returned a
1-carbon chain for ring-adjacent C(=S), giving "O-alkyl methanethioate".
Now uses _aryl_sulfonyl_prefix for both benzene and heteroaromatic.
"""
import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Heteroaromatic O-thioate esters
    ("CCOC(=S)c1ccccn1",  "O-ethyl pyridine-2-carbothioate"),
    ("CCOC(=S)c1cccnc1",  "O-ethyl pyridine-3-carbothioate"),
    ("CCOC(=S)c1ccncc1",  "O-ethyl pyridine-4-carbothioate"),
    ("CCOC(=S)c1cccs1",   "O-ethyl thiophene-2-carbothioate"),
    ("CCOC(=S)c1ccco1",   "O-ethyl furan-2-carbothioate"),
    # Benzene now also fixed
    ("CCOC(=S)c1ccccc1",  "O-ethyl benzenecarbothioate"),
    # Regression: aliphatic O-thioate esters unchanged
    ("CCOC(=S)C",         "O-ethyl ethanethioate"),
    ("COC(=S)CC",         "O-methyl propanethioate"),
])
def test_phase535_heteroaromatic_o_thioate_ester(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
