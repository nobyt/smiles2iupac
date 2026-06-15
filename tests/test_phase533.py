"""Phase 533: Heteroaromatic carbothioate ester (S-ester) naming (IUPAC 2013 P-65.1.6).

_name_thioester had no ring detection at all; _collect_acid_chain returned
a 1-carbon chain for ring-adjacent carbonyl C, giving "S-alkyl methanethioate"
for both benzene and heteroaromatic substrates.
Now both use _aryl_sulfonyl_prefix → "S-alkyl {ring}carbothioate".
"""
import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Heteroaromatic S-thioate esters
    ("CCSC(=O)c1ccccn1",  "S-ethyl pyridine-2-carbothioate"),
    ("CCSC(=O)c1cccnc1",  "S-ethyl pyridine-3-carbothioate"),
    ("CCSC(=O)c1ccncc1",  "S-ethyl pyridine-4-carbothioate"),
    ("CSC(=O)c1cccs1",    "S-methyl thiophene-2-carbothioate"),
    ("CCSC(=O)c1ccco1",   "S-ethyl furan-2-carbothioate"),
    ("CSC(=O)c1ccc[nH]1", "S-methyl 1H-pyrrole-2-carbothioate"),
    # Benzene now also fixed
    ("CCSC(=O)c1ccccc1",  "S-ethyl benzenecarbothioate"),
    ("CSC(=O)c1ccccc1",   "S-methyl benzenecarbothioate"),
    # Regression: aliphatic S-thioate esters unchanged
    ("CCSC(=O)C",         "S-ethyl ethanethioate"),
    ("CSC(=O)CC",         "S-methyl propanethioate"),
])
def test_phase533_heteroaromatic_s_thioate_ester(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
