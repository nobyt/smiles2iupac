"""Phase 536: Heteroaromatic carbodithioate ester naming (IUPAC 2013 P-65.1.6).

_name_s_dithioate_ester had no ring detection; _collect_acid_chain returned
a 1-carbon chain, giving "S-alkyl methanedithioate". Now uses
_aryl_sulfonyl_prefix for both benzene and heteroaromatic.
"""
import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Heteroaromatic S-dithioate esters
    ("CCSC(=S)c1ccccn1",  "S-ethyl pyridine-2-carbodithioate"),
    ("CCSC(=S)c1cccnc1",  "S-ethyl pyridine-3-carbodithioate"),
    ("CSC(=S)c1cccs1",    "S-methyl thiophene-2-carbodithioate"),
    ("CCSC(=S)c1ccco1",   "S-ethyl furan-2-carbodithioate"),
    # Benzene now also fixed
    ("CCSC(=S)c1ccccc1",  "S-ethyl benzenecarbodithioate"),
    # Regression: aliphatic S-dithioate esters unchanged
    ("CSC(=S)C",          "S-methyl ethanedithioate"),
    ("CCSC(=S)CC",        "S-ethyl propanedithioate"),
])
def test_phase536_heteroaromatic_s_dithioate_ester(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
