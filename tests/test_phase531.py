"""Phase 531: Heteroaromatic imidate ester naming (IUPAC 2013 P-65.1.2.4).

Rings bearing -C(=NH)OR were returning "alkyl methanimidate" because
_name_imidate_ester used _collect_acid_chain which yielded a 1-carbon chain.
Also fixes benzene: "alkyl benzenecarboximidate" (spec.benzene_name was None).
"""
import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Heteroaromatic imidate esters
    ("CCOC(=N)c1ccccn1",  "ethyl pyridine-2-carboximidate"),
    ("CCOC(=N)c1cccnc1",  "ethyl pyridine-3-carboximidate"),
    ("CCOC(=N)c1ccncc1",  "ethyl pyridine-4-carboximidate"),
    ("CCOC(=N)c1cccs1",   "ethyl thiophene-2-carboximidate"),
    ("CCOC(=N)c1ccco1",   "ethyl furan-2-carboximidate"),
    ("COC(=N)c1ccc[nH]1", "methyl 1H-pyrrole-2-carboximidate"),
    # Regression: benzene now produces systematic name
    ("CCOC(=N)c1ccccc1",  "ethyl benzenecarboximidate"),
    ("COC(=N)c1ccccc1",   "methyl benzenecarboximidate"),
    # Regression: aliphatic imidate esters unchanged
    ("CCOC(=N)C",         "ethyl ethanimidate"),
    ("COC(=N)CC",         "methyl propanimidate"),
])
def test_phase531_heteroaromatic_imidate_ester(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
