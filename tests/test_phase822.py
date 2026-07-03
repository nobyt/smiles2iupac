"""Phase 822: tetrazolo[1,5-a]pyrazine and tetrazolo[1,5-d][1,2,4]triazine α-ol/thiol → tautomers."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("Oc1cncc2nnnn12",   "tetrazolo[1,5-a]pyrazin-5(1H)-one"),
    ("Sc1cncc2nnnn12",   "tetrazolo[1,5-a]pyrazin-5(1H)-thione"),
    ("Oc1cn2nnnc2cn1",   "tetrazolo[1,5-a]pyrazin-6(4H)-one"),
    ("Sc1cn2nnnc2cn1",   "tetrazolo[1,5-a]pyrazin-6(4H)-thione"),
    ("Oc1nccn2nnnc12",   "tetrazolo[1,5-a]pyrazin-8(5H)-one"),
    ("Sc1nccn2nnnc12",   "tetrazolo[1,5-a]pyrazin-8(5H)-thione"),
    ("Oc1nncc2nnnn12",   "tetrazolo[1,5-d][1,2,4]triazin-5(1H)-one"),
    ("Sc1nncc2nnnn12",   "tetrazolo[1,5-d][1,2,4]triazin-5(1H)-thione"),
    ("Oc1nncn2nnnc12",   "tetrazolo[1,5-d][1,2,4]triazin-8(5H)-one"),
    ("Sc1nncn2nnnc12",   "tetrazolo[1,5-d][1,2,4]triazin-8(5H)-thione"),
    ("c1cn2nnnc2cn1",    "tetrazolo[1,5-a]pyrazine"),
    ("c1nncn2nnnc12",    "tetrazolo[1,5-d][1,2,4]triazine"),
])
def test_phase822(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
