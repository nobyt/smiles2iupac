"""Phase 390: Cyclic N in carbamate → {alkyl} {ring}-{N-locant}-carboxylate (IUPAC 2013 P-65.6).

When the N of a carbamate is part of a ring (piperidine, pyrrolidine,
morpholine, etc.) the compound is named as an N-heterocycle carboxylate,
not as an N,N-disubstituted carbamate.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # piperidine (retained name, N at 1)
    ("CCOC(=O)N1CCCCC1",    "ethyl piperidine-1-carboxylate"),
    ("COC(=O)N1CCCCC1",     "methyl piperidine-1-carboxylate"),
    # pyrrolidine (HW name, N at 1)
    ("COC(=O)N1CCCC1",      "methyl pyrrolidine-1-carboxylate"),
    ("CCOC(=O)N1CCCC1",     "ethyl pyrrolidine-1-carboxylate"),
    # morpholine (retained name, N at 4)
    ("CCOC(=O)N1CCOCC1",    "ethyl morpholine-4-carboxylate"),
    ("COC(=O)N1CCOCC1",     "methyl morpholine-4-carboxylate"),
])
def test_phase390_cyclic_carbamate(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
