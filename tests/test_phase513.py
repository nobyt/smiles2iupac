"""Phase 513: fix purine tautomer indicated-H positions (7H vs 9H vs 6H)
(IUPAC 2013 P-31.1.6: purine numbering, N9 adjacent to C4/N3-side ring junction).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Standard purine input (NH at N9 position, canonical -> 9H-purine)
    ("c1ncnc2[nH]cnc12",  "9H-purine"),
    # Canonical 9H-purine form (NH adjacent to C4 bridgehead)
    ("c1ncc2nc[nH]c2n1",  "9H-purine"),
    # 7H-purine tautomer (NH adjacent to C5 bridgehead)
    ("c1ncc2[nH]cnc2n1",  "7H-purine"),
    # Adenine (amino-substituted 9H-purine) unaffected
    ("Nc1ncnc2[nH]cnc12", "adenine"),
])
def test_phase513_purine_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
