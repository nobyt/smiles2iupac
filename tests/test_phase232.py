"""Phase 232: 1,3,5-trioxane retained name (IUPAC 2013 P-31.1.3.4)."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("C1OCOCO1", "1,3,5-trioxane"),
    # regression: 1,4-dioxane still works
    ("C1COCCO1", "1,4-dioxane"),
])
def test_phase232_trioxane(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
