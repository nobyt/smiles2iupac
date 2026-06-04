"""Phase 424: Tetracene (naphthacene).

IUPAC 2013 P-31.1.3: retained name for the linear tetracyclic PAH C18H12.
Chrysene (angular C18H12) was already supported from Phase 138.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # tetracene (naphthacene) — C18H12, four fused 6-membered rings in a linear row
    ("c1ccc2cc3cc4ccccc4cc3cc2c1",      "tetracene"),
    # regression: chrysene unchanged (Phase 138) — angular arrangement
    ("c1ccc2cc3c(ccc4ccccc43)cc2c1",    "chrysene"),
    # regression: triphenylene unchanged (Phase 138) — symmetric triangle
    ("c1ccc2c(c1)c1ccccc1c1ccccc21",    "triphenylene"),
    # regression: pyrene unchanged (Phase 138)
    ("c1cc2ccc3cccc4ccc(c1)c2c34",       "pyrene"),
    # regression: anthracene unchanged — three rings linear
    ("c1ccc2cc3ccccc3cc2c1",             "anthracene"),
    # regression: benzene unchanged
    ("c1ccccc1",                           "benzene"),
])
def test_phase424_tetracene(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
