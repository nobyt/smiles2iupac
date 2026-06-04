"""Phase 425: Perylene.

IUPAC 2013 P-31.1.3: retained name for the pentacyclic PAH C20H12
(already in Phase 138 dict; this phase adds regression coverage).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # perylene — C20H12, five fused 6-membered rings (peri-fusion)
    ("C1=Cc2cccc3ccc4ccc5cccc1c5c4c23",  "perylene"),
    # regression: pyrene unchanged (Phase 138) — four rings
    ("c1cc2ccc3cccc4ccc(c1)c2c34",        "pyrene"),
    # regression: chrysene unchanged (Phase 138) — four rings, angular
    ("c1ccc2cc3c(ccc4ccccc43)cc2c1",      "chrysene"),
    # regression: tetracene unchanged (Phase 424) — four rings, linear
    ("c1ccc2cc3cc4ccccc4cc3cc2c1",        "tetracene"),
    # regression: naphthalene unchanged — two rings
    ("c1ccc2ccccc2c1",                     "naphthalene"),
    # regression: benzene unchanged
    ("c1ccccc1",                            "benzene"),
])
def test_phase425_perylene(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
