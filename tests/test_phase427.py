"""Phase 427: 1H-phenalene.

IUPAC 2013 P-31.1.3: retained name for the tricyclic peri-fused 3×6
hydrocarbon ring system with one sp3 CH2 (C13H10).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1H-phenalene — C13H10, three 6-membered rings sharing one sp3 CH2 bridge
    ("C1=Cc2cccc3cccc(c23)C1",   "1H-phenalene"),
    # regression: acenaphthylene unchanged (Phase 421) — similar 3-ring system with 5-membered ring
    ("C1=Cc2cccc3cccc1c23",       "acenaphthylene"),
    # regression: acenaphthene unchanged (Phase 421)
    ("c1cc2c3c(cccc3c1)CC2",      "acenaphthene"),
    # regression: naphthalene unchanged — two rings
    ("c1ccc2ccccc2c1",             "naphthalene"),
    # regression: fluorene unchanged (Phase 134) — three rings with CH2 bridge
    ("c1ccc2c(c1)Cc1ccccc1-2",    "fluorene"),
    # regression: benzene unchanged
    ("c1ccccc1",                    "benzene"),
])
def test_phase427_phenalene(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
