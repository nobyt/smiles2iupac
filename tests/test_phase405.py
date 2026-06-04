"""Phase 405: Phthalimide derivatives → isoindole-1,3(2H)-dione names.

When a 5-membered pyrrolidine-dione ring is fused with an aromatic ring
(the 2 junction atoms are aromatic), use the isoindole-1,3(2H)-dione parent
and give N its IUPAC locant 2.

IUPAC 2013 P-31.1.3: isoindole-1,3(2H)-dione (phthalimide).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # N-H phthalimide (already handled by __init__ lookup → unchanged)
    ("O=C1NC(=O)c2ccccc21",           "phthalimide"),
    # N-methylphthalimide
    ("CN1C(=O)c2ccccc2C1=O",          "2-methylisoindole-1,3(2H)-dione"),
    # N-ethylphthalimide
    ("CCN1C(=O)c2ccccc2C1=O",         "2-ethylisoindole-1,3(2H)-dione"),
    # regression: succinimide (no fused aromatic) unchanged
    ("O=C1CCC(=O)N1",                 "pyrrolidine-2,5-dione"),
    # regression: maleimide (no fused aromatic) unchanged
    ("O=C1C=CC(=O)N1",                "1H-pyrrole-2,5-dione"),
    # regression: piperidine-2,6-dione (glutarimide) unchanged
    ("O=C1CCCC(=O)N1",                "piperidine-2,6-dione"),
    # regression: 1-methylindoline unchanged
    ("CN1CCc2ccccc21",                "1-methylindoline"),
    # regression: benzene unchanged
    ("c1ccccc1",                       "benzene"),
])
def test_phase405_phthalimide_dione(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
