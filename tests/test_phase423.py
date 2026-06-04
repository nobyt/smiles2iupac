"""Phase 423: Phenoxathiin and triphenylene.

IUPAC 2013 P-31.1.3: retained names for the O/S bridged dibenzene ring
system and the symmetric tetracyclic hydrocarbon.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # phenoxathiin — dibenzo[b,e][1,4]oxathiin, O and S bridge two benzene rings
    ("c1ccc2c(c1)Oc1ccccc1S2",          "phenoxathiin"),
    # triphenylene — three benzene rings fused in a symmetrical triangle
    ("c1ccc2c(c1)c1ccccc1c1ccccc21",    "triphenylene"),
    # regression: xanthene unchanged (Phase 134) — O bridge, no S
    ("c1ccc2c(c1)Cc1ccccc1O2",           "xanthene"),
    # regression: thioxanthene unchanged (Phase 134) — S bridge, no O
    ("c1ccc2c(c1)Cc1ccccc1S2",           "thioxanthene"),
    # regression: pyrene unchanged (Phase 422)
    ("c1cc2ccc3cccc4ccc(c1)c2c34",       "pyrene"),
    # regression: benzene unchanged
    ("c1ccccc1",                           "benzene"),
])
def test_phase423_phenoxathiin_triphenylene(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
