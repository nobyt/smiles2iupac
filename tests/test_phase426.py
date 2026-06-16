"""Phase 426: Thianthrene and benzo[a]pyrene (corrected).

IUPAC 2013 P-31.1.3: retained name for the S,S-bridged dibenzene
ring system and the 5-ring PAH C20H12.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # thianthrene — C12H8S2, two benzene rings bridged by two S atoms
    ("c1ccc2c(c1)Sc1ccccc1S2",            "thianthrene"),
    # benzo[a]pyrene — C20H12, five fused 6-membered rings (correct SMILES)
    ("c1ccc2c(c1)cc1ccc3cccc4ccc2c1c34",  "benzo[a]pyrene"),
    # regression: phenoxathiin unchanged (Phase 423) — O and S bridges
    ("c1ccc2c(c1)Oc1ccccc1S2",             "phenoxathiin"),
    # regression: thioxanthene unchanged (Phase 134) — one S bridge + CH2
    ("c1ccc2c(c1)Cc1ccccc1S2",             "thioxanthene"),
    # regression: pyrene unchanged
    ("c1cc2ccc3cccc4ccc(c1)c2c34",          "pyrene"),
    # regression: benzene unchanged
    ("c1ccccc1",                              "benzene"),
])
def test_phase426_thianthrene_benzoapyrene(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
