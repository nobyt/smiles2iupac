"""Phase 404: N-substituted non-aromatic fused ring systems (IUPAC 2013 P-31.1.2).

Add _FUSED_LOCANT_MAP entries for indoline and 1,2,3,4-tetrahydroquinoline so
that _try_fused_hetero_retained can handle substituted versions.  Also fix
_collect_hetero_substituents to skip ring-junction bonds (fused aromatic
ring atoms should not appear as phantom phenyl substituents).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # N-methylindoline (N at pos 1 → numeric locant)
    ("CN1CCc2ccccc21",                "1-methyl-2,3-dihydroindole"),
    # C-substituted indoline (methyl on benzene ring)
    ("Cc1ccc2c(c1)NCC2",              "6-methyl-2,3-dihydro-1H-indole"),
    ("Cc1ccc2c(c1)CCN2",              "5-methyl-2,3-dihydro-1H-indole"),
    # N-methyl-1,2,3,4-tetrahydroquinoline
    ("CN1CCCc2ccccc21",               "1-methyl-1,2,3,4-tetrahydroquinoline"),
    # regression: unsubstituted indoline unchanged
    ("c1ccc2c(c1)CCN2",               "2,3-dihydro-1H-indole"),
    # regression: 1,2,3,4-tetrahydroquinoline unchanged
    ("c1ccc2c(c1)CCCN2",              "1,2,3,4-tetrahydroquinoline"),
    # regression: 1-methylindole (aromatic) — Phase 845: N1 substituted drops 1H-
    ("Cn1ccc2ccccc21",                "1-methylindole"),
    # regression: piperidine unchanged
    ("C1CCNCC1",                      "piperidine"),
    # regression: N-methylpiperidine unchanged
    ("CN1CCCCC1",                     "N-methylpiperidine"),
])
def test_phase404_fused_ring_n_sub(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
