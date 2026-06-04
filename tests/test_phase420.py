"""Phase 420: 2H-Pyran-2-one (alpha-pyrone) and 4H-pyran-4-one (gamma-pyrone).

IUPAC 2013: retained/systematic names for the monocyclic O-containing
lactone and vinylogous ester ring systems.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 2H-pyran-2-one (alpha-pyrone) — O1 adjacent to C2=O
    ("O=c1cccco1",              "2H-pyran-2-one"),
    # 4H-pyran-4-one (gamma-pyrone) — O1 at para to C4=O
    ("O=c1ccocc1",              "4H-pyran-4-one"),
    # regression: coumarin unchanged (2H-pyran-2-one fused with benzene)
    ("O=c1ccc2ccccc2o1",        "coumarin"),
    # regression: chromone unchanged (4H-pyran-4-one fused with benzene)
    ("O=c1ccoc2ccccc12",        "chromone"),
    # regression: benzene unchanged
    ("c1ccccc1",                 "benzene"),
    # regression: naphthalen-2(1H)-one unchanged (Phase 419)
    ("O=C1C=Cc2ccccc2C1",       "naphthalen-2(1H)-one"),
])
def test_phase420_pyranones(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
