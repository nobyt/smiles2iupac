"""Phase 403: Maleimide → 1H-pyrrole-2,5-dione (IUPAC 2013 P-31.1.7).

When a 5-membered N-ring with two flanking C=O groups (imide pattern) is
detected via partial-unsaturation, rename from the dihydro/trihydropyrrole
form to the 1H-pyrrole parent with indicated-H notation.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # maleimide (N-H)
    ("O=C1C=CC(=O)N1",                "1H-pyrrole-2,5-dione"),
    # N-phenylmaleimide
    ("O=C1C=CC(=O)N1c1ccccc1",        "1-phenyl-1H-pyrrole-2,5-dione"),
    # regression: succinimide (no C=C in ring) → pyrrolidine-2,5-dione
    ("O=C1CCC(=O)N1",                 "pyrrolidine-2,5-dione"),
    # regression: glutarimide (6-membered) → piperidine-2,6-dione
    ("O=C1CCCC(=O)N1",                "piperidine-2,6-dione"),
    # regression: caprolactam (1 C=O only) → azepan-2-one
    ("O=C1CCCCCN1",                   "azepan-2-one"),
    # regression: benzene unchanged
    ("c1ccccc1",                       "benzene"),
])
def test_phase403_maleimide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
