"""Phase 264: benzoic anhydride, isobenzofuran-1,3-dione, furan-2,5-dione (IUPAC 2013).

  O=C(OC(=O)c1ccccc1)c1ccccc1 → benzoic anhydride          (aryl anhydride fix)
  O=C1OC(=O)c2ccccc21          → isobenzofuran-1,3-dione    (phthalic anhydride)
  O=C1C=CC(=O)O1               → furan-2,5-dione            (maleic anhydride)
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # benzoic anhydride (symmetric aryl anhydride)
    ("O=C(OC(=O)c1ccccc1)c1ccccc1",  "benzoic anhydride"),
    # mixed aryl/alkyl anhydride
    ("O=C(OC(=O)C)c1ccccc1",          "acetic benzoic anhydride"),
    # isobenzofuran-1,3-dione (phthalic anhydride) — fused bicyclic retained name
    ("O=C1OC(=O)c2ccccc21",           "isobenzofuran-1,3-dione"),
    # furan-2,5-dione (maleic anhydride) — unsaturated cyclic anhydride retained name
    ("O=C1C=CC(=O)O1",                "furan-2,5-dione"),
    # regression: acetic anhydride unchanged
    ("CC(=O)OC(=O)C",                 "acetic anhydride"),
    # regression: propionic anhydride unchanged
    ("CCC(=O)OC(=O)CC",               "propanoic anhydride"),
    # regression: succinic anhydride (saturated, non-aromatic ring)
    ("O=C1CCC(=O)O1",                 "oxolane-2,5-dione"),
])
def test_phase264_benzoic_anhydride(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
