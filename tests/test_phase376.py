"""Phase 376: Silyl ether naming (IUPAC 2013 P-68.4).

R-O-Si bonds were previously ignored, naming the compound as if only the
C-Si bonds existed (e.g., methoxytrimethylsilane → trimethylsilane).
Added silyl_ether_org detection and a naming function that forms alkoxy
prefixes from the O-substituents using the standard -oxy suffix.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Monosilyl ethers
    ("[Si](C)(C)(C)OC",         "methoxytrimethylsilane"),
    ("CO[Si](C)(C)C",           "methoxytrimethylsilane"),
    ("[Si](C)(C)(C)OCC",        "ethoxytrimethylsilane"),
    ("C[Si](C)(C)OCC",          "ethoxytrimethylsilane"),
    ("[Si](C)(C)(C)OCCC",       "trimethylpropoxysilane"),
    # Branched alkoxy
    ("[Si](C)(C)(C)OC(C)(C)C",  "(2-methylpropan-2-oxy)trimethylsilane"),
    # Disilyl ether (two alkoxy groups)
    ("C[Si](C)(OCC)OCC",        "diethoxydimethylsilane"),
    # Regressions: silanol unchanged
    ("[Si](C)(C)(C)O",          "trimethylsilanol"),
    # Regressions: silane unchanged
    ("C[Si](C)C",               "trimethylsilane"),
    ("[Si](C)(C)(C)C",          "tetramethylsilane"),
])
def test_phase376_silyl_ether(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
