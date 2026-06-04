"""Phase 398: Pyrazinone and pyridazinone naming (IUPAC 2013 P-31.1.7).

Added two _RETAINED_NAMES entries:
  ("NH","C","C","N","C","C") → ("pyrazine",   True)  — pyrazin-2(1H)-one
  ("NH","C","C","C","C","N") → ("pyridazine",  True)  — pyridazin-3(2H)-one

The Phase 397 multi-N locant algorithm handles both rings correctly.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # pyrazinone
    ("O=C1NC=CN=C1",    "pyrazin-2(1H)-one"),
    # pyridazinone
    ("O=C1C=CC=NN1",    "pyridazin-3(2H)-one"),
    # regression: plain parent rings
    ("c1cnccn1",         "pyrazine"),
    ("c1ccnnc1",         "pyridazine"),
    ("c1ccncn1",         "pyrimidine"),
    ("c1ccccn1",         "pyridine"),
    # regression: Phase 397 pyrimidinones
    ("O=C1NC=CC=N1",    "pyrimidin-2(1H)-one"),
    ("O=C1NC=NC=C1",    "pyrimidin-4(3H)-one"),
    # regression: Phase 396 pyridinone
    ("O=C1NC=CC=C1",    "pyridin-2(1H)-one"),
    # regression: 1H-imidazole (no exo C=O)
    ("c1cnc[nH]1",       "1H-imidazole"),
])
def test_phase398_azinones(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
