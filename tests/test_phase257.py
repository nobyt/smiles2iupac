"""Phase 257: hydrazide suffix fix and cyclic imide PIN update (IUPAC 2013).

  CC(=O)NN      → ethanohydrazide    (was ethanehydrazide — suffix -ohydrazide)
  CCC(=O)NN     → propanohydrazide
  C=CC(=O)NN    → prop-2-enohydrazide
  O=C1CCC(=O)N1 → pyrrolidine-2,5-dione  (PIN; succinimide is not a PIN)
  O=C1CCCC(=O)N1→ piperidine-2,6-dione   (PIN; glutarimide is not a PIN)
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # hydrazide suffix: -ohydrazide (IUPAC 2013 P-66.3.3.1)
    ("CC(=O)NN",      "ethanohydrazide"),
    ("CCC(=O)NN",     "propanohydrazide"),
    ("CCCC(=O)NN",    "butanohydrazide"),
    ("CCCCC(=O)NN",   "pentanohydrazide"),
    ("C=CC(=O)NN",    "prop-2-enohydrazide"),
    # N-substituted hydrazides
    ("CC(=O)NNC",     "N'-methylethanohydrazide"),
    ("CC(=O)NNCC",    "N'-ethylethanohydrazide"),
    # cyclic imide PINs (IUPAC 2013 P-66.8.3)
    ("O=C1CCC(=O)N1",  "pyrrolidine-2,5-dione"),
    ("O=C1NC(=O)CC1",  "pyrrolidine-2,5-dione"),  # same compound
    ("O=C1CCCC(=O)N1", "piperidine-2,6-dione"),
    ("O=C1NC(=O)CCC1", "piperidine-2,6-dione"),   # same compound
    # regression: phthalimide retained
    ("O=C1NC(=O)c2ccccc21", "phthalimide"),
    # regression: lactams (single C=O) unchanged
    ("O=C1CCCN1",  "pyrrolidin-2-one"),
    ("O=C1CCCCN1", "piperidin-2-one"),
    # regression: benzohydrazide retained
    ("O=C(NN)c1ccccc1", "benzohydrazide"),
])
def test_phase257_hydrazide_imide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
