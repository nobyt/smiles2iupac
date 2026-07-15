"""Phase 866: selenourea and tellurourea (Se/Te analogues of thiourea).

NC(=[Se])N and NC(=[Te])N collided with the seleno/telluramide namers and were
wrongly given "methaneselenoamide" / "methaneteluramide", dropping the second
nitrogen. The codebase already names NC(=S)N "thiourea"; this generalises
_name_thiourea_if_match to the =Se / =Te chalcogens:

  NC(=[Se])N   -> selenourea
  NC(=[Te])N   -> tellurourea
  CNC(=[Se])NC -> N,N'-dimethylselenourea

Selenoamides/telluramides (one nitrogen) are unaffected -- the urea pattern
requires exactly two nitrogens on the C(=X) centre.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # selenourea / tellurourea and substituted variants
    ("NC(=[Se])N",     "selenourea"),
    ("NC(=[Te])N",     "tellurourea"),
    ("CNC(=[Se])N",    "N-methylselenourea"),
    ("CNC(=[Se])NC",   "N,N'-dimethylselenourea"),
    ("CNC(=[Te])N",    "N-methyltellurourea"),
    ("CCNC(=[Se])N",   "N-ethylselenourea"),
    # regression: urea / thiourea unchanged
    ("NC(=O)N",        "urea"),
    ("NC(=S)N",        "thiourea"),
    ("CNC(=S)N",       "N-methylthiourea"),
    ("CNC(=S)NC",      "N,N'-dimethylthiourea"),
    # regression: seleno/thioamides (one N) not swept up as ureas
    ("CC(=[Se])N",     "ethaneselenoamide"),
    ("CC(=S)N",        "ethanethioamide"),
])
def test_phase866_seleno_tellurourea(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
