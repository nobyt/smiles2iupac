"""Phase 867: thio/seleno/telluro carbamic acids (+ telluramide typo fix).

The chalcogen analogues of carbamic acid (H2N-C(=O)-OH) were mis-named:
NC(=S)O -> "methanethioamide", NC(=O)S -> "methanethioic S-acid",
NC(=[Se])O -> "methaneselenoamide". Only the exact unsubstituted NC(=S)S was
handled (a hardcoded lookup), not the general or N-substituted forms.

General detection now covers RnN-C(=X)(YH) with a single heavy chalcogen
element. As with the phosphorus/sulfur acids (Phase 864/865) the tautomers
(H2N-C(=S)-OH <-> H2N-C(=O)-SH) share one name, so NO positional tag is used;
the S/Se/Te count gives thio/dithio, the element gives thio/seleno/telluro:

  NC(=S)O / NC(=O)S   -> carbamothioic acid
  NC(=S)S             -> carbamodithioic acid
  NC(=[Se])O          -> carbamoselenoic acid
  NC(=[Se])[SeH]      -> carbamodiselenoic acid
  NC(=[Te])O          -> carbamotelluroic acid
  CNC(=S)O            -> N-methylcarbamothioic acid

Also fixes a pre-existing typo: telluramide was spelled "teluramide".
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # carbamothioic (tautomers share one name)
    ("NC(=S)O",       "carbamothioic acid"),
    ("NC(=O)S",       "carbamothioic acid"),
    ("NC(=S)S",       "carbamodithioic acid"),
    # carbamoselenoic
    ("NC(=[Se])O",    "carbamoselenoic acid"),
    ("NC(=O)[SeH]",   "carbamoselenoic acid"),
    ("NC(=[Se])[SeH]", "carbamodiselenoic acid"),
    # carbamotelluroic
    ("NC(=[Te])O",    "carbamotelluroic acid"),
    # N-substituted
    ("CNC(=S)O",      "N-methylcarbamothioic acid"),
    ("CNC(=[Se])O",   "N-methylcarbamoselenoic acid"),
    ("CN(C)C(=S)O",   "N,N-dimethylcarbamothioic acid"),
    # regression: plain carbamic acid unchanged
    ("NC(=O)O",       "carbamic acid"),
    ("CNC(=O)O",      "N-methylcarbamic acid"),
    # telluramide typo fix (was "teluramide")
    ("CC(=[Te])N",    "ethanetelluramide"),
    ("[Te]=CN",       "methanetelluramide"),
    # regression: one-N amides and thiourea not swept up
    ("CC(=S)N",       "ethanethioamide"),
    ("CC(=[Se])N",    "ethaneselenoamide"),
    ("NC(=S)N",       "thiourea"),
])
def test_phase867_carbamo_chalcogenoic_acids(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
