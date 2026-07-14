"""Phase 859: primary/secondary phosphine/arsane/bismuthane/stibane oxides & sulfides.

The tertiary cases R3El=X (El = P/As/Bi/Sb, X = O/S) were handled in
Phase 187/857/858, but the primary (RElH2=X) and secondary (R2ElH=X)
analogues were broken:

  - For P, R2PH=O / RPH2=O were captured by the phosphinate/phosphonate
    *ester* elif branches whose bodies only appended when an ester oxygen
    was present, so nothing was appended and the compound fell through to
    the generic pathway producing garbage like "(P)methane".
  - R2PH=S / RPH2=S fell into the plain "phosphane" branch, silently
    dropping the =S entirely ("dimethylphosphane").
  - The As/Bi/Sb oxide & sulfide branches gated on len(c_neighbors) >= 3,
    so n<3 dropped the chalcogen the same way.

Fix (Phase 859):
  - Moved the ester-oxygen requirement into the phosphonate_ester /
    phosphinate_ester elif *conditions* so real phosphine oxides reach
    the phosphine_oxide branch.
  - Loosened phosphine_oxide/phosphine_sulfide and the arsane/bismuthane/
    stibane oxide/sulfide detectors from >= 3 to >= 1 carbon substituents.

Naming is count-generic via _name_by_c_substituents, so the same
"{alkyl(s)}phosphane oxide" pattern extends automatically.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # secondary phosphine oxide / sulfide (R2PH=X)
    ("CP(C)=O",          "dimethylphosphane oxide"),
    ("CP(C)=S",          "dimethylphosphane sulfide"),
    ("CC[PH](CC)=O",     "diethylphosphane oxide"),
    ("CC[PH](CC)=S",     "diethylphosphane sulfide"),
    # primary phosphine oxide / sulfide (RPH2=X)
    ("CP=O",             "methylphosphane oxide"),
    ("CP=S",             "methylphosphane sulfide"),
    # secondary / primary arsane oxide / sulfide
    ("C[AsH](C)=O",      "dimethylarsane oxide"),
    ("C[AsH](C)=S",      "dimethylarsane sulfide"),
    ("C[As]=O",          "methylarsane oxide"),
    ("C[As]=S",          "methylarsane sulfide"),
    # secondary bismuthane / stibane oxide / sulfide
    ("C[Bi](C)=O",       "dimethylbismuthane oxide"),
    ("C[Bi](C)=S",       "dimethylbismuthane sulfide"),
    ("C[Sb](C)=O",       "dimethylstibane oxide"),
    ("C[Sb](C)=S",       "dimethylstibane sulfide"),
    # regression: tertiary unchanged (Phase 187/857/858)
    ("C[P](C)(C)=O",     "trimethylphosphane oxide"),
    ("C[P](C)(C)=S",     "trimethylphosphane sulfide"),
    ("C[As](C)(C)=O",    "trimethylarsane oxide"),
    ("C[As](C)(C)=S",    "trimethylarsane sulfide"),
    ("C[Bi](C)(C)=O",    "trimethylbismuthane oxide"),
    ("C[Sb](C)(C)=O",    "trimethylstibane oxide"),
    # regression: phosphinate / phosphonate esters must NOT be swallowed
    ("CP(=O)(OC)OC",     "dimethyl methylphosphonate"),
    ("CCP(C)(=O)OC",     "methyl ethylmethylphosphinate"),
    # regression: phosphonic / phosphinic acids unchanged
    ("CP(=O)(O)O",       "methylphosphonic acid"),
    ("CP(C)(=O)O",       "dimethylphosphinic acid"),
    # regression: plain phosphanes/arsanes/bismuthanes/stibanes unchanged
    ("CP(C)C",           "trimethylphosphane"),
    ("C[As](C)C",        "trimethylarsane"),
    ("C[Bi](C)(C)C",     "tetramethylbismuthane"),
    ("C[Sb](C)C",        "trimethylstibane"),
])
def test_phase859_primary_secondary_pnictogen_chalcogenides(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
