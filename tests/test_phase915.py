"""Phase 915: _name_branched_substituent (the recursive branched-alkyl
namer in substituent.py, used for isopropyl/sec-butyl-shaped substituents
whose attachment point is an internal, non-terminal carbon) never checked
chirality either -- the last unfixed gap flagged at the end of the
Phase 912/913/914 stereo investigation.

OCCCC(CCCCC)[C@H](C)CC -> wrongly "4-(butan-2-yl)nonan-1-ol"
    (should be "4-[(2R)-butan-2-yl]nonan-1-ol")

Fixed the same way as Phase 912/914: compute the stereo descriptor once
(assign_stereochemistry against the substituent's own through-path
`main_path` + `main_locant`, already computed by this function for
locant-numbering purposes) and prepend it to all 6 return paths (root at
chain end vs. internal position, each with/without its own further
substituents, each with/without unsaturation).

This closes the entire "912/913/914/915 chiral stereocenter dropped"
investigation -- every substituent-naming and chain-building code path in
the codebase has now been verified to preserve R/S (and E/Z, where
applicable) stereo descriptors.

All fixed names round-trip verified via OPSIN + RDKit canonical-SMILES
matching.
"""

from smiles2iupac import smiles_to_iupac


class TestBranchedSubstituentStereo:
    def test_internal_root_chiral_center_on_alcohol_chain(self):
        assert (
            smiles_to_iupac("OCCCC(CCCCC)[C@H](C)CC")
            == "4-[(2R)-butan-2-yl]nonan-1-ol"
        )

    def test_internal_root_chiral_center_on_benzene(self):
        assert smiles_to_iupac("c1ccccc1[C@H](C)CC") == "[(2R)-butan-2-yl]benzene"

    def test_internal_root_chiral_longer_chain(self):
        assert (
            smiles_to_iupac("OCCCC(CCCCCC)[C@H](C)CC")
            == "4-[(2R)-butan-2-yl]decan-1-ol"
        )


class TestBranchedSubstituentAchiralRegression:
    def test_achiral_branched_substituent_unaffected(self):
        assert smiles_to_iupac("c1ccccc1C(C)CC") == "(butan-2-yl)benzene"
