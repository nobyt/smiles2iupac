"""Phase 923: three related substituent-naming bugs found via the same
branching-substituent sweep that led to Phase 921/922, all in
substituent.py's generic sulfur-substituent branch and the shared
enclosing-marks logic in name_assembler.py.

1. -SO3H / -SO2Cl / -SO2NH2 (and the sulfinyl analogues -SO2H / -SOCl /
   -SONH2) as a NON-PRINCIPAL substituent all silently collapsed to a
   bare "sulfonyl"/"sulfinyl", dropping the halogen/hydroxyl/amino
   entirely:

       OC(=O)c1ccc(cc1)S(=O)(=O)Cl -> wrongly "4-sulfonylbenzoic acid"
           (chlorine dropped entirely -- should be
           "4-(chlorosulfonyl)benzoic acid")
       OC(=O)c1ccc(cc1)S(=O)(=O)N -> wrongly "4-sulfonylbenzoic acid"
           (amino dropped entirely -- should be "4-sulfamoylbenzoic acid")
       OC(=O)c1ccc(cc1)S(=O)(=O)O -> wrongly "4-sulfonylbenzoic acid"
           (should be "4-sulfobenzoic acid" -- "sulfonyl" additionally
           implies a DIFFERENT, unrelated functional group, -SO2-)

   Root cause: substituent.py's S-atom branch only ever checked for an
   OTHER CARBON neighbor (the c_other list, for genuine alkyl/aryl-
   sulfonyl substituents like "methylsulfonyl") before falling back to a
   hardcoded "sulfonyl"/"sulfinyl" whenever that list was empty -- never
   inspecting what ELSE (halogen, N, O-H) was actually attached to the S.
   Fixed by checking, in order, for a halogen (-> "{halo}sulfonyl" /
   "{halo}sulfinyl"), an unsubstituted N (-> "sulfamoyl" / "sulfinamoyl"),
   and an O-H (-> "sulfo" / "sulfino"), falling back to the old bare
   "sulfonyl"/"sulfinyl" only when none of those match (e.g. the rare
   N-substituted-sulfonamide-as-a-bare-substituent case, left as a known
   gap since it needs full N-substituent-prefix construction).

2. Acyl HALIDES as a non-principal substituent (-C(=O)Cl etc.) fell into
   the SAME "no other carbon neighbor" fallback used for plain -CHO,
   getting misnamed "formyl" and silently dropping the halogen:

       CCC(C(=O)Cl)(C(=O)Cl)C(=O)Cl -> wrongly involved "diformyl..."
           substituents (chlorine dropped from 2 of the 3 acyl groups)

   Fixed by checking for a halogen neighbor on the acyl root BEFORE the
   "single-carbon acyl = formyl" shortcut in substituent.py's acyl-as-
   substituent branch (mirroring the Phase916 fix that added the O/N/S
   heteroatom check ahead of the same shortcut) -> "{halo}carbonyl".

3. Composite "X-carbonyl" substituent names (chlorocarbonyl,
   methoxycarbonyl) were never flagged as needing bis/tris + enclosing
   marks when there were 2+ of them, so "dichlorocarbonyl" /
   "dimethoxycarbonyl" was produced -- genuinely ambiguous with "(di-
   chloro)carbonyl" (a carbonyl carbon bearing 2 chlorines, a different
   structure). Fixed with a new `_needs_bis_tris_multiplier()` helper in
   name_assembler.py: a superset of `_needs_bis_tris()` that additionally
   flags any name ending in "carbonyl", used ONLY at the n>1 (bis/tris
   vs. plain di/tri) decision points -- deliberately NOT folded into the
   base `_needs_bis_tris()`, since a single (non-multiplied)
   "methoxycarbonyl" is unambiguous and Phase916 already established
   (and this phase preserves) that it must NOT get enclosing parens at
   n=1 ("4-methoxycarbonylbenzoic acid", not
   "4-(methoxycarbonyl)benzoic acid").

All fixed names round-trip verified via OPSIN + RDKit canonical-SMILES
matching. An intermediate version of fix #3 folded the "carbonyl" check
directly into `_needs_bis_tris()`, which broke 4 pre-existing Phase916
tests by adding spurious parens at n=1 -- caught by the full suite before
any commit, and fixed by splitting into the narrower
`_needs_bis_tris_multiplier()` used only where multiplication is actually
happening.
"""

from smiles2iupac import smiles_to_iupac


class TestSulfonylFamilySubstituentDropping:
    def test_sulfonic_acid_as_substituent(self):
        assert (
            smiles_to_iupac("OC(=O)c1ccc(cc1)S(=O)(=O)O")
            == "4-sulfobenzoic acid"
        )

    def test_sulfonyl_chloride_as_substituent(self):
        assert (
            smiles_to_iupac("OC(=O)c1ccc(cc1)S(=O)(=O)Cl")
            == "4-(chlorosulfonyl)benzoic acid"
        )

    def test_sulfonamide_as_substituent(self):
        assert (
            smiles_to_iupac("OC(=O)c1ccc(cc1)S(=O)(=O)N")
            == "4-sulfamoylbenzoic acid"
        )

    def test_methylsulfonyl_regression(self):
        assert (
            smiles_to_iupac("OC(=O)c1ccc(cc1)S(=O)(=O)C")
            == "4-(methylsulfonyl)benzoic acid"
        )

    def test_multiple_sulfo_substituents(self):
        assert (
            smiles_to_iupac("OC(=O)C(S(=O)(=O)O)(S(=O)(=O)O)CC")
            == "2,2-disulfobutanoic acid"
        )

    def test_multiple_chlorosulfonyl_substituents(self):
        assert (
            smiles_to_iupac("OC(=O)C(S(=O)(=O)Cl)(S(=O)(=O)Cl)CC")
            == "2,2-bis(chlorosulfonyl)butanoic acid"
        )


class TestSulfinylFamilySubstituentDropping:
    def test_sulfinic_acid_as_substituent(self):
        assert (
            smiles_to_iupac("OC(=O)c1ccc(cc1)S(=O)O") == "4-sulfinobenzoic acid"
        )

    def test_sulfinyl_chloride_as_substituent(self):
        assert (
            smiles_to_iupac("OC(=O)c1ccc(cc1)S(=O)Cl")
            == "4-(chlorosulfinyl)benzoic acid"
        )

    def test_sulfinamide_as_substituent(self):
        assert (
            smiles_to_iupac("OC(=O)c1ccc(cc1)S(=O)N")
            == "4-sulfinamoylbenzoic acid"
        )

    def test_methylsulfinyl_regression(self):
        assert (
            smiles_to_iupac("OC(=O)c1ccc(cc1)S(=O)C")
            == "4-(methylsulfinyl)benzoic acid"
        )

    def test_thiol_thioether_regression(self):
        assert smiles_to_iupac("OC(=O)c1ccc(cc1)S") == "4-sulfanylbenzoic acid"
        assert (
            smiles_to_iupac("OC(=O)c1ccc(cc1)SC")
            == "4-(methylsulfanyl)benzoic acid"
        )


class TestAcylHalideAsSubstituent:
    def test_triple_acid_chloride_branching(self):
        assert (
            smiles_to_iupac("CCC(C(=O)Cl)(C(=O)Cl)C(=O)Cl")
            == "2,2-bis(chlorocarbonyl)butanoyl chloride"
        )

    def test_triple_acid_fluoride_branching(self):
        assert (
            smiles_to_iupac("CCC(C(=O)F)(C(=O)F)C(=O)F")
            == "2,2-bis(fluorocarbonyl)butanoyl fluoride"
        )

    def test_plain_acetyl_chloride_regression(self):
        assert smiles_to_iupac("CC(=O)Cl") == "acetyl chloride"

    def test_benzoyl_chloride_regression(self):
        assert smiles_to_iupac("ClC(=O)c1ccccc1") == "benzoyl chloride"


class TestCarbonylCompositeBisTris:
    def test_single_methoxycarbonyl_no_parens(self):
        # Phase916, must remain unchanged: single instance, no bis/tris.
        assert (
            smiles_to_iupac("OC(=O)c1ccc(C(=O)OC)cc1")
            == "4-methoxycarbonylbenzoic acid"
        )

    def test_double_methoxycarbonyl_needs_bis(self):
        assert (
            smiles_to_iupac("CCC(C(=O)OC)(C(=O)OC)C(=O)OC")
            == "methyl 2,2-bis(methoxycarbonyl)butanoate"
        )

    def test_single_chlorocarbonyl_no_parens(self):
        assert smiles_to_iupac("CCC(C(=O)Cl)(CC)") == "2-ethylbutanoyl chloride"
