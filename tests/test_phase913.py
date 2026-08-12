"""Phase 913: a chiral (R/S) stereocenter on a standalone namer's OWN main
chain was silently dropped whenever that chain had no C=C/C#C double bond.

Phase 912 fixed stereo loss inside SUBSTITUENT branches (via the shared
_name_carbon_substituent helper). This phase found a related but distinct
bug: many standalone group_namers.py functions (sulfonic/sulfenic acid,
sulfonyl/sulfenyl halides, carboxylate/dicarboxylate anions, anhydrides,
thioamide, amidine, and ~20 more) already HAD their own local
assign_stereochemistry() call for their own chain -- but it was gated
behind `if _ene_xxx:` (an "is there a C=C/C#C in this chain" check), so it
only ever ran when the chain ALSO happened to contain a double bond. A
plain R/S stereocenter with no double bond present silently lost its
descriptor entirely, e.g.:

    C[C@H](Cl)C(=S)N -> wrongly "2-chloropropanethioamide"
        (should be "(2S)-2-chloropropanethioamide")
    C[C@H](Cl)S(=O)(=O)O -> wrongly "1-chloroethanesulfonic acid"
        (should be "(S)-1-chloroethanesulfonic acid")

Found by systematically grepping every assign_stereochemistry(graph, call
site in group_namers.py (29 total) and checking its nearest enclosing `if`
-- ALL 29 were gated this way. Fixed by removing/hoisting the gate so the
R/S computation always runs (E/Z computation for any double bond present
is unaffected, since assign_stereochemistry already computes both from
the same chain path in one call).

A SECOND, compounding bug was found while verifying each of the 29 fixes:
in ~6 functions (_name_dicarboxylate, _name_carboxylate, the anhydride
family's _acid_stem_name, sulfonic/sulfenic acid, sulfonyl chloride), the
stereo prefix WAS now being computed unconditionally, but was still only
actually spliced into the ene/yne-branch's return string -- the plain
(no-double-bond) branch's own return statement never referenced the
prefix variable at all, so the fix was a no-op for exactly the cases it
was meant to catch. Fixed by adding the prefix to every return path.

All fixed names round-trip verified via OPSIN + RDKit canonical-SMILES
matching (including chirality tags).

NOTE: ~25 more group_namers.py functions (see project memory) build their
own chain but never call assign_stereochemistry AT ALL (not even the
ene-gated version) -- e.g. _name_acyl_azide, _name_sulfenamide,
_name_sulfenate_ester, _name_sulfonimidamide. These still silently drop
stereo unconditionally and are a queued item for a future session, same
root cause but a larger, ground-up addition rather than an unblocking fix.
"""

from smiles2iupac import smiles_to_iupac


class TestMainChainStereoNoDoubleBond:
    def test_thioamide(self):
        assert smiles_to_iupac("C[C@H](Cl)C(=S)N") == "(2S)-2-chloropropanethioamide"

    def test_amidine(self):
        assert smiles_to_iupac("C[C@H](Cl)C(=N)N") == "(2S)-2-chloropropanimidamide"

    def test_sulfonic_acid(self):
        assert (
            smiles_to_iupac("C[C@H](Cl)S(=O)(=O)O") == "(S)-1-chloroethanesulfonic acid"
        )

    def test_sulfenic_acid(self):
        assert smiles_to_iupac("C[C@H](Cl)SO") == "(S)-1-chloroethanesulfenic acid"

    def test_sulfonyl_chloride(self):
        assert (
            smiles_to_iupac("C[C@H](Cl)S(=O)(=O)Cl")
            == "(S)-1-chloroethanesulfonyl chloride"
        )

    def test_carboxylate_anion(self):
        assert smiles_to_iupac("C[C@H](Cl)C(=O)[O-]") == "(2S)-2-chloropropanoate"

    def test_dicarboxylate_dianion_retained_name_branch(self):
        # Exercises the retained-name path (adipate), which is NOT the
        # ene/yne branch -- this specifically regression-tests the "computed
        # but not spliced into the plain-branch return" compounding bug.
        assert (
            smiles_to_iupac("[O-]C(=O)C[C@H](Cl)CCC(=O)[O-]") == "(3R)-3-chloroadipate"
        )

    def test_anhydride(self):
        assert (
            smiles_to_iupac("C[C@H](Cl)C(=O)OC(=O)C")
            == "(2S)-2-chloropropanoic acetic anhydride"
        )


class TestMainChainStereoRegressionWithDoubleBond:
    # These already worked before Phase 913 (the ene-gated path) -- must
    # stay correct now that the gate is gone.
    def test_carboxylic_acid_with_ene_and_stereocenter(self):
        assert (
            smiles_to_iupac("C/C=C/[C@H](Cl)C(=O)O")
            == "(2S,3E)-2-chloropent-3-enoic acid"
        )

    def test_achiral_regression(self):
        assert smiles_to_iupac("CC(Cl)C(=S)N") == "2-chloropropanethioamide"
        assert smiles_to_iupac("CS(=O)(=O)O") == "methanesulfonic acid"

    def test_symmetric_dicarboxylate_correctly_shows_no_stereo(self):
        # C2 here has two IDENTICAL -C(=O)[O-] arms, so despite the SMILES
        # @ marker it is genuinely not a stereocenter (RDKit assigns no CIP
        # code) -- confirms the fix doesn't fabricate spurious descriptors.
        assert smiles_to_iupac("[O-]C(=O)[C@H](Cl)C(=O)[O-]") == "2-chloromalonate"
