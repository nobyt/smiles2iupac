"""Phase 918: _aryl_sulfonyl_prefix's heteroaromatic branch -- a SHARED
helper used by many callers (thioamide, sulfonic/sulfinic/sulfenic acid
family, diazonium, and more) to build a ring-name-plus-locant prefix for
an exocyclic principal group attached to a heteroaromatic ring -- had the
exact same "silently drops every other ring substituent" bug found in
Phase 917's separate, only-loosely-parallel _name_ester implementation.

    NC(=S)c1ccc(Cl)nc1 -> wrongly "pyridine-3-carbothioamide"
        (should be "6-chloropyridine-3-carbothioamide" -- the chlorine
        vanishes without a trace)

Unlike the sibling pure-benzene branch a few lines above it in the same
function (which correctly collects ring substituents via
collect_ring_substituents), the heteroaromatic branch returned a bare
"{ring}-{locant}-" string and never inspected any other ring atom's
neighbors. Note this branch's locant-DIRECTION logic (picking whichever
rotation gives the lower locant to the exocyclic attachment point) was
already correct -- only the substituent-collection was missing, unlike
Phase 917's _name_ester bug which had both problems.

Fixed by collecting substituents via _collect_hetero_substituents
(excluding the exocyclic attachment atom) and prepending the resulting
prefix, mirroring the adjacent pure-benzene branch's own pattern.

Because _aryl_sulfonyl_prefix is shared, this single fix simultaneously
corrects thioamide, sulfonic acid, and diazonium naming for substituted
heteroaromatic rings (and any other caller reaching this code path) --
verified all three below.

All fixed names round-trip verified via OPSIN + RDKit canonical-SMILES
matching.
"""

from smiles2iupac import smiles_to_iupac


class TestArylSulfonylPrefixHeteroaromaticSubstituents:
    def test_thioamide_chloropyridine(self):
        assert (
            smiles_to_iupac("NC(=S)c1ccc(Cl)nc1")
            == "6-chloropyridine-3-carbothioamide"
        )

    def test_sulfonic_acid_chloropyridine(self):
        assert (
            smiles_to_iupac("OS(=O)(=O)c1ccc(Cl)nc1")
            == "6-chloropyridine-3-sulfonic acid"
        )

    def test_diazonium_chloropyridine(self):
        assert (
            smiles_to_iupac("[N+](#N)c1ccc(Cl)nc1")
            == "6-chloropyridine-3-diazonium"
        )


class TestArylSulfonylPrefixRegressions:
    def test_thioamide_plain_pyridine(self):
        assert smiles_to_iupac("NC(=S)c1ccccn1") == "pyridine-2-carbothioamide"

    def test_thioamide_pyridine_4_position(self):
        assert smiles_to_iupac("NC(=S)c1ccncc1") == "pyridine-4-carbothioamide"

    def test_sulfonic_acid_plain_pyridine(self):
        assert smiles_to_iupac("c1ccc(nc1)S(=O)(=O)O") == "pyridine-2-sulfonic acid"

    def test_diazonium_plain_pyridine(self):
        assert smiles_to_iupac("c1ccc(nc1)[N+]#N") == "pyridine-2-diazonium"

    def test_benzothioamide_still_uses_separate_path(self):
        # Pure benzene rings bypass _aryl_sulfonyl_prefix's heteroaromatic
        # branch entirely (handled by a different, already-correct path).
        assert smiles_to_iupac("NC(=S)c1ccccc1") == "benzothioamide"

    def test_benzenesulfonic_acid_with_substituent(self):
        assert (
            smiles_to_iupac("Cc1ccc(cc1)S(=O)(=O)O") == "4-methylbenzenesulfonic acid"
        )
