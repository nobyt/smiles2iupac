"""Phase 917: esters on a heteroaromatic ring (pyridine, thiophene, furan,
etc.) had two compounding bugs in _name_ester's heteroaromatic branch.

Bug 1 -- wrong ring-numbering DIRECTION, giving an unnecessarily high
locant to the ester's own attachment point:

    c1ccc(nc1)C(=O)OC -> wrongly "methyl pyridine-6-carboxylate"
        (should be "methyl pyridine-2-carboxylate" -- methyl picolinate,
        a common, real compound; "6" and "2" name the IDENTICAL molecule,
        confirmed via RDKit canonical-SMILES round-trip through OPSIN, but
        IUPAC requires the LOWER locant)

Root cause: this branch called heterocycle_handler._find_best_start to
pick the ring's numbering direction, but that function only tie-breaks on
the ring's OWN heteroatom positions / any DIRECTLY-ring-attached C=O --
it has no awareness of an exocyclic principal-group substituent one bond
further out (the ester's own carbonyl carbon). The analogous carboxylic-
-acid code path in name_heterocycle already had the correct fix (flip the
rotation if the reversed locant is lower) via its _EXOCYCLIC_SUFFIX
branch; esters aren't in that dict and never reached the fix. Applied the
identical flip-check directly in _name_ester's heteroaromatic branch.

Bug 2 -- ANY other ring substituent (halogen, alkyl, hydroxy, amino, ...)
was silently dropped entirely, independent of bug 1:

    COC(=O)c1ccc(Cl)nc1 -> wrongly "methyl pyridine-3-carboxylate"
        (should be "methyl 6-chloropyridine-3-carboxylate" -- the
        chlorine vanishes without a trace)

Root cause: unlike the sibling pure-benzene branch just above it (which
correctly collects other ring substituents via name_substituent +
_build_prefix), this heteroaromatic branch built ONLY the bare
"{ring}-{locant}-carboxylate" string and never looked at any other ring
atom's neighbors at all. Fixed by collecting substituents via
heterocycle_handler._collect_hetero_substituents (excluding the ester's
own carbonyl carbon) and splicing the resulting prefix in, mirroring the
pure-benzene branch's own pattern.

All fixed names round-trip verified via OPSIN + RDKit canonical-SMILES
matching.
"""

from smiles2iupac import smiles_to_iupac


class TestHeteroaromaticEsterLocantDirection:
    def test_methyl_picolinate_lowest_locant(self):
        assert smiles_to_iupac("c1ccc(nc1)C(=O)OC") == "methyl pyridine-2-carboxylate"

    def test_methyl_picolinate_alt_smiles(self):
        assert smiles_to_iupac("COC(=O)c1ccccn1") == "methyl pyridine-2-carboxylate"

    def test_methyl_nicotinate_unaffected(self):
        # Position 3 is already symmetric either direction (locant 3 both
        # ways for a monosubstituted 6-ring numbered from N=1), so this
        # case was never broken -- included as a direction-logic regression.
        assert smiles_to_iupac("COC(=O)c1cccnc1") == "methyl pyridine-3-carboxylate"


class TestHeteroaromaticEsterSubstituentsPreserved:
    def test_chloropyridine_ester(self):
        assert (
            smiles_to_iupac("COC(=O)c1ccc(Cl)nc1")
            == "methyl 6-chloropyridine-3-carboxylate"
        )

    def test_methylpyridine_ester(self):
        assert (
            smiles_to_iupac("COC(=O)c1ccc(C)nc1")
            == "methyl 6-methylpyridine-3-carboxylate"
        )

    def test_hydroxypyridine_ester(self):
        assert (
            smiles_to_iupac("Oc1ccc(C(=O)OC)nc1")
            == "methyl 5-hydroxypyridine-2-carboxylate"
        )

    def test_aminopyridine_ester(self):
        assert (
            smiles_to_iupac("COC(=O)c1ccc(N)nc1")
            == "methyl 6-aminopyridine-3-carboxylate"
        )


class TestHeteroaromaticEsterOtherRingTypes:
    def test_thiophene_ester(self):
        assert smiles_to_iupac("COC(=O)c1cccs1") == "methyl thiophene-2-carboxylate"

    def test_bromothiophene_ester(self):
        assert (
            smiles_to_iupac("COC(=O)c1ccc(Br)s1")
            == "methyl 5-bromothiophene-2-carboxylate"
        )

    def test_furan_ester(self):
        assert smiles_to_iupac("COC(=O)c1ccoc1") == "methyl furan-3-carboxylate"


class TestBenzeneEsterRegression:
    def test_plain_methyl_benzoate(self):
        assert smiles_to_iupac("c1ccc(cc1)C(=O)OC") == "methyl benzoate"

    def test_chlorobenzoate_still_correct(self):
        assert smiles_to_iupac("COC(=O)c1ccccc1Cl") == "methyl 2-chlorobenzoate"
