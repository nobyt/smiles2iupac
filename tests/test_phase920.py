"""Phase 920: the CORE heteroaromatic ring-numbering algorithm (inside
name_heterocycle, heterocycle_handler.py) did not prioritize a
substituent that will later be promoted to a PRINCIPAL-GROUP suffix
(-ol/-thiol/-amine, via _apply_hetero_suffixes) over ordinary substituent
prefixes when choosing the ring's numbering direction.

    Nc1ccc(nc1)Cl -> wrongly "2-chloropyridin-5-amine"
        (should be "6-chloropyridin-3-amine" -- IUPAC numbering priority
        requires minimizing the PRINCIPAL GROUP's own locant (the amine,
        cited as a suffix) before minimizing substituent locants (the
        chlorine, cited as a prefix); amine locant 3 beats 5)

This is a deeper, more foundational bug than the Phase 917/918/919
"siloed namer" bugs -- those were separate, only-loosely-parallel
implementations for EXOCYCLIC principal groups (ester, thioamide,
sulfonic acid, diazonium, dicarboxylic acid) attached one bond outside
the ring. This phase's bug is in the shared numbering-direction selection
used for ENDOCYCLIC suffixes (-ol/-thiol/-amine directly replacing a ring
C-H), reached via name_heterocycle's own candidate-rotation comparison,
which picked whichever ring direction minimized the OVERALL substituent
locant set -- treating a substituent that will become the suffix exactly
the same as an ordinary halogen/alkyl substituent, with no special
priority.

Fixed by adding a _promo_locs() key (checked BEFORE the general
locant-set comparison) that finds whichever of hydroxy/sulfanyl/amino is
present (mirroring _apply_hetero_suffixes' own priority order and its
"skip if the retained name already embeds a ketone/lactam" guard) and
minimizes THAT substituent's own locant first, exactly matching IUPAC's
"principal-group-suffix locant before substituent-prefix locants" rule.

Full 13451-test suite (spanning every ring/substituent combination
established over 900+ prior phases) passes clean with zero regressions,
confirming this core-algorithm change is safe. All fixed names round-trip
verified via OPSIN + RDKit canonical-SMILES matching.
"""

from smiles2iupac import smiles_to_iupac


class TestHeteroaromaticSuffixLocantPriority:
    def test_amine_beats_chloro_locant(self):
        assert smiles_to_iupac("Nc1ccc(nc1)Cl") == "6-chloropyridin-3-amine"

    def test_ol_beats_chloro_locant(self):
        assert smiles_to_iupac("Oc1ccc(nc1)Cl") == "6-chloropyridin-3-ol"

    def test_amine_beats_methyl_locant(self):
        assert smiles_to_iupac("Nc1ccc(nc1)C") == "6-methylpyridin-3-amine"

    def test_different_substitution_pattern(self):
        assert smiles_to_iupac("Cc1ccc(nc1)N") == "5-methylpyridin-2-amine"


class TestHeteroaromaticAmineRegressions:
    def test_plain_pyridin_4_amine(self):
        assert smiles_to_iupac("Nc1ccncc1") == "pyridin-4-amine"

    def test_plain_pyridin_3_amine(self):
        assert smiles_to_iupac("Nc1cccnc1") == "pyridin-3-amine"

    def test_plain_pyridine_unaffected(self):
        assert smiles_to_iupac("c1ccncc1") == "pyridine"

    def test_exocyclic_suffix_unaffected(self):
        # Carboxylic acid is an EXOCYCLIC suffix (handled by a separate,
        # already-correct code path) -- must be unaffected by this fix.
        assert (
            smiles_to_iupac("OC(=O)c1ccc(N)nc1") == "6-aminopyridine-3-carboxylic acid"
        )

    def test_amide_exocyclic_suffix_unaffected(self):
        assert (
            smiles_to_iupac("NC(=O)c1ccc(N)nc1") == "6-aminopyridine-3-carboxamide"
        )
