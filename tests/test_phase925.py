"""Phase 925: two related, deep pre-existing bugs in the amine-naming
pipeline, found while investigating the mixed primary/tertiary-amine
merge case deliberately deferred at the end of Phase 924.

1. `aggregate_groups()` merged secondary/tertiary amine instances into
   diamine/triamine purely by instance count, exactly like the primary-
   only cases Phase921-924 fixed -- but for secondary/tertiary amines,
   `atom_indices` is `[n_idx] + c_neighbors`, which includes ALL of that
   nitrogen's carbon substituents (its OWN N-alkyl groups), not just a
   single "chain-worthy" carbon. Once merged, the N atom itself becomes
   part of the principal-group exclusion set passed to
   `collect_substituents`, which makes that N (and everything hanging
   off it, including its N-alkyl substituents) completely UNREACHABLE
   as an ordinary substituent -- it silently vanishes:

       CNCCNC -> wrongly "ethane-1,2-diamine"
           (BOTH N-methyl groups silently dropped -- this is the name of
           a different, simpler molecule, H2N-CH2CH2-NH2, confirmed
           structurally distinct via RDKit canonical SMILES)

   This is worse than an invalid name: it's a valid name for the WRONG
   molecule. Unlike the Phase921/922/924 branching/geminal cases (where
   the "extra" instance is always a plain terminal -OH/-CH2NH2/=O that
   collect_substituents can trivially re-discover once released), a
   secondary/tertiary amine's OWN alkyl substituents need proper
   N-prefix construction, which the merged-suffix pathway has no
   mechanism for at all.

   Fixed at the true root: `aggregate_groups()` (functional_group.py)
   now skips the diamine/triamine merge entirely whenever any instance
   being merged is not a plain primary amine (`len(atom_indices) != 2`,
   mirroring the existing amide shared-N merge guard in the same
   function). This routes such molecules through the pre-existing,
   already-correct single-amine-plus-N-substituent-prefix pathway
   (`_name_secondary_tertiary_amine`) instead of ever reaching the
   broken merged-suffix path.

   CNCCNC -> "N-methyl-2-methylaminoethanamine" (verified round-trip)

2. This fix, by routing more molecules through the amine-as-substituent
   naming code (`_name_nitrogen_substituent`), exposed a SEPARATE,
   independently pre-existing bug dating back to Phase 520: an
   asymmetric N,N-disubstituted amino substituent (different R and R',
   e.g. -N(CH3)(CH2CH3)) was named with a COMMA-JOINED format inside one
   pair of parens, "(ethyl,methyl)amino" -- established and tested since
   Phase520, but never actually verified to round-trip correctly.
   Confirmed via OPSIN + RDKit that this comma format does NOT parse
   back to the right structure: OPSIN silently merges the two comma-
   separated names into a single different alkyl group ("propyl")
   instead of recognizing two independent N-substituents.

       CN(CC)CC(=O)O -> wrongly "2-[(ethyl,methyl)amino]acetic acid"
           (OPSIN parses this to C(C)CNCC(=O)O, i.e. N-propylglycine --
           a DIFFERENT molecule, confirmed via RDKit)

   Fixed by individually bracketing EACH substituent name and
   concatenating them with NO separator: "(ethyl)(methyl)amino".
   Verified via OPSIN + RDKit to round-trip correctly. Always bracketing
   both parts (not just ones containing digits/hyphens, which was the
   original asymmetric-case rule) avoids any ambiguity with a bare
   concatenation like "ethylmethylamino" that could misread as one
   composite substituent name.

Two pre-existing tests had stale (confirmed-wrong) expectations updated
in this phase: test_phase520.py's asymmetric-amino case, and
test_phase63.py's mixed-amine "triamine" case (which was already
OPSIN-unparseable before this session -- a purely latent, pre-existing
bug, unrelated to anything touched in 921-924, just newly surfaced as a
test failure once this phase's fix changed the specific wrong output).
test_phase924.py's own test asserting the OLD broken value as
"deliberately unfixed" was also updated now that it IS fixed.
"""

from smiles2iupac import smiles_to_iupac


class TestSecondaryAmineNoLongerMerges:
    def test_two_secondary_amines_both_methyls_preserved(self):
        assert (
            smiles_to_iupac("CNCCNC") == "N-methyl-2-methylaminoethanamine"
        )

    def test_two_secondary_amines_branched_chain(self):
        assert (
            smiles_to_iupac("CNCC(C)NC")
            == "N-methyl-2-methylaminopropan-1-amine"
        )

    def test_tertiary_plus_two_primary_amines(self):
        assert (
            smiles_to_iupac("NCCCN(CCCN)CCC")
            == "3-[(3-aminopropyl)(propyl)amino]propan-1-amine"
        )


class TestPureAmineMergeStillWorks:
    def test_two_primary_amines_still_merge_to_diamine(self):
        assert smiles_to_iupac("NCCN") == "ethane-1,2-diamine"

    def test_geminal_triamine_still_fixed(self):
        assert smiles_to_iupac("CCC(N)(N)N") == "propane-1,1,1-triamine"

    def test_branching_triamine_still_fixed(self):
        assert (
            smiles_to_iupac("CCC(CN)(CN)CN")
            == "2-(aminomethyl)-2-ethylpropane-1,3-diamine"
        )


class TestAsymmetricDisubstitutedAmino:
    def test_ethyl_methyl_amino_substituent(self):
        assert (
            smiles_to_iupac("CN(CC)CC(=O)O")
            == "2-[(ethyl)(methyl)amino]acetic acid"
        )

    def test_symmetric_dimethylamino_unaffected(self):
        assert (
            smiles_to_iupac("CN(C)CC(=O)O") == "2-(dimethylamino)acetic acid"
        )

    def test_asymmetric_aniline_substituent(self):
        assert smiles_to_iupac("CN(CC)c1ccccc1") == "N-ethyl-N-methylaniline"


class TestPlainAmineRegressions:
    def test_simple_secondary_amine(self):
        assert smiles_to_iupac("CNCC") == "N-methylethanamine"

    def test_simple_tertiary_amine_same_r(self):
        assert smiles_to_iupac("CN(C)CC") == "N,N-dimethylethanamine"

    def test_n_methylaniline(self):
        assert smiles_to_iupac("CNc1ccccc1") == "N-methylaniline"

    def test_n_n_dimethylaniline(self):
        assert smiles_to_iupac("CN(C)c1ccccc1") == "N,N-dimethylaniline"
