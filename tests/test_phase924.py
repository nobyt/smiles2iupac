"""Phase 924: extended the Phase 921/922 branching/geminal chain-fit fix
to the amine family (diamine/triamine), found via the continuing sweep
of the other `_multi_map` group types.

1. Geminal triamine (3 -NH2 on ONE carbon, the amine analogue of
   ethane-1,1-diol):

       CCC(N)(N)N -> wrongly "propane-1,2-triamine"
           (OPSIN rejects it: "Mismatch between locant and multiplier
           counts (2 and 3)" -- correct: "propane-1,1,1-triamine")

   Root cause: unlike diol/triol/dithiol, diamine/triamine had NO
   geminal-instance special case in __init__.py's suffix-locant
   computation, so aggregate_groups()'s atom-dedup (the shared carbon
   only appears once in the merged atom_indices) meant only ONE locant
   was ever collected for N instances that all share the same carbon,
   collapsing to `None` and hitting the hardcoded triamine fallback.
   Fixed by extending the existing diol/triol geminal-walk branch to
   also cover diamine/triamine (walking each N to its parent C, exactly
   like the O-walk used for diol/triol).

2. Branching triamine (3 -CH2NH2 arms on a quaternary carbon, the amine
   analogue of trimethylolpropane):

       CCC(CN)(CN)CN -> wrongly "2-ethylpropane-1,3-triamine"
           (same invalid-name signature: only 2 locants for triamine)

   Fixed by extending the Phase921/922 `_polyol_families` downgrade
   dict (already generalized in 922 to cover the ketone family) to also
   cover diamine/triamine, reusing the same recount-and-rebuild-pgrp
   logic with an added heteroatom-symbol parameter ("O" for
   alcohol/ketone, "N" for amine) so no per-family code duplication was
   needed.

3. Found and avoided a REAL regression during this fix: a molecule with
   THREE amine nitrogens where one is a TERTIARY amine bearing 2 EXTRA
   carbon substituents beyond the "on-chain" one (e.g.
   `NCCCN(CCCN)CCC`, a tris-aminopropylamine-type structure) also merges
   into "triamine" purely by count, but the tertiary nitrogen's extra
   substituents are NOT plain -NH2/-NH- groups that can become ordinary
   substituents the way an off-chain -CH2OH or -CH2NH2 branch can --
   they need to be named as N-alkyl prefixes on the resulting amine
   parent (correct name, confirmed via OPSIN: "N-(3-aminopropyl)-N-
   propylpropane-1,3-diamine"), which this fix's mechanism cannot
   produce. The FIRST version of this fix applied the same downgrade
   logic unconditionally and silently DROPPED the tertiary nitrogen's
   extra branches entirely (worse than the PRE-EXISTING bug, which at
   least produced an invalid-but-honest-looking name that OPSIN would
   reject rather than a valid-but-wrong, silently-incomplete one).
   Caught by the pre-existing full suite (`test_phase63.py`) before any
   commit. Fixed with a purity guard (`_is_pure_primary_amine_merge`):
   the new geminal-walk and downgrade logic only applies when EVERY
   nitrogen in the merged group has exactly one carbon neighbor (i.e.
   every instance is a plain primary amine); mixed primary+secondary/
   tertiary merges fall through to the untouched, pre-existing
   (already broken -- "propane-1,3-triamine" was ALREADY an invalid,
   OPSIN-unparseable name before this session, confirmed via OPSIN)
   generic code path, deliberately left as a known, pre-existing,
   out-of-scope bug for a future session rather than half-fixed here.

All fixed names round-trip verified via OPSIN + RDKit canonical-SMILES
matching. Full suite (13491 tests with this file) run clean twice --
once to confirm the initial (buggy, later reverted) version's source
change didn't regress anything else, and once after the purity guard
was added.
"""

from smiles2iupac import smiles_to_iupac


class TestGeminalTriamine:
    def test_triamine_on_one_carbon(self):
        assert smiles_to_iupac("CCC(N)(N)N") == "propane-1,1,1-triamine"

    def test_diamine_on_one_carbon_regression(self):
        assert smiles_to_iupac("NC(N)CC") == "propane-1,1-diamine"


class TestBranchingTriamine:
    def test_triamine_on_quaternary_center_downgrades_to_diamine(self):
        assert (
            smiles_to_iupac("CCC(CN)(CN)CN")
            == "2-(aminomethyl)-2-ethylpropane-1,3-diamine"
        )

    def test_four_instance_case_does_not_merge_regression(self):
        # No _multi_map entry for ("amine", 4) -- never merges, so this
        # sidesteps the bug entirely (must remain unaffected by this fix).
        assert (
            smiles_to_iupac("NCC(CN)(CN)CN")
            == "3-amino-2,2-bis(aminomethyl)propan-1-amine"
        )


class TestPlainAmineRegressions:
    def test_plain_diamine(self):
        assert smiles_to_iupac("NCCN") == "ethane-1,2-diamine"

    def test_plain_primary_amine(self):
        assert smiles_to_iupac("CC(N)C") == "propan-2-amine"


class TestMixedPrimaryTertiaryAmineMergeUnaffected:
    def test_tris_aminopropylamine_type_unchanged_by_this_fix(self):
        # Pre-existing bug, NOT fixed by this phase (deliberately out of
        # scope -- see module docstring point 3): this name is already
        # invalid (OPSIN: "Mismatch between locant and multiplier counts
        # (2 and 3)") both before and after this phase's changes. The
        # purpose of this test is only to confirm the NEW geminal/
        # branching-downgrade logic added in this phase does not further
        # corrupt this pre-existing wrong output into something worse
        # (e.g. silently dropping the tertiary nitrogen's substituents).
        assert smiles_to_iupac("NCCCN(CCCN)CCC") == "propane-1,3-triamine"
