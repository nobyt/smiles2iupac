"""Phase 922: extended the Phase 921 branching-polyol chain-fit fix to the
KETONE family (dione/trione/tetraone/pentaone), which has the identical
merge-by-count-only bug.

    CCC(C(C)=O)(C(C)=O)C(C)=O -> wrongly "3-ethylpentane-2,4-trione"
        ("trione" implies 3 locants, but only "2,4" (2 locants) given --
        OPSIN rejects this: "Mismatch between locant and multiplier counts
        (2 and 3)". Correct: "3-acetyl-3-ethylpentane-2,4-dione")

Same root cause as Phase 921: functional_group.py's _multi_map merges
("ketone", 3): "trione" purely by count, with no check that a single
chain can hold all 3 ketone carbons. On a branching (neopentyl-type)
center with 3+ acyl arms, at most 2 can be chain continuations (a path
uses at most 2 of a branching atom's bonds), so the merged "trione"
suffix and the chain's actual 2 on-chain ketone locants disagreed --
name_assembler.py's hardcoded trione locant fallback then papered over
the mismatch with a locant COUNT that didn't match the suffix word,
producing an invalid name.

Fixed by generalizing the Phase921 __init__.py downgrade block (which
recounts on-chain instances after chain_finder's overlap-maximizing chain
selection, then rebuilds pgrp with the correct smaller merge type) from a
`_polyol_multi_types` dict covering only the alcohol family into a
`_polyol_families` dict also covering dione/trione/tetraone/pentaone.
Ketone doesn't have alcohol's geminal-diol atom-dedup complication (a
carbon can't bear two separate ketone C=O groups), so the existing
O-walk-to-parent-C counting technique carries over unchanged and correctly
counts each ketone instance once.

All fixed names round-trip verified via OPSIN + RDKit canonical-SMILES
matching. Full suite (13468 tests) re-run clean after this change with
zero regressions.
"""

from smiles2iupac import smiles_to_iupac


class TestBranchingKetoneDowngrade:
    def test_trione_on_quaternary_center_downgrades_to_dione(self):
        assert (
            smiles_to_iupac("CCC(C(C)=O)(C(C)=O)C(C)=O")
            == "3-acetyl-3-ethylpentane-2,4-dione"
        )

    def test_trione_with_methyl_branch(self):
        assert (
            smiles_to_iupac("CC(=O)C(C)(C(C)=O)C(C)=O")
            == "3-acetyl-3-methylpentane-2,4-dione"
        )

    def test_tetraone_on_fully_substituted_center_downgrades_to_dione(self):
        assert (
            smiles_to_iupac("CC(=O)C(C(C)=O)(C(C)=O)C(C)=O")
            == "3,3-diacetylpentane-2,4-dione"
        )


class TestBranchingKetoneRegressions:
    def test_plain_ketone(self):
        assert smiles_to_iupac("CC(C)=O") == "propan-2-one"

    def test_unbranched_dione(self):
        assert smiles_to_iupac("CC(=O)CC(C)=O") == "pentane-2,4-dione"
