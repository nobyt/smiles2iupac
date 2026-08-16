"""Phase 929: 6 dedicated "dual carbonyl/imine-anchor" namer functions
all shared the same chain-selection bug -- found via a fresh systematic
sweep testing geminal (2-instance) constructions of every remaining
`_multi_map` group type, continuing the methodology from the 921-928
arc but targeting a different bug shape.

Root cause: `_name_diimine`, `_name_diacid_halide`, `_name_dithioamide`,
`_name_diselenoamide`, `_name_diester`, and `_name_dicarboxylate` each
independently built their acyclic (non-ring) chain by picking whichever
of two candidate chains -- one grown from anchor C1 via
`_collect_acid_chain`, one grown from anchor C2 -- happened to be
LONGER, with no requirement that either candidate chain actually reach
the OTHER anchor. When the central carbon bears an unrelated substituent
LONGER than the path to the other anchor (e.g. a hexyl branch), the
"longest arm" walk wanders off into that substituent instead, producing
a chain containing only ONE of the two required anchor atoms:

    CCCCCCC(C(=O)Cl)(C(=O)Cl)C -> wrongly "2-methyloctanedioyl dichloride"
        (the chain absorbed the hexyl branch into an 8-carbon "octane"
        stem, treating the SECOND acid chloride as if it weren't there,
        rather than correctly stopping the chain at both -C(=O)Cl carbons
        and reporting hexyl as an ordinary "2-hexyl" substituent)

`_name_diimine` additionally hit a locant-computation crash-adjacent bug
from the same root: when the chosen (wrong) chain didn't contain the
second imine carbon at all, `.get(c, 0)` silently defaulted to locant 0,
producing invalid names like "butane-0,1-diimine" (confirmed unparseable
via OPSIN). `_name_diimine` also never collected chain substituents at
all, so even a CORRECTLY-anchored chain with a branching substituent
(e.g. both ethyl and methyl on the central carbon) lost them entirely.

By contrast, `_name_dioic_acid`/`_name_diamide`/`_name_dinitrile`/
`_name_dial` (structurally identical "dual anchor" shape) were NOT
affected, because those 4 dedicated functions are ring-only handlers
that return `None` for acyclic molecules, falling through to the
already-correct generic `_name_acyclic` pipeline (which benefits from
Phase921's chain_finder.py overlap-maximizing fix). This asymmetry --
same functional shape, different code path, only some affected -- is
exactly why a systematic sweep across the group-type table was needed
to find this rather than a single spot-check.

Fixed by extracting a shared `_find_dual_anchor_chain(graph, c1, c2,
excluded, get_atom)` helper: BFS the actual carbon-only path connecting
c1 and c2 first, THEN extend each end outward via `_collect_acid_chain`
to find the true longest chain that is GUARANTEED to contain both
anchors. All 6 functions now use this helper instead of independently
re-implementing (and each subtly getting wrong) the same "pick the
longer of two candidate chains" logic. `_name_diimine` additionally
gained a `collect_substituents` call it never had.

All fixed names round-trip verified via OPSIN + RDKit canonical-SMILES
matching. Full suite (13542 tests) run clean in isolation before this
test file was written (matches the "isolate a source-only regression
signal" discipline used for every core/widely-shared change in this
session, even though each individual function here is narrower than
Phase927/928's shared substituent scanner).
"""

from smiles2iupac import smiles_to_iupac


class TestDiimineDualAnchorFix:
    def test_branching_diimine_previously_invalid_locant(self):
        assert (
            smiles_to_iupac("CC(C=N)(C=N)CC")
            == "2-ethyl-2-methylpropane-1,3-diimine"
        )

    def test_branching_diimine_with_halogen_substituent(self):
        assert (
            smiles_to_iupac("ClCC(C=N)(C=N)CC")
            == "2-chloromethyl-2-ethylpropane-1,3-diimine"
        )

    def test_plain_diimine_regression(self):
        assert smiles_to_iupac("N=CCC=N") == "propane-1,3-diimine"

    def test_longer_diimine_regression(self):
        assert smiles_to_iupac("N=CCCC=N") == "butane-1,4-diimine"

    def test_carbodiimide_unaffected(self):
        assert smiles_to_iupac("CN=C=NC") == "N,N'-dimethylmethanediimine"


class TestDiacidHalideDualAnchorFix:
    def test_branching_diacid_halide(self):
        assert (
            smiles_to_iupac("CCCCCCC(C(=O)Cl)(C(=O)Cl)C")
            == "2-hexyl-2-methylpropanedioyl dichloride"
        )

    def test_plain_diacid_halide_regression(self):
        assert smiles_to_iupac("ClC(=O)CC(=O)Cl") == "propanedioyl dichloride"

    def test_benzene_diacid_halide_regression(self):
        assert (
            smiles_to_iupac("O=C(Cl)c1ccccc1C(=O)Cl")
            == "benzene-1,2-dicarbonyl dichloride"
        )


class TestDithioamideDualAnchorFix:
    def test_branching_dithioamide(self):
        assert (
            smiles_to_iupac("CCCCCCC(C(=S)N)(C(=S)N)C")
            == "2-hexyl-2-methylpropanedithioamide"
        )

    def test_plain_dithioamide_regression(self):
        assert smiles_to_iupac("NC(=S)CC(=S)N") == "propanedithioamide"


class TestDiesterDualAnchorFix:
    def test_branching_diester(self):
        assert (
            smiles_to_iupac("CCCCCCC(C(=O)OC)(C(=O)OC)C")
            == "dimethyl 2-hexyl-2-methylmalonate"
        )

    def test_plain_diester_regression(self):
        assert smiles_to_iupac("CCOC(=O)CC(=O)OCC") == "diethyl malonate"

    def test_oxalate_diester_regression(self):
        assert smiles_to_iupac("COC(=O)C(=O)OC") == "dimethyl oxalate"


class TestDicarboxylateDualAnchorFix:
    def test_branching_dicarboxylate(self):
        assert (
            smiles_to_iupac("CCCCCCC(C(=O)[O-])(C(=O)[O-])C")
            == "2-hexyl-2-methylmalonate"
        )

    def test_plain_dicarboxylate_regression(self):
        assert smiles_to_iupac("[O-]C(=O)CC(=O)[O-]") == "malonate"

    def test_oxalate_dicarboxylate_regression(self):
        assert smiles_to_iupac("[O-]C(=O)C(=O)[O-]") == "oxalate"
