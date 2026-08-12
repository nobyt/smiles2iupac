"""Phase 912: chiral (R/S) stereocenters inside a SUBSTITUENT branch were
silently dropped everywhere in the codebase.

_name_carbon_substituent (substituent.py) is the single shared low-level
helper used to name every carbon substituent branch anywhere -- N-alkyl
groups on amines/amides/ureas/hydrazines/guanidines, R-groups on
ammonium/phosphane/silane/borane/sulfonium, ester O-alkyl groups,
disulfide/diselenide R-groups, organomercury R-groups, and ordinary branch
substituents on a ordinary alcohol/acid/ketone main chain. It never once
consulted `chiral_tag` (the RDKit-derived CIP R/S code already computed
in molecule_analyzer.py and used correctly by the top-level acyclic chain
pipeline) -- so ANY molecule whose stereocenter happened to fall inside a
*substituent* rather than directly on the numbered principal chain lost
its (R)/(S) descriptor entirely and silently, e.g.:

    C[C@H](Cl)[N+](C)(C)C -> wrongly "(1-chloroethyl)trimethylazanium"
        (should be "[(S)-1-chloroethyl]trimethylazanium")
    OCCC(CCC)[C@H](Cl)C -> wrongly "3-(1-chloroethyl)hexan-1-ol"
        (should be "3-[(R)-1-chloroethyl]hexan-1-ol")

Fixed by computing the stereo descriptor list ONCE (reusing
`assign_stereochemistry`, the same function the main pipeline already
uses, called against the substituent's own `chain_path`) right after the
"simple linear alkyl substituent" chain_path is established inside
`_name_carbon_substituent`, and prepending it to every return path in
that branch (carboxyalkyl, alkynyl, alkenyl -- which already had its OWN
duplicate per-call computation, now deduplicated to reuse the single
shared value -- halo-, hydroxy-, amino-, sulfanyl-, oxo-substituted alkyl,
and the plain unsubstituted alkyl fallback).

All fixed names round-trip verified via OPSIN + RDKit canonical-SMILES
matching (including chirality tags, not just connectivity).

NOTE: two related, larger, NOT-yet-fixed gaps found during this
investigation and deliberately left for a future session (see project
memory): (1) `_name_branched_substituent` (the recursive branched-alkyl
namer, e.g. for isopropyl-shaped substituents) still has zero stereo
handling; (2) any BYPASS/standalone namer that builds its OWN principal
chain without going through `_name_carbon_substituent` at all -- e.g.
_name_thioamide, _name_substituted_amidine, and likely most of the
sulfonic/phosphonic/ester-family namers touched across Phases 902-909 --
still silently drops a stereocenter that lands directly on THEIR chain
(confirmed: `C[C@H](Cl)C(=S)N` -> wrongly "2-chloropropanethioamide",
no (R)/(S) at all). Both are structurally the same root cause (a
standalone chain-building code path that never calls
`assign_stereochemistry`) but are a much larger surface area than this
phase's fix (which only covers the shared `_name_carbon_substituent`
substituent-branch code path).
"""

from smiles2iupac import smiles_to_iupac


class TestSubstituentBranchStereoPreserved:
    def test_ammonium_r_group_stereocenter(self):
        assert (
            smiles_to_iupac("C[C@H](Cl)[N+](C)(C)C")
            == "[(S)-1-chloroethyl]trimethylazanium"
        )

    def test_phosphane_r_group_stereocenter(self):
        assert (
            smiles_to_iupac("C[C@H](Cl)P(C)C") == "[(S)-1-chloroethyl]dimethylphosphane"
        )

    def test_disulfide_r_group_stereocenter(self):
        assert (
            smiles_to_iupac("C[C@H](Cl)SSC") == "[(S)-1-chloroethyl] methyl disulfide"
        )

    def test_hydrazine_r_group_stereocenter(self):
        assert smiles_to_iupac("C[C@H](Cl)NN") == "(S)-1-chloroethylhydrazine"

    def test_organomercury_r_group_stereocenter(self):
        assert (
            smiles_to_iupac("C[C@H](Cl)[Hg]C") == "[(R)-1-chloroethyl]methylmercury"
        )

    def test_urea_n_substituent_stereocenter(self):
        assert smiles_to_iupac("C[C@H](Cl)NC(=O)N") == "N-[(S)-1-chloroethyl]urea"

    def test_nitrosamine_n_substituent_stereocenter(self):
        assert (
            smiles_to_iupac("C[C@H](Cl)N(C)N=O")
            == "N-[(S)-1-chloroethyl]-N-nitrosomethanamine"
        )

    def test_true_branch_off_main_alcohol_chain(self):
        # The chiral center is NOT on the numbered parent chain at all here --
        # it's two bonds into a genuine branch substituent.
        assert (
            smiles_to_iupac("OCCC(CCC)[C@H](Cl)C") == "3-[(R)-1-chloroethyl]hexan-1-ol"
        )

    def test_branch_stereocenter_on_acid_chain_context(self):
        assert smiles_to_iupac("C[C@H](Br)CC(=O)O") == "(3S)-3-bromobutanoic acid"

    def test_other_enantiomer(self):
        assert (
            smiles_to_iupac("CC[C@@H](F)[N+](C)(C)C")
            == "[(R)-1-fluoropropyl]trimethylazanium"
        )


class TestMainChainStereoRegression:
    def test_alanine_still_works(self):
        # Baseline: stereocenter directly ON the general acyclic pipeline's
        # own principal chain (not inside a helper-named substituent branch)
        # must be unaffected by this phase's change.
        assert smiles_to_iupac("C[C@H](N)C(=O)O") == "(2S)-2-aminopropanoic acid"

    def test_achiral_substituent_unaffected(self):
        assert smiles_to_iupac("C[C@H](Cl)C") == "2-chloropropane"
        assert smiles_to_iupac("CC(Cl)[N+](C)(C)C") == "(1-chloroethyl)trimethylazanium"
