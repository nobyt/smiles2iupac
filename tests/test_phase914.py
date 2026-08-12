"""Phase 914: ~25 standalone group_namers.py functions built their own
chain but never called assign_stereochemistry AT ALL (not even the
ene-gated version fixed in Phase 913) -- so any stereocenter on their own
chain was silently dropped unconditionally, e.g.:

    C[C@H](Cl)S(=O)(=O)N=[N+]=[N-] -> wrongly "1-chloroethanesulfonyl azide"
        (should be "(S)-1-chloroethanesulfonyl azide")
    C[C@H](Cl)C(=O)N=[N+]=[N-] -> wrongly "2-chloropropanoyl azide"
        (should be "(2S)-2-chloropropanoyl azide")

Found by re-running the Phase 913 cross-reference script (functions that
build a chain via _collect_acid_chain/_chain_through_pivot but never
reference assign_stereochemistry anywhere in their body) after Phase 913
closed the ene-gating bug -- the 25 functions it found were a genuinely
different, unaddressed set: _name_acyl_azide, _name_acyl_peroxide,
_name_acylhydrazone, _name_chalcogen_oxyacid (covers selenonic/seleninic/
selenenic/telluronic/tellurinic/tellurenic acid), _name_di_xisocyanate,
_name_diimine, _name_disulfonamide, _name_disulfonic_acid, _name_nitroso,
_name_o_substituted_oxime, _name_o_thioester, _name_peroxyester,
_name_s_dithioate_ester, _name_selenohydrazide,
_name_substituted_urea_if_match, _name_sulfenamide, _name_sulfenate_ester,
_name_sulfenyl_halide, _name_sulfinyl_chloride, _name_sulfinylhydrazide,
_name_sulfonate_anion, _name_sulfonimidamide, _name_sulfonohydrazide,
_name_sulfonyl_azide, _name_thiohydrazide.

Added a shared `_chain_stereo_prefix(graph, chain)` helper (wraps
assign_stereochemistry + PrincipalChain construction + the "(R,S,...)-"
join, returning "" when there's no stereo) to avoid repeating the same
6-line boilerplate in all 25 call sites, then spliced its result into
every return path of each function (mirroring the "compute once, use in
every branch" discipline Phase 913 established for the compounding bug it
found).

Also found and fixed a SEPARATE gap in the SAME investigation: the
acyl-as-substituent branch inside _name_carbon_substituent
(substituent.py) -- used whenever a whole acyl group like "2-chloro-
propanoyl" is itself a substituent, e.g. an O-acyl sulfenate ester or an
N-acyl imide -- has its OWN independent substituent-collection logic
(never routes through the plain-alkyl branch Phase 912 already fixed) and
also never checked chirality, e.g.
C[C@H](Cl)C(=O)OSC -> wrongly "2-chloropropanoyl methanesulfenate".

Also fixed a pre-existing, stereo-unrelated formatting bug found while
verifying _name_acylhydrazone: when the acid-chain side has its own
substituent (so its formatted name starts with a digit locant, e.g.
"2,2-dichloropropanohydrazide"), it was concatenated directly onto the
N'-alkylidene prefix with no separator, producing an ambiguous run-on
like "benzylidene2,2-dichloropropanohydrazide" instead of
"benzylidene-2,2-dichloropropanohydrazide".

All fixed names round-trip verified via OPSIN + RDKit canonical-SMILES
matching (except the O-substituted-oxime case, which is a pre-existing,
unrelated OPSIN grammar limitation already documented in Phase 903's
memory notes -- confirmed failing identically on the untouched/achiral
form too, not a regression from this phase).
"""

from smiles2iupac import smiles_to_iupac


class TestPreviouslyUnhandledChainStereo:
    def test_acyl_azide(self):
        assert smiles_to_iupac("C[C@H](Cl)C(=O)N=[N+]=[N-]") == "(2S)-2-chloropropanoyl azide"

    def test_acyl_peroxide(self):
        assert (
            smiles_to_iupac("C[C@H](Cl)C(=O)OOC(=O)C")
            == "(2S)-2-chloropropanoyl ethanoyl peroxide"
        )

    def test_chalcogen_oxyacid_seleninic(self):
        assert smiles_to_iupac("C[C@H](Cl)[Se](=O)O") == "(R)-1-chloroethaneseleninic acid"

    def test_nitroso(self):
        assert smiles_to_iupac("C[C@H](Cl)N=O") == "(S)-1-chloro-1-nitrosoethane"

    def test_peroxyester(self):
        assert (
            smiles_to_iupac("C[C@H](Cl)C(=O)OOC")
            == "methyl (2S)-2-chloropropaneperoxoate"
        )

    def test_s_dithioate_ester(self):
        assert (
            smiles_to_iupac("C[C@H](Cl)C(=S)SC")
            == "S-methyl (2S)-2-chloropropanedithioate"
        )

    def test_substituted_urea(self):
        assert smiles_to_iupac("C[C@H](Cl)NC(=O)N") == "N-[(S)-1-chloroethyl]urea"

    def test_sulfonyl_azide(self):
        assert (
            smiles_to_iupac("C[C@H](Cl)S(=O)(=O)N=[N+]=[N-]")
            == "(S)-1-chloroethanesulfonyl azide"
        )

    def test_sulfonohydrazide(self):
        assert (
            smiles_to_iupac("C[C@H](Cl)S(=O)(=O)NN") == "(S)-1-chloroethanesulfonohydrazide"
        )

    def test_sulfinylhydrazide(self):
        assert smiles_to_iupac("C[C@H](Cl)S(=O)NN") == "(S)-1-chloroethanesulfinylhydrazide"

    def test_sulfenamide(self):
        assert smiles_to_iupac("C[C@H](Cl)SN") == "(S)-1-chloroethanesulfenamide"

    def test_sulfonimidamide(self):
        assert (
            smiles_to_iupac("C[C@H](Cl)S(=O)(=N)N") == "(S)-1-chloroethanesulfonimidamide"
        )

    def test_thiohydrazide(self):
        assert smiles_to_iupac("C[C@H](Cl)C(=S)NN") == "(2S)-2-chloropropanethiohydrazide"

    def test_selenohydrazide(self):
        assert (
            smiles_to_iupac("C[C@H](Cl)C(=[Se])NN") == "(2S)-2-chloropropaneselenohydrazide"
        )


class TestAcylAsSubstituentStereo:
    def test_o_thioester_acyl_side(self):
        assert (
            smiles_to_iupac("C[C@H](Cl)C(=O)OSC")
            == "(2S)-2-chloropropanoyl methanesulfenate"
        )

    def test_n_acyl_imide_side(self):
        assert (
            smiles_to_iupac("CC(=O)N(C)C(=O)[C@H](Cl)C")
            == "N-[(2R)-2-chloropropanoyl]-N-methylacetamide"
        )


class TestAcylhydrazoneSeparatorFix:
    def test_substituted_acid_chain_needs_hyphen(self):
        # Pre-existing, stereo-independent bug: confirmed broken even with
        # NO stereocenter present (achiral CC(Cl)(Cl)... input).
        assert (
            smiles_to_iupac("CC(Cl)(Cl)C(=O)N/N=C/c1ccccc1")
            == "(E)-N'-benzylidene-2,2-dichloropropanohydrazide"
        )

    def test_plain_unsubstituted_regression(self):
        assert smiles_to_iupac("CC(=O)N/N=C/c1ccccc1") == "(E)-N'-benzylideneacetohydrazide"

    def test_stereocenter_plus_substituted_chain(self):
        assert (
            smiles_to_iupac("C[C@H](Cl)C(=O)N/N=C/c1ccccc1")
            == "(2S)-(E)-N'-benzylidene-2-chloropropanohydrazide"
        )


class TestAchiralRegressions:
    def test_plain_forms_unaffected(self):
        assert smiles_to_iupac("CS(=O)(=O)N=[N+]=[N-]") == "methanesulfonyl azide"
        assert smiles_to_iupac("CC(=O)OSC") == "acetyl methanesulfenate"
