"""Phase 916: an ester, amide, or thioester group used as a SUBSTITUENT
(not the molecule's principal characteristic group) was silently
misidentified as "formyl" (a plain aldehyde, -CHO), dropping the
O-alkyl/N-substituent/S-alkyl side entirely.

    OC(=O)c1ccc(C(=O)OC)cc1 -> wrongly "4-formylbenzoic acid"
        (should be "4-methoxycarbonylbenzoic acid" -- monomethyl
        terephthalate, a real, common PET-production intermediate,
        confidently renamed to a structurally unrelated aldehyde)
    OC(=O)c1ccc(C(=O)N)cc1 -> wrongly "4-formylbenzoic acid"
        (should be "4-carbamoylbenzoic acid")
    OC(=O)c1ccc(C(=O)SC)cc1 -> wrongly "4-formylbenzoic acid"
        (should be "4-(methylsulfanylcarbonyl)benzoic acid")

Root cause: _name_carbon_substituent's acyl-as-substituent branch (used
whenever a root carbon is itself a carbonyl carbon being named as a
substituent, e.g. an acyl group on a ring) walked ONLY the carbon-only DFS
(_collect_substituent_carbons) to determine acyl chain length -- this DFS
cannot cross an O/N/S atom, so for -C(=O)-O-CH3 it only ever "saw" the
single carbonyl carbon itself (n_acyl == 1), triggering the
"n_acyl == 1 -> formyl" shortcut meant for a genuine aldehyde-derived
acyl group (-C(=O)-H). The O-methyl/N-H2/S-methyl side was never even
inspected.

Fixed by checking, BEFORE the benzoyl/formyl shortcuts, whether the
carbonyl carbon's other (non-excluded, non-double-bonded-O) neighbor is
an O/N/S heteroatom carrying its own carbon substituent (or, for N, is a
bare amide nitrogen) -- and building the correct substituent-prefix form
in each case: "(R-oxy)carbonyl" for esters (e.g. "methoxycarbonyl"),
"carbamoyl"/"N-alkylcarbamoyl"/"N,N-dialkylcarbamoyl" for amides, and
"(R-sulfanyl)carbonyl" for thioesters.

Found and fixed two follow-up bugs while verifying: (1) the carbamoyl
N-substituent scan initially forgot to exclude the acyl carbon itself
from "N's carbon neighbors" (root_idx is trivially adjacent to the amide
N), causing infinite-recursion-flavored garbage like
"N-benzoylcarbamoyl"; (2) "N-methylcarbamoyl"/"N,N-dimethylcarbamoyl"
need to always be parenthesized when used as a substituent, since the
embedded "N-" collides with a numeric ring locant immediately in front of
it (e.g. "4-N-methylcarbamoylbenzoic acid" fails to parse in OPSIN with
"Cannot find in scope fragment with atom with locant N4"), unlike the
digit-only ambiguity check (_needs_bis_tris) used elsewhere in the
codebase, which doesn't know about this specific "N-" collision.

All fixed names round-trip verified via OPSIN + RDKit canonical-SMILES
matching.
"""

from smiles2iupac import smiles_to_iupac


class TestAcylAsSubstituentHeteroatomAwareness:
    def test_ester_as_substituent_methyl(self):
        assert (
            smiles_to_iupac("OC(=O)c1ccc(C(=O)OC)cc1")
            == "4-methoxycarbonylbenzoic acid"
        )

    def test_ester_as_substituent_ethyl(self):
        assert (
            smiles_to_iupac("OC(=O)c1ccc(C(=O)OCC)cc1")
            == "4-ethoxycarbonylbenzoic acid"
        )

    def test_ester_as_substituent_other_direction(self):
        assert (
            smiles_to_iupac("COC(=O)c1ccc(C(=O)O)cc1")
            == "4-methoxycarbonylbenzoic acid"
        )

    def test_amide_as_substituent_plain(self):
        assert smiles_to_iupac("OC(=O)c1ccc(C(=O)N)cc1") == "4-carbamoylbenzoic acid"

    def test_amide_as_substituent_n_methyl(self):
        assert (
            smiles_to_iupac("OC(=O)c1ccc(C(=O)NC)cc1")
            == "4-(N-methylcarbamoyl)benzoic acid"
        )

    def test_amide_as_substituent_n_n_dimethyl(self):
        assert (
            smiles_to_iupac("OC(=O)c1ccc(C(=O)N(C)C)cc1")
            == "4-(N,N-dimethylcarbamoyl)benzoic acid"
        )

    def test_thioester_as_substituent(self):
        assert (
            smiles_to_iupac("OC(=O)c1ccc(C(=O)SC)cc1")
            == "4-(methylsulfanylcarbonyl)benzoic acid"
        )

    def test_ester_as_substituent_on_amide_parent(self):
        assert (
            smiles_to_iupac("NC(=O)c1ccc(C(=O)OC)cc1")
            == "4-methoxycarbonylbenzamide"
        )


class TestAcylAsSubstituentRegressions:
    def test_plain_ester_as_principal_group(self):
        assert smiles_to_iupac("c1ccc(cc1)C(=O)OC") == "methyl benzoate"

    def test_plain_amide_as_principal_group(self):
        assert smiles_to_iupac("c1ccc(cc1)C(=O)N") == "benzamide"

    def test_plain_ketone_as_substituent_unaffected(self):
        assert smiles_to_iupac("c1ccc(cc1)C(=O)C") == "1-phenylethanone"

    def test_plain_aldehyde_still_formyl(self):
        # A genuine aldehyde-as-substituent (root has ONLY H besides =O,
        # no O/N/S heteroatom neighbor) must still resolve to "formyl".
        assert smiles_to_iupac("OC(=O)c1ccc(C=O)cc1") == "4-formylbenzoic acid"

    def test_symmetric_diester_unaffected(self):
        # Both esters are principal groups here (dicarboxylate-style
        # naming), never reaching the acyl-as-substituent branch at all.
        assert (
            smiles_to_iupac("COC(=O)c1ccccc1C(=O)OCC")
            == "ethyl methyl benzene-1,2-dicarboxylate"
        )

    def test_aliphatic_ester_chain_unaffected(self):
        # The aliphatic (non-ring-substituent) case already worked before
        # this phase -- confirms the fix is scoped to the substituent
        # branch and doesn't regress the already-correct chain path.
        assert smiles_to_iupac("OC(=O)CCC(=O)OC") == "4-methoxy-4-oxobutanoic acid"
