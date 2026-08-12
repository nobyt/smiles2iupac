"""Phase 911: systematic sweep for substituent-prefix mis-alphabetization.

Triggered by noticing _name_by_c_substituents (ammonium/phosphane/silane/
phosphinic-acid-family shared namer) sorted pre-formatted substituent name
strings like "2,2,2-trifluoroethyl" with a plain sorted(), which alphabetizes
by the leading digit ('2' < 'm') instead of the substituent word ('methyl' <
'trifluoroethyl') -- e.g. C[N+](C)(C)CC(F)(F)F -> wrongly
"2,2,2-trifluoroethyltrimethylazanium" instead of
"trimethyl(2,2,2-trifluoroethyl)azanium" (also missing the disambiguating
parens a locant-bearing substituent needs).

Grepping every `sorted(` call in group_namers.py touching substituent-name
strings (as opposed to integer locants) found the identical bug shape
repeated across ~50 call sites spanning many unrelated functional-group
families: ammonium/phosphane/silane/borane (_name_by_c_substituents),
phosphinate/phosphonothioate/phosphonate/phosphate/phosphite esters,
the phosphoramidic/phosphonimidic acid family, organomercury,
disulfide/diselenide/selenoxide/sulfoxide (_dual_c_group_prefix's four
copies), disilane/disiloxane/disilazane, sulfamide, hydrazine compounds,
carbodiimide, guanidine/amidine, secondary/tertiary amines and amides,
nitrosamine, thioamide, cyanamide, sulfate/sulfite esters, and acyl
peroxide. Fixed via one shared `_substituent_alpha_key` helper (strips a
leading "(" and any "[\\d,]+-" locant prefix before comparing -- but does
NOT strip "di/tri/tetra" the way name_assembler._build_prefix's own
alpha_key does, since that also incorrectly strips "tri" from a composite
name like "trifluoroethyl" that only coincidentally starts with a
multiplier syllable).

Also found and fixed two deeper bugs surfaced while regression-testing:
- _name_secondary_tertiary_amine's own chain substituents were joined with
  a naive one-per-atom loop instead of _build_prefix-style grouping, so a
  CF3 group exploded into "2-fluoro-2-fluoro-2-fluoro" instead of merging
  to "2,2,2-trifluoro" (CNCC(F)(F)F -> wrongly
  "2-fluoro-2-fluoro-2-fluoro-N-methylethanamine").
- The amine/amide functions' own local _alpha_key/_alpha_key_amide closures
  (which merge chain_sub_parts + n_prefix_parts into one alphabetized list)
  ran the digit-locant strip BEFORE stripping the "N-"/"N,N-" prefix, so it
  never matched (the digit lives after "N-", not at the very start of the
  part string) -- fixed by reordering strips. The di/tri-multiplier strip
  itself needed a MORE careful split than the general
  _substituent_alpha_key helper: chain_sub_parts entries (no "N-" prefix)
  always self-generate their own multiplying prefix in THIS function's own
  grouping code (e.g. "2,2,2-trifluoro" = 3x plain "fluoro" merged), so
  stripping "tri" there is always correct and alphabetizes as "fluoro"
  (confirmed by an existing Phase905 test: CNC(=O)C(F)(F)F ==
  "2,2,2-trifluoro-N-methylacetamide", fluoro < methyl) -- but a bare
  "N-(2,2,2-trifluoroethyl)" IS a single composite substituent whose name
  merely starts with "tri", and must alphabetize as "trifluoroethyl"
  (methyl < trifluoroethyl). Guarded so only genuinely-multiplied N-parts
  ("N,N-" with a comma) get the strip, chain parts always do.
- _name_nitrosamine's parent-chain-length heuristic (_chain_len) also
  failed to recognize a substituted chain name as being that chain's
  length at all (returned 0), only fixed for the narrow case where a
  digit-locant is followed directly by the stem with no other prefix word.

All new/changed names round-trip verified via OPSIN + RDKit canonical-
SMILES matching.
"""

from smiles2iupac import smiles_to_iupac


class TestAmmoniumPhosphanesSilanesBoranes:
    def test_ammonium_mixed_alkyl(self):
        assert smiles_to_iupac("C[N+](C)(C)CC(F)(F)F") == "trimethyl(2,2,2-trifluoroethyl)azanium"

    def test_ammonium_bis_identical_complex(self):
        assert (
            smiles_to_iupac("C[N+](C)(CC(F)(F)F)CC(F)(F)F")
            == "dimethylbis(2,2,2-trifluoroethyl)azanium"
        )

    def test_ammonium_tris(self):
        assert (
            smiles_to_iupac("C[N+](CC(F)(F)F)(CC(F)(F)F)CC(F)(F)F")
            == "methyltris(2,2,2-trifluoroethyl)azanium"
        )

    def test_ammonium_plain_regression(self):
        assert smiles_to_iupac("C[N+](C)(C)C") == "tetramethylazanium"

    def test_phosphane_mixed(self):
        assert smiles_to_iupac("CP(C)CC(F)(F)F") == "dimethyl(2,2,2-trifluoroethyl)phosphane"

    def test_phosphane_branched_no_locant(self):
        assert smiles_to_iupac("CC(C)P(C)C") == "dimethyl(propan-2-yl)phosphane"

    def test_silane_mixed(self):
        assert smiles_to_iupac("C[Si](C)(C)CC(F)(F)F") == "trimethyl(2,2,2-trifluoroethyl)silane"

    def test_borane_halo_alkyl_mixed(self):
        assert smiles_to_iupac("CB(Cl)CC(F)(F)F") == "chloro(methyl)(2,2,2-trifluoroethyl)borane"

    def test_borane_alkyl_only_mixed(self):
        assert smiles_to_iupac("CB(C)CC(F)(F)F") == "dimethyl(2,2,2-trifluoroethyl)borane"


class TestPhosphorusEsterFamily:
    def test_phosphinate_ester_mixed(self):
        assert (
            smiles_to_iupac("COP(=O)(C)CC(F)(F)F")
            == "methyl methyl(2,2,2-trifluoroethyl)phosphinate"
        )

    def test_phosphine_imine_mixed(self):
        assert (
            smiles_to_iupac("CN=P(C)CC(F)(F)F")
            == "N-methylmethyl(2,2,2-trifluoroethyl)phosphanimine"
        )

    def test_phosphate_triester_mixed(self):
        assert (
            smiles_to_iupac("COP(=O)(OC)OCC(F)(F)F")
            == "dimethyl 2,2,2-trifluoroethyl phosphate"
        )

    def test_phosphorodiamidic_mixed(self):
        assert (
            smiles_to_iupac("CNP(=O)(NCC(F)(F)F)O")
            == "N-methyl-N'-(2,2,2-trifluoroethyl)phosphorodiamidic acid"
        )

    def test_phosphonimidic_mixed(self):
        assert (
            smiles_to_iupac("CP(=NCC(F)(F)F)(O)O")
            == "P-methyl-N-(2,2,2-trifluoroethyl)phosphonimidic acid"
        )

    def test_phosphonimidic_plain_regression(self):
        assert smiles_to_iupac("CP(=N)(O)O") == "P-methylphosphonimidic acid"


class TestOrganomercury:
    def test_mixed_no_halogen(self):
        assert smiles_to_iupac("C[Hg]CC(F)(F)F") == "methyl(2,2,2-trifluoroethyl)mercury"


class TestChalcogenPairFamily:
    def test_disulfide_mixed(self):
        assert smiles_to_iupac("CSSCC(F)(F)F") == "methyl (2,2,2-trifluoroethyl) disulfide"

    def test_diselenide_mixed(self):
        assert smiles_to_iupac("C[Se][Se]CC(F)(F)F") == "methyl (2,2,2-trifluoroethyl) diselenide"

    def test_selenoxide_mixed(self):
        assert smiles_to_iupac("C[Se](=O)CC(F)(F)F") == "methyl (2,2,2-trifluoroethyl) selenoxide"

    def test_sulfoxide_mixed(self):
        assert smiles_to_iupac("C[S](=O)CC(F)(F)F") == "methyl (2,2,2-trifluoroethyl) sulfoxide"

    def test_disulfide_plain_regression(self):
        assert smiles_to_iupac("CSSC") == "dimethyl disulfide"


class TestAmineChainSubstituentGrouping:
    def test_secondary_amine_cf3_groups_not_exploded(self):
        assert smiles_to_iupac("CNCC(F)(F)F") == "2,2,2-trifluoro-N-methylethanamine"

    def test_tertiary_amine_nn_cf3_groups_not_exploded(self):
        assert smiles_to_iupac("CN(C)CC(F)(F)F") == "2,2,2-trifluoro-N,N-dimethylethanamine"

    def test_amine_three_way_mixed_alpha_order(self):
        assert (
            smiles_to_iupac("CN(CC)CC(F)(F)F")
            == "N-methyl-N-(2,2,2-trifluoroethyl)ethanamine"
        )

    def test_amine_plain_regression(self):
        assert smiles_to_iupac("CNC") == "N-methylmethanamine"
        assert smiles_to_iupac("CN(C)C") == "N,N-dimethylmethanamine"


class TestAmideAlphaOrder:
    def test_tertiary_amide_mixed(self):
        assert (
            smiles_to_iupac("CC(=O)N(C)CC(F)(F)F")
            == "N-methyl-N-(2,2,2-trifluoroethyl)acetamide"
        )

    def test_amide_plain_regression(self):
        assert smiles_to_iupac("CC(=O)NC") == "N-methylacetamide"
        assert smiles_to_iupac("CC(=O)N(C)C") == "N,N-dimethylacetamide"

    def test_thioamide_mixed(self):
        assert smiles_to_iupac("CC(=S)NCC(F)(F)F") == "N-2,2,2-trifluoroethylethanethioamide"


class TestSulfateSulfiteEsters:
    def test_sulfate_mixed(self):
        assert (
            smiles_to_iupac("COS(=O)(=O)OCC(F)(F)F")
            == "methyl 2,2,2-trifluoroethyl sulfate"
        )

    def test_sulfite_mixed(self):
        assert smiles_to_iupac("COS(=O)OCC(F)(F)F") == "methyl 2,2,2-trifluoroethyl sulfite"

    def test_sulfate_plain_regression(self):
        assert smiles_to_iupac("COS(=O)(=O)OC") == "dimethyl sulfate"


class TestAcylPeroxideAlphaOrder:
    def test_mixed_acyl_peroxide(self):
        assert (
            smiles_to_iupac("CC(=O)OOC(=O)CC(F)(F)F")
            == "ethanoyl 3,3,3-trifluoropropanoyl peroxide"
        )


class TestNitrosamine:
    def test_mixed_alpha_order_and_parent_selection(self):
        # Also exercises the _chain_len fix: without it "methyl" (chain_len=1)
        # incorrectly beat "2,2,2-trifluoroethyl" (chain_len bugged to 0) for
        # parent-chain selection.
        assert (
            smiles_to_iupac("CN(CC(F)(F)F)N=O")
            == "N-nitroso-N-2,2,2-trifluoroethylmethanamine"
        )

    def test_plain_regression(self):
        assert smiles_to_iupac("CN(N=O)C") == "N-methyl-N-nitrosomethanamine"
