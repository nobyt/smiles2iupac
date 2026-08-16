"""Phase 926: a crash (missing FUNCTIONAL_GROUPS entry for a merged
group type) and the amide-family twin of Phase 925's amine-merge
substituent-dropping bug, both found via the continuing sweep of
_multi_map-adjacent group types.

1. `_multi_map` has an entry `("amidine", 2): "diamidine"`, but
   `constants.py`'s `FUNCTIONAL_GROUPS` (which `SUFFIX_MAP` is derived
   from) never had a matching "diamidine" entry -- any molecule with 2
   amidine groups crashed with `KeyError: 'diamidine'` in
   `name_assembler.py`'s `assemble_name`. Fixed by adding a "diamidine"
   `FunctionalGroupSpec` (suffix "diimidamide", doubling amidine's own
   irregular "imidamide" suffix, following the same C1-and-Cn-anchored
   pattern already used for dinitrile/dioic_acid/dial) plus the matching
   `_build_name_body` branch in name_assembler.py.

       CC(C(=N)N)(C(=N)N)CC -> crashed
           -> now "2-ethyl-2-methylpropanediimidamide"

2. Same root mechanism as Phase925's amine-merge bug, but amide's
   `atom_indices` shape ([carbonyl_C, O, N], ALWAYS 3 elements
   regardless of N-substitution -- unlike amine, whose shape varies)
   meant the Phase925-style "check atom_indices length" guard doesn't
   apply here. N-substituted (secondary/tertiary) amides were merging
   into "diamide" purely by count, and since the merged group's N atoms
   land in the collect_substituents exclusion set, each N's OWN N-alkyl
   substituents became unreachable and silently vanished:

       CNC(=O)C(C(=O)NC)CC -> wrongly involved dropped N-methyls
           (both instances lost their N-methyl group entirely)

   Fixed with a new `aggregate_groups()` guard (functional_group.py):
   checks whether each amide instance's N has any carbon neighbor OTHER
   than its own carbonyl carbon (i.e. is N-substituted); if any instance
   is, skip the diamide merge entirely, routing through the pre-existing
   `_name_secondary_tertiary_amide` pathway instead.

   This surfaced a THIRD, previously-unreachable bug: when the two
   N-substituted amides are chain-connected (not just branch-connected),
   `_name_secondary_tertiary_amide`'s acid-side chain collector
   (`_collect_acid_chain`) walked straight through the SECOND amide's
   own carbonyl carbon as if it were an ordinary chain carbon, silently
   absorbing the entire second amide into the "parent" acid chain and
   producing a structurally wrong name (verified via RDKit --
   "N-methyl-3-oxo-3-propanamidopropanamide" for
   `CNC(=O)CC(=O)NC" parsed to a different molecule than the input).
   Fixed by excluding any OTHER amide's carbonyl carbon (detected via
   the existing `_is_amide` check) from the acid-chain walk, so the
   chain correctly stops there and the second amide is picked up as an
   ordinary "(N-methylcarbamoyl)"-style substituent by the existing,
   already-correct carbamoyl-substituent naming machinery instead.

All fixed names round-trip verified via OPSIN + RDKit canonical-SMILES
matching. Full suite (13511 tests) run clean, 0 failures.

NOT fixed this phase (found while probing around this bug, documented
in project memory for a dedicated future session): `_name_carbon_
substituent`'s "simple linear alkyl substituent" branch checks for
halogen / hydroxy / amino / sulfanyl / oxo substituents on a chain
SEQUENTIALLY with early returns, so when 2+ of these types co-occur on
the SAME substituent-of-a-substituent chain, only the first-checked
type survives -- e.g. `CC(=O)NCC(=O)N` (an N-substituent chain with a
terminal amide, i.e. both oxo AND amino on the same carbon) produces
"N-(2-aminoethyl)acetamide", silently dropping the oxo and changing the
functional group entirely (amide -> amine). This is a materially larger
fix (combining 5 independent early-return checks into one accumulating
pass) than anything else in this phase and was deliberately not
attempted here.
"""

from smiles2iupac import smiles_to_iupac


class TestDiamidineCrash:
    def test_branching_diamidine_no_longer_crashes(self):
        assert (
            smiles_to_iupac("CC(C(=N)N)(C(=N)N)CC")
            == "2-ethyl-2-methylpropanediimidamide"
        )

    def test_unbranched_diamidine(self):
        assert smiles_to_iupac("NC(=N)CCC(=N)N") == "butanediimidamide"


class TestSecondaryTertiaryAmideNoLongerMerges:
    def test_branch_connected_secondary_diamide(self):
        assert (
            smiles_to_iupac("CNC(=O)C(C(=O)NC)CC")
            == "2-(N-methylcarbamoyl)-N-methylbutanamide"
        )

    def test_chain_connected_symmetric_secondary_diamide(self):
        assert (
            smiles_to_iupac("CNC(=O)CC(=O)NC")
            == "2-(N-methylcarbamoyl)-N-methylacetamide"
        )

    def test_chain_connected_secondary_diamide_longer(self):
        assert (
            smiles_to_iupac("CNC(=O)CCC(=O)NC")
            == "3-(N-methylcarbamoyl)-N-methylpropanamide"
        )


class TestAmideMergeRegressions:
    def test_n_acylamide_shared_n_regression(self):
        assert smiles_to_iupac("CC(=O)NC(=O)C") == "N-acetylacetamide"

    def test_plain_primary_diamide_still_merges(self):
        assert smiles_to_iupac("NC(=O)CC(=O)N") == "propanediamide"

    def test_simple_secondary_amide(self):
        assert smiles_to_iupac("CC(=O)NC") == "N-methylacetamide"

    def test_simple_tertiary_amide(self):
        assert smiles_to_iupac("CC(=O)N(C)C") == "N,N-dimethylacetamide"

    def test_anilide_regression(self):
        assert smiles_to_iupac("CC(=O)Nc1ccccc1") == "N-phenylacetamide"
