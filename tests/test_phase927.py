"""Phase 927: fixed the high-value bug flagged (but deliberately
deferred, given the scope/risk) at the end of Phase 926 -- a
substituent-of-a-substituent chain carrying 2+ of
{halogen, hydroxy, amino, sulfanyl, oxo} silently lost all but one.

`substituent.py`'s `_name_carbon_substituent` "simple linear alkyl
substituent" branch (the one that builds names like "2-hydroxyethyl",
"aminomethyl", "trifluoromethyl" for a chain that is itself a
substituent, e.g. an amide's N-alkyl group) checked for these 5 terminal
substituent types as 5 INDEPENDENT blocks, each of which returned
IMMEDIATELY upon finding ANY match of its own type -- checked in the
fixed priority order halogen > hydroxy > amino > sulfanyl > oxo. When 2+
of these types co-occurred anywhere on the same chain, only the
FIRST-checked type survived; everything else silently vanished from the
name:

    CC(=O)NCC(=O)N -> wrongly "N-(2-aminoethyl)acetamide"
        (the chain terminal has BOTH oxo AND amino -- i.e. it's a
        carboxamide terminus, -CH2-C(=O)NH2 -- but "amino" was checked
        before "oxo" ever got a chance, so the oxo silently vanished and
        the substituent read as a plain amine instead of an amide)

    CC(=O)NCC(O)Cl -> wrongly "N-(2-chloroethyl)acetamide"
        (halogen checked before hydroxy; hydroxy silently dropped)

This is structurally different from every other bug fixed in the
921-926 arc: those were all narrow, `_multi_map`/merge-related issues
confined to specific functional-group families. This one is not
confined to any group family at all -- it can affect ANY molecule where
a substituent-of-a-substituent chain happens to carry 2+ of these 5
types anywhere along its length, an essentially unbounded combinatorial
surface (found via probing the exact `_name_secondary_tertiary_amide`
code path unblocked by Phase926's amide-merge-prevention fix, but the
underlying bug long predates Phase926 and is unrelated to amides
specifically).

Fixed by combining the 5 independent early-return checks into ONE
accumulating pass: for every chain position, collect ALL matching
substituent types (not just the first found) into one list, THEN build
the combined prefix string ONCE at the end -- grouping and sorting
exactly as each individual block already did, just applied across all 5
types together instead of one type in isolation. Every individual
per-type detection condition (hydroxy requires a free O-H with no other
carbon i.e. not an ether; sulfanyl the same for S-H; amino requires an
unsubstituted, non-ring, non-nitrile N; oxo requires a C=O double bond;
halogen is a simple symbol match) was preserved EXACTLY as before --
only the control flow (early-return-per-type -> accumulate-then-format-
once) changed.

All fixed names round-trip verified via OPSIN + RDKit canonical-SMILES
matching. Given this touches an extremely widely-shared, core
substituent-naming function used by essentially every branch/N-
substituent/O-substituent chain across the entire codebase, the full
pre-existing suite (13521 tests) was run in ISOLATION first, before this
test file was written, specifically to separate "did this change break
anything already-established" from "does the new test file also pass" --
passed clean, 0 failures, before this file was added.

NOT fixed this phase, noted for a smaller future follow-up: the ene/yne
(double/triple bond) branch of the same function, which sits ABOVE this
phase's combined block and ALSO returns early, has the identical bug
shape for a chain that has BOTH a multiple bond AND a heteroatom
substituent (e.g. `CC(=O)NCC=CCO` -> wrongly "N-(but-2-en-1-yl)-
acetamide", dropping the terminal hydroxyl). Smaller in scope than the
bug fixed here (one remaining early-return branch rather than five), but
deliberately not folded into this phase to avoid extending an
already-large, already-verified change mid-flight.
"""

from smiles2iupac import smiles_to_iupac


class TestCoOccurringSubstituentTypes:
    def test_oxo_and_amino_terminal_amide_shape(self):
        assert (
            smiles_to_iupac("CC(=O)NCC(=O)N")
            == "N-(2-amino-2-oxoethyl)acetamide"
        )

    def test_halogen_and_hydroxy(self):
        assert (
            smiles_to_iupac("CC(=O)NCC(O)Cl")
            == "N-(2-chloro-2-hydroxyethyl)acetamide"
        )

    def test_amino_and_hydroxy(self):
        assert (
            smiles_to_iupac("CC(=O)NCC(N)O")
            == "N-(2-amino-2-hydroxyethyl)acetamide"
        )

    def test_sulfanyl_and_hydroxy(self):
        assert (
            smiles_to_iupac("CC(=O)NCC(S)O")
            == "N-(2-hydroxy-2-sulfanylethyl)acetamide"
        )

    def test_halogen_and_hydroxy_different_positions(self):
        assert (
            smiles_to_iupac("CC(=O)NCC(Cl)CO")
            == "N-(2-chloro-3-hydroxypropyl)acetamide"
        )


class TestSingleTypeRegressions:
    def test_halogen_only(self):
        assert smiles_to_iupac("CC(=O)NCCCCl") == "N-(3-chloropropyl)acetamide"

    def test_multiple_identical_halogens(self):
        assert (
            smiles_to_iupac("CC(=O)NC(Cl)CCl") == "N-(1,2-dichloroethyl)acetamide"
        )

    def test_trifluoromethyl(self):
        assert (
            smiles_to_iupac("CC(=O)NC(F)(F)F") == "N-(trifluoromethyl)acetamide"
        )

    def test_plain_secondary_amide(self):
        assert smiles_to_iupac("CC(=O)NC") == "N-methylacetamide"


class TestBranchingPolyolRegressionsUnaffectedByThisChange:
    # These went through the SAME code path (hydroxy-detection) fixed in
    # Phase921/922 -- must remain unaffected by combining the checks.
    def test_pentaerythritol(self):
        assert (
            smiles_to_iupac("OCC(CO)(CO)CO")
            == "2,2-bis(hydroxymethyl)propane-1,3-diol"
        )

    def test_trimethylolpropane(self):
        assert (
            smiles_to_iupac("CCC(CO)(CO)CO")
            == "2-ethyl-2-(hydroxymethyl)propane-1,3-diol"
        )

    def test_gem_triamine(self):
        assert smiles_to_iupac("CCC(N)(N)N") == "propane-1,1,1-triamine"
