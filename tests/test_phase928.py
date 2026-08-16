"""Phase 928: fixed the small residual finding noted at the end of
Phase 927 -- the ene/yne (multiple-bond) branch of
`_name_carbon_substituent`'s linear-chain naming logic had the exact
same "early return drops co-occurring substituents" bug shape that
Phase927 fixed for the halogen/hydroxy/amino/sulfanyl/oxo checks, just
in a smaller, separate branch that sits ABOVE (is checked before) them.

    CC(=O)NCC=CCO -> wrongly "N-(but-2-en-1-yl)acetamide"
        (the chain has a C=C double bond AND a terminal -CH2OH; the
        ene branch detected the double bond and returned immediately,
        never reaching Phase927's combined heteroatom-substituent scan
        at all -- the hydroxyl silently vanished)

Fixed by no longer returning immediately upon finding a multiple bond:
instead, the "yl"-suffix word (plain "propyl"-style "yl", or the
ene/yne-specific "but-2-en-1-yl"/"but-2-yn-1-yl" style) is computed
into a local variable, and the function falls through into Phase927's
already-fixed heteroatom-substituent accumulation, which now uses that
computed suffix (instead of a hardcoded plain "{prefix}yl") for both of
its return points. Every existing ene/yne locant-computation rule
(unlocanted "ethenyl"/"ethynyl" for a 2-carbon chain, explicit
"prop-2-en-1-yl"-style locants otherwise) was preserved exactly --
only the control flow changed from early-return to fall-through.

All fixed names round-trip verified via OPSIN + RDKit canonical-SMILES
matching. Full suite (13533 tests) run clean before this test file was
added (same "isolate the source change's regression signal" discipline
used for Phase927, given this touches the same core, widely-shared
function), then again with this file included.
"""

from smiles2iupac import smiles_to_iupac


class TestEneYneWithCoOccurringSubstituents:
    def test_ene_with_hydroxy(self):
        assert (
            smiles_to_iupac("CC(=O)NCC=CCO")
            == "N-(4-hydroxybut-2-en-1-yl)acetamide"
        )

    def test_ene_with_halogen(self):
        assert (
            smiles_to_iupac("CC(=O)NCC=CCCl")
            == "N-(4-chlorobut-2-en-1-yl)acetamide"
        )

    def test_ene_with_oxo(self):
        assert (
            smiles_to_iupac("CC(=O)NCC=CC(=O)C")
            == "N-(4-oxopent-2-en-1-yl)acetamide"
        )

    def test_yne_with_hydroxy(self):
        assert (
            smiles_to_iupac("CC(=O)NCC#CCO")
            == "N-(4-hydroxybut-2-yn-1-yl)acetamide"
        )


class TestEneYneRegressions:
    def test_plain_ene_no_heteroatom(self):
        assert smiles_to_iupac("CC(=O)NCC=CC") == "N-(but-2-en-1-yl)acetamide"

    def test_short_ene_unlocanted(self):
        assert smiles_to_iupac("CC(=O)NCC=C") == "N-(prop-2-en-1-yl)acetamide"

    def test_ethenyl(self):
        assert smiles_to_iupac("CC(=O)NC=C") == "N-ethenylacetamide"

    def test_ethynyl(self):
        assert smiles_to_iupac("CC(=O)NC#C") == "N-ethynylacetamide"

    def test_main_chain_alkene_unaffected(self):
        assert smiles_to_iupac("C=CC(Cl)C") == "3-chlorobut-1-ene"
