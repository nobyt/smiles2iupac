"""Phase 938: isotopic labeling (IUPAC 2013 P-82) on substituent branches --
the last remaining item from docs/nomenclature_gap_audit.md's isotope
finding, after chains (933), rings (936), and the amide/amine/amidine/
imine special paths' own principal chains (937).

Before this phase, an isotope label on an atom that was itself part of a
*substituent* (as opposed to the principal chain/ring) was silently
dropped: `smiles_to_iupac("CC(=O)NCC[2H]")` returned plain
"N-ethylacetamide".

Scope: the two substituent-naming branches in `_name_carbon_substituent`
(`substituent.py`) that had an obvious, cheap locant map to reuse --
the simple linear-chain substituent path (methyl/ethyl/propyl.../
2-chloroethyl etc., root=1) and the plain (non-hetero) cycloalkyl
substituent path (cyclohexyl etc., root=1, ring walked out from root).
Isotope labels on more complex substituent shapes (branched substituents,
aryl substituents, heterocyclic substituents) remain out of scope and
continue to silently drop the label, unchanged from before.

Non-obvious formatting point, confirmed via OPSIN: the isotope descriptor
in a substituent name is itself parenthesized (`format_isotope_descriptor`
always returns a `(...)`-wrapped string), so when the substituent's own
name already needs enclosing marks (any name containing a digit needs
`_needs_bis_tris`-style wrapping at the call site), the *outer* wrapping
escalates to square brackets to avoid nested `(...)（...)` ambiguity --
e.g. "N-[(2-2H)ethyl]acetamide", not "N-((2-2H)ethyl)acetamide". This
escalation was already handled by pre-existing bracket-nesting logic
elsewhere in the codebase; this phase didn't need to add it.

A second coupled bug, same shape as Phase936's ring fix: the pre-existing
"single substituent on an otherwise-unmarked 2-carbon chain (ethane/
ethene) omits its own locant" rules in `assemble_name` (name_assembler.py)
didn't know about isotope descriptors either -- `"ClCC[2H]"` produced
`"chloro(2-2H)ethane"`, which OPSIN parses but flags APPEARS_AMBIGUOUS.
Fixed by adding `and not isotope_descriptor` to those rules' guards,
mirroring Phase936's ring-side fix exactly.
"""
from smiles2iupac import smiles_to_iupac


class TestPhase938IsotopeLinearSubstituent:
    def test_deuterated_n_ethyl_substituent(self):
        assert smiles_to_iupac("CC(=O)NCC[2H]") == "N-[(2-2H)ethyl]acetamide"

    def test_deuterated_n_methyl_substituent(self):
        assert smiles_to_iupac("CC(=O)NC[2H]") == "N-[(1-2H)methyl]acetamide"

    def test_deuterated_substituent_combined_with_halogen(self):
        assert smiles_to_iupac("CC(=O)NCC(Cl)[2H]") == "N-[2-chloro(2-2H)ethyl]acetamide"

    def test_deuterated_substituent_on_amine(self):
        assert smiles_to_iupac("CNC[2H]") == "N-[(1-2H)methyl]methanamine"

    def test_deuterated_substituent_among_multiple_n_subs(self):
        assert smiles_to_iupac("CC(=O)N(C)CC[2H]") == "N-[(2-2H)ethyl]-N-methylacetamide"

    def test_unlabeled_substituent_unaffected(self):
        assert smiles_to_iupac("CC(=O)NCC") == "N-ethylacetamide"


class TestPhase938IsotopeCycloalkylSubstituent:
    def test_deuterated_cyclohexyl_at_attachment_point(self):
        assert smiles_to_iupac("CC(=O)NC1([2H])CCCCC1") == "N-[(1-2H)cyclohexyl]acetamide"

    def test_deuterated_cyclohexyl_away_from_attachment(self):
        assert smiles_to_iupac("CC(=O)NC1CCC([2H])CC1") == "N-[(4-2H)cyclohexyl]acetamide"

    def test_unlabeled_cyclohexyl_substituent_unaffected(self):
        assert smiles_to_iupac("CC(=O)NC1CCCCC1") == "N-cyclohexylacetamide"


class TestPhase938IsotopeChainLocantOmissionRegression:
    def test_labeled_chloroethane_keeps_explicit_locant(self):
        # Without the isotope, this omits the "1-" ("chloroethane"); the
        # isotope descriptor pins the chain's numbering, so the
        # substituent locant must stay explicit to avoid ambiguity.
        assert smiles_to_iupac("ClCC[2H]") == "1-chloro(2-2H)ethane"

    def test_labeled_chloroethene_keeps_explicit_locant(self):
        assert smiles_to_iupac("[2H]C=CCl") == "1-chloro(2-2H)ethene"

    def test_unlabeled_chloroethane_still_omits_locant(self):
        assert smiles_to_iupac("ClCC") == "chloroethane"

    def test_unlabeled_chloroethene_still_omits_locant(self):
        assert smiles_to_iupac("C=CCl") == "chloroethene"
