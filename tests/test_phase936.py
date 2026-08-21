"""Phase 936: isotopic labeling on ring atoms (IUPAC 2013 P-82), extending
Phase 933's isotope descriptor mechanism from the acyclic path to rings.

Before this phase, isotope labels on ring atoms were silently dropped, same
shape as the pre-933 acyclic gap: `smiles_to_iupac("[2H]C1CCCCC1")` returned
plain "cyclohexane" with no indication of the deuterium label.

Scope: only isotope labels on the principal ring's atoms (reusing
`find_principal_ring`'s locant_map, mirroring `_name_acyclic`'s use of
`find_principal_chain`'s locant_map from Phase 933). Isotope labels on
substituent branches attached to a ring (e.g. an isotope-labeled ring that
itself appears as a substituent, as in an N-cyclohexyl group) remain out of
scope and are silently dropped, unchanged from before -- consistent with
Phase 933 leaving substituent-branch isotopes as a documented open gap.

Format and prefix ordering verified against OPSIN 2.9.0: same as the
acyclic case, the isotope descriptor sits after ordinary substituent
prefixes and before the ring parent name.

Non-obvious fix bundled into this phase: the pre-existing "single
substituent on an otherwise unmarked ring omits its '1-' locant" rule
(e.g. "methylcyclohexane", not "1-methylcyclohexane") had to be suppressed
whenever an isotope descriptor is also present. The ring's numbering is no
longer a free choice once the isotope descriptor pins a specific locant, so
omitting the substituent's locant produces a name OPSIN flags as
APPEARS_AMBIGUOUS ("methyl(4-2H)cyclohexane") even though it still
resolves correctly; "1-methyl(4-2H)cyclohexane" is unambiguous and is what
this phase now emits.
"""
from smiles2iupac import smiles_to_iupac


class TestPhase936RingIsotopeLabeling:
    def test_deuterated_cyclohexane(self):
        assert smiles_to_iupac("[2H]C1CCCCC1") == "(1-2H)cyclohexane"

    def test_deuterated_benzene(self):
        assert smiles_to_iupac("[2H]C1=CC=CC=C1") == "(1-2H)benzene"

    def test_carbon13_cyclohexane(self):
        assert smiles_to_iupac("C1[13CH2]CCCC1") == "(2-13C)cyclohexane"

    def test_deuterium_on_suffix_bearing_ring_atom(self):
        assert smiles_to_iupac("OC1([2H])CCCCC1") == "(1-2H)cyclohexanol"

    def test_deuterium_elsewhere_on_ring_with_suffix(self):
        assert smiles_to_iupac("OC1CCCC([2H])C1") == "(5-2H)cyclohexanol"

    def test_deuterated_cyclohexanone(self):
        assert smiles_to_iupac("O=C1CCC([2H])CC1") == "(4-2H)cyclohexanone"

    def test_deuterated_ring_diol(self):
        assert smiles_to_iupac("OC1([2H])CCC(O)CC1") == "(1-2H)cyclohexane-1,4-diol"

    def test_deuterated_ring_carboxamide(self):
        assert smiles_to_iupac("NC(=O)C1([2H])CCCCC1") == "(1-2H)cyclohexanecarboxamide"

    def test_substituent_locant_not_omitted_when_isotope_present(self):
        # Without the isotope label this would omit the "1-" for methyl
        # ("methylcyclohexane"); the isotope descriptor pins the ring's
        # numbering, so the substituent locant must stay explicit.
        assert smiles_to_iupac("CC1CCC([2H])CC1") == "1-methyl(4-2H)cyclohexane"

    def test_unlabeled_single_substituent_ring_still_omits_locant(self):
        # Regression guard: the pre-existing locant-omission rule still
        # applies when there is no isotope descriptor to pin the numbering.
        assert smiles_to_iupac("CC1CCCCC1") == "methylcyclohexane"

    def test_unlabeled_cyclohexane_unaffected(self):
        assert smiles_to_iupac("C1CCCCC1") == "cyclohexane"

    def test_ring_substituent_isotope_out_of_scope_at_time_of_this_phase(self):
        # At the time THIS phase (936) landed, an isotope-labeled ring
        # appearing as a substituent (not the principal ring itself) was
        # not yet covered -- the label was silently dropped. Phase 938
        # later added support for exactly this case (plain cycloalkyl
        # substituents); see test_phase938.py for the current behavior.
        # Kept here only as a historical marker of Phase 936's own scope.
        assert smiles_to_iupac("CC(=O)NC1([2H])CCCCC1") == "N-[(1-2H)cyclohexyl]acetamide"
