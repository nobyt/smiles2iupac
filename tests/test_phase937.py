"""Phase 937: isotopic labeling (IUPAC 2013 P-82) on the 4 dedicated
early-return naming paths in `_name_acyclic` (secondary/tertiary amide,
secondary/tertiary amine, substituted amidine, N-substituted imine).

Before this phase, these 4 paths -- listed as "still open" in
docs/nomenclature_gap_audit.md since Phase 933 -- silently dropped any
isotope label, same shape as the pre-933/936 chain/ring gaps:
`smiles_to_iupac("CC(=O)NC[2H]")` returned plain "N-methylacetamide".

Scope: only isotope labels on each path's own principal-chain atoms
(the acid-side chain for amide/amidine, the parent alkyl chain for amine,
`find_principal_chain`'s chain for imine -- which already routes through
`assemble_name`, so it only needed the new isotope_descriptor param wired
through). Isotope labels on N-substituents (the "other" branch, e.g. the
N-methyl in an N-methylamide) were out of scope for THIS phase, consistent
with Phase933/936's substituent-branch boundary at the time -- Phase 938
(a later, separate phase) added substituent-branch support for the simple
linear-chain/cycloalkyl cases, so the N-substituent tests below now assert
the Phase938-fixed behavior rather than the drop; see test_phase938.py.

All fixed names verified via OPSIN 2.9.0 CLI round-trip + RDKit InChI
structural-equivalence comparison against the original SMILES before
being added here.
"""
from smiles2iupac import smiles_to_iupac


class TestPhase937IsotopeAmide:
    def test_deuterated_acid_chain(self):
        assert smiles_to_iupac("CC([2H])C(=O)NC") == "N-methyl(2-2H)propanamide"

    def test_deuterated_formamide_carbonyl_carbon(self):
        assert smiles_to_iupac("[2H]C(=O)NC") == "N-methyl(1-2H)formamide"

    def test_unlabeled_amide_unaffected(self):
        assert smiles_to_iupac("CC(=O)NC") == "N-methylacetamide"

    def test_n_substituent_isotope_out_of_scope_at_time_of_this_phase(self):
        # At the time THIS phase (937) landed, N-substituent isotopes were
        # not yet covered. Phase 938 later added support for the simple
        # linear-chain substituent case; see test_phase938.py.
        assert smiles_to_iupac("CC(=O)NC[2H]") == "N-[(1-2H)methyl]acetamide"


class TestPhase937IsotopeAmine:
    def test_deuterated_parent_chain(self):
        assert smiles_to_iupac("[2H]CNC") == "N-methyl(1-2H)methanamine"

    def test_deuterated_parent_chain_ethyl(self):
        assert smiles_to_iupac("CC([2H])NC") == "N-methyl(1-2H)ethanamine"

    def test_unlabeled_amine_unaffected(self):
        assert smiles_to_iupac("CNC") == "N-methylmethanamine"


class TestPhase937IsotopeAmidine:
    def test_deuterated_acid_chain(self):
        assert smiles_to_iupac("[2H]CC(=N)N") == "(2-2H)ethanimidamide"

    def test_deuterated_acid_chain_branched(self):
        assert smiles_to_iupac("CC([2H])C(=N)N") == "(2-2H)propanimidamide"

    def test_unlabeled_amidine_unaffected(self):
        assert smiles_to_iupac("CC(=N)N") == "ethanimidamide"

    def test_n_substituent_isotope_out_of_scope_at_time_of_this_phase(self):
        # At the time THIS phase (937) landed, N-substituent isotopes were
        # not yet covered. Phase 938 later added support for the simple
        # linear-chain substituent case; see test_phase938.py.
        assert smiles_to_iupac("CC(=N)NC[2H]") == "N-[(1-2H)methyl]ethanimidamide"


class TestPhase937IsotopeImine:
    def test_deuterated_principal_chain(self):
        assert smiles_to_iupac("[2H]C(=NC)C") == "N-methyl(1-2H)ethanimine"

    def test_deuterated_principal_chain_longer(self):
        assert smiles_to_iupac("CC([2H])C(=NC)C") == "N-methyl(3-2H)butan-2-imine"

    def test_unlabeled_imine_unaffected(self):
        assert smiles_to_iupac("C(=NC)C") == "N-methylethanimine"
