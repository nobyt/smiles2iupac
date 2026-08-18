"""Phase 933: isotopic labeling (IUPAC 2013 P-82), Track B's finding #1.

Before this phase, isotope-labeled atoms ([2H], [13C], etc.) were silently
and completely dropped from the output name -- confirmed via
docs/nomenclature_gap_audit.md's B2 probing:
`smiles_to_iupac("[2H]C([2H])([2H])C(=O)O")` returned plain "acetic acid",
with no indication the input was deuterium-labeled at all.

Scope: only isotope labels on principal-chain atoms (the acyclic
`_name_acyclic` pathway, reusing `find_principal_chain`'s locant_map).
Isotope labels on ring atoms or on substituent branches are out of scope
for this phase and are silently ignored, same as before -- a genuinely
new descriptor class was judged too large to implement with full ring
coverage in one phase; this establishes the mechanism and the acyclic
case, which is also the overwhelmingly common real-world case for
isotope-labeled reagents/solvents (CDCl3-shaped, D-labeled amino acids,
13C-labeled simple metabolites, etc.).

Format and prefix ordering verified against OPSIN 2.9.0 (this project's
usual round-trip oracle): isotope descriptor goes AFTER ordinary
substituent prefixes and BEFORE the parent name/suffix (not at the very
front of the whole name) -- "2-chloro(2,2-2H2)acetic acid" parses
correctly in OPSIN; "(2,2-2H2)2-chloroacetic acid" does not ("Could not
find the atom with locant 2"). Stereo descriptors still come before
everything else: "(R)-(1-2H)ethanamine".
"""
from smiles2iupac import smiles_to_iupac


class TestPhase933IsotopeLabeling:
    def test_deuterated_acetic_acid(self):
        assert smiles_to_iupac("[2H]C([2H])([2H])C(=O)O") == "(2,2,2-2H3)acetic acid"

    def test_carbon13_propanoic_acid(self):
        assert smiles_to_iupac("C[13CH2]C(=O)O") == "(2-13C)propanoic acid"

    def test_single_deuterium_ethane(self):
        assert smiles_to_iupac("CC[2H]") == "(2-2H)ethane"

    def test_deuterium_on_suffix_bearing_carbon(self):
        assert smiles_to_iupac("OC[2H]") == "(1-2H)methanol"

    def test_stereo_and_isotope_combined(self):
        assert smiles_to_iupac("C[C@@H]([2H])N") == "(R)-(1-2H)ethanamine"

    def test_isotope_and_ordinary_substituent_combined(self):
        assert smiles_to_iupac("[2H]C([2H])(Cl)C(=O)O") == "2-chloro(2,2-2H2)acetic acid"

    def test_unlabeled_ethanol_unaffected(self):
        assert smiles_to_iupac("OCC") == "ethanol"

    def test_unlabeled_acetic_acid_unaffected(self):
        assert smiles_to_iupac("CC(=O)O") == "acetic acid"
