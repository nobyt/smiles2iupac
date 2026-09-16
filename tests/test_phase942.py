"""Phase 942: Polycyclic von Baeyer Hydrocarbons (IUPAC 2013 P-23.2.3 - P-23.2.5).

Systematic nomenclature for saturated 3+-ring fused and bridged hydrocarbons
(tricyclo[...], tetracyclo[...]).
"""
from smiles2iupac import smiles_to_iupac


class TestPhase942PolycyclicVonBaeyer:
    def test_bridged_tricyclic_11_atoms(self):
        assert smiles_to_iupac("C1CC2CC1C1CCCCC21") == "tricyclo[6.2.1.0²,⁷]undecane"

    def test_ortho_fused_tricyclic_15_atoms(self):
        assert smiles_to_iupac("C1CCC2CCCC3CCCCC3C2C1") == "tricyclo[9.4.0.0²,⁷]pentadecane"

    def test_fused_bridged_tricyclic_11_atoms(self):
        assert smiles_to_iupac("C1CC2CCC1C1CCCC21") == "tricyclo[5.2.2.0²,⁶]undecane"
