"""Phase 947: cage retained-name disambiguation by ring-skeleton canonicalization.

_try_cage_retained identified adamantane (and cubane) purely by degree sequence
(10 atoms, 4x deg-3 + 6x deg-2), which degenerate C10 cages such as the
twistane-like tricyclo[5.2.1.0^4,10]decane share -- so those were mis-named
"adamantane".  The matcher now confirms the ring skeleton canonicalizes to the
real adamantane/cubane before returning the retained name; otherwise it falls
through to von Baeyer nomenclature.

All expected names below are OPSIN round-trip verified.
"""
import pytest
from smiles2iupac import smiles_to_iupac


class TestPhase947CageDisambiguation:
    def test_degenerate_c10_cage_not_adamantane(self):
        # same degree sequence as adamantane but a different cage topology
        assert smiles_to_iupac("C1CC2CCC3CCC1C23") == "tricyclo[5.2.1.0⁴,¹⁰]decane"

    def test_real_adamantane_still_named(self):
        assert smiles_to_iupac("C12CC3CC(C1)CC(C2)C3") == "adamantane"

    def test_real_adamantane_alt_smiles(self):
        assert smiles_to_iupac("C1C2CC3CC1CC(C2)C3") == "adamantane"

    def test_cubane_still_named(self):
        assert smiles_to_iupac("C12C3C4C1C5C2C3C45") == "cubane"

    def test_substituted_adamantane_still_named(self):
        assert smiles_to_iupac("ClC12CC3CC(C1)CC(C2)C3") == "1-chloroadamantane"
