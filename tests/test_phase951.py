"""Phase 951: von Baeyer secondary bridges with internal (non-terminal) atoms.

Phase 950 made _try_polycyclic_von_baeyer reject any candidate whose secondary
bridge had internal atoms (since the old edge-by-edge collector mis-split such
a bridge into multiple spurious 0-length entries), falling through to the loud
"not supported" ValueError. That was a safe stopgap but left a real class of
tricyclic+ hydrocarbons -- any secondary bridge longer than a single direct
bond -- unsupported.

This phase properly detects secondary bridges of any length: after fixing the
main ring + main bridge numbering, the remaining ring atoms are grouped into
connected components via the leftover edges. Each component must form a
simple path (every atom has total degree, interior + attachment, of exactly 2)
whose two ends attach to two distinct already-numbered atoms; its length is
the number of atoms in the component, and its locants are the two attachment
points. Both zero-length (direct edge) and multi-atom secondary bridges are
now handled uniformly, for tricyclo and beyond (tetracyclo needs 2 secondary
bridges, pentacyclo 3, ...).

All expected names below are OPSIN round-trip verified.
"""
from smiles2iupac import smiles_to_iupac


class TestPhase951InternalAtomSecondaryBridges:
    def test_length_one_secondary_bridge(self):
        assert smiles_to_iupac("C1CC2CC3CC1CC(C2)C3") == "tricyclo[4.3.1.1³,⁸]undecane"

    def test_length_two_secondary_bridge(self):
        assert smiles_to_iupac("C1CC2CC3CC1CCC(C2)C3") == "tricyclo[5.3.1.1⁴,⁹]dodecane"


class TestPhase951MultipleSecondaryBridges:
    def test_pentacyclic_three_secondary_bridges(self):
        assert smiles_to_iupac("C1C2C3C4C1C1C2C3C41") == "pentacyclo[5.2.0.0²,⁵.0³,⁹.0⁴,⁸]nonane"


class TestPhase951ZeroLengthBridgesUnaffected:
    def test_bridged_tricyclic_still_works(self):
        assert smiles_to_iupac("C1CC2CCC1C1CCCCC21") == "tricyclo[6.2.2.0²,⁷]dodecane"

    def test_ortho_fused_tricyclic_still_works(self):
        assert smiles_to_iupac("C1CCC2CCCC3CCCCC3C2C1") == "tricyclo[9.4.0.0²,⁷]pentadecane"

    def test_tetracyclic_still_works(self):
        assert smiles_to_iupac("C1CC2CCC3CCC4CCC1C2C34") == "tetracyclo[8.2.2.0⁴,¹².0⁷,¹¹]tetradecane"

    def test_bicyclo_octane_still_works(self):
        assert smiles_to_iupac("C1CC2CCC1CC2") == "bicyclo[2.2.2]octane"

    def test_norbornane_still_works(self):
        assert smiles_to_iupac("C1CC2CCC1C2") == "bicyclo[2.2.1]heptane"


class TestPhase951SpiroSystemsStillRejected:
    def test_spiro_fused_system_still_raises_not_supported(self):
        import pytest
        with pytest.raises(ValueError, match="not supported"):
            smiles_to_iupac("C1CC2(CCC1)CC1CCCCC1C2")
