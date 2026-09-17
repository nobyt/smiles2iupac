"""Phase 950: reject von Baeyer candidates with the wrong secondary-bridge count.

_try_polycyclic_von_baeyer's candidate scorer computed
`len(sec_bridges) != n_rings - 2` (a tricyclo must have exactly 1 secondary
bridge, a tetracyclo exactly 2, ...) but the mismatch branch was a no-op
(`pass`) -- it never rejected the candidate. A secondary bridge with one or
more internal (non-terminal) atoms gets its edges enumerated individually by
the "direct edge" collector, producing multiple spurious 0-length secondary
bridge entries that share a locant instead of one bridge entry with the
correct bridge length. This produced invalid, OPSIN-unparseable descriptors
such as "tricyclo[4.3.1.0^3,11.0^8,11]undecane" (2 secondary-bridge terms
under a "tricyclo" prefix, which permits exactly 1).

The scorer now rejects (`continue`) any candidate whose secondary-bridge count
doesn't match n_rings - 2, so such systems fall through to the existing loud
"not supported" ValueError guard instead of returning a confidently wrong
name. Genuine von Baeyer systems (whose true secondary bridge happens to have
zero internal atoms) are unaffected.
"""
import pytest
from smiles2iupac import smiles_to_iupac


class TestPhase950WrongBridgeCountRejected:
    def test_length_one_secondary_bridge_raises_not_supported(self):
        with pytest.raises(ValueError, match="not supported"):
            smiles_to_iupac("C1CC2CC3CC1CC(C2)C3")


class TestPhase950GenuineSystemsUnaffected:
    def test_bridged_tricyclic_still_works(self):
        assert smiles_to_iupac("C1CC2CCC1C1CCCCC21") == "tricyclo[6.2.2.0²,⁷]dodecane"

    def test_ortho_fused_tricyclic_still_works(self):
        assert smiles_to_iupac("C1CCC2CCCC3CCCCC3C2C1") == "tricyclo[9.4.0.0²,⁷]pentadecane"

    def test_tetracyclic_still_works(self):
        assert smiles_to_iupac("C1CC2CCC3CCC4CCC1C2C34") == "tetracyclo[8.2.2.0⁴,¹².0⁷,¹¹]tetradecane"

    def test_bicyclo_octane_still_works(self):
        assert smiles_to_iupac("C1CC2CCC1CC2") == "bicyclo[2.2.2]octane"
