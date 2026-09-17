"""Phase 948: reject spiro (cut-vertex) systems in von Baeyer nomenclature.

A pure von Baeyer fused/bridged polycycle is biconnected (its ring subgraph has
no articulation point).  A spiro junction is a cut vertex.  Feeding a spiro+fused
system (e.g. a cyclohexane spiro-fused to a fused bicyclic) to
_try_polycyclic_von_baeyer previously produced a degenerate, OPSIN-unparseable
descriptor such as "tricyclo[4.3.0.0^8,11.0^8,12.0^10,11.0^10,14...]tetradecane"
(7 secondary bridges under a "tricyclo" prefix, which permits exactly 1).

The handler now detects the cut vertex and bails, so such systems raise the loud
"not supported" guard instead of returning a confidently wrong name.  Genuine
von Baeyer and simple spiro systems are unaffected.
"""
import pytest
from smiles2iupac import smiles_to_iupac


class TestPhase948SpiroFusedRejected:
    def test_spiro_fused_system_raises_not_supported(self):
        with pytest.raises(ValueError, match="not supported"):
            smiles_to_iupac("C1CC2(CCC1)CC1CCCCC1C2")


class TestPhase948GenuineSystemsUnaffected:
    def test_bridged_tricyclic_still_works(self):
        assert smiles_to_iupac("C1CC2CCC1C1CCCCC21") == "tricyclo[6.2.2.0²,⁷]dodecane"

    def test_ortho_fused_tricyclic_still_works(self):
        assert smiles_to_iupac("C1CCC2CCCC3CCCCC3C2C1") == "tricyclo[9.4.0.0²,⁷]pentadecane"

    def test_tetracyclic_still_works(self):
        assert smiles_to_iupac("C1CC2CCC3CCC4CCC1C2C34") == "tetracyclo[8.2.2.0⁴,¹².0⁷,¹¹]tetradecane"

    def test_simple_spiro_still_works(self):
        assert smiles_to_iupac("C1CC2(CC1)CCCCC2") == "spiro[4.5]decane"
