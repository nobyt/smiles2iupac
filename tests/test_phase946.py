"""Phase 946: hyphen separator between mixed skeletal-replacement prefixes (P-15.4/P-55).

Phase 941 joined the per-element 'a' prefixes with "" so a chain containing
more than one kind of heteroatom produced a malformed, run-together name such
as "2,5,11-trioxa8-azadodecane" (no separator before the "8-aza" locant).
The prefixes must be joined so a hyphen separates a prefix from the following
locant: "2,5,11-trioxa-8-azadodecane".

All expected names below are OPSIN round-trip verified.
"""
from smiles2iupac import smiles_to_iupac


class TestPhase946MixedHeteroatomChains:
    def test_trioxa_aza(self):
        assert smiles_to_iupac("COCCOCCNCCOC") == "2,5,11-trioxa-8-azadodecane"

    def test_dioxa_dithia(self):
        assert smiles_to_iupac("CSCCOCCOCCSC") == "5,8-dioxa-2,11-dithiadodecane"

    def test_dioxa_diaza(self):
        assert smiles_to_iupac("NCCOCCOCCNCC") == "4,7-dioxa-1,10-diazadodecane"

    def test_oxa_thia_aza(self):
        assert smiles_to_iupac("COCCSCCNCCOC") == "2,11-dioxa-5-thia-8-azadodecane"


class TestPhase946SingleHeteroatomUnaffected:
    def test_tetraoxa_no_spurious_hyphen(self):
        assert smiles_to_iupac("COCCOCCOCCOC") == "2,5,8,11-tetraoxadodecane"

    def test_tetraaza_no_spurious_hyphen(self):
        assert smiles_to_iupac("NCCNCCNCCN") == "1,4,7,10-tetraazadecane"
