"""Phase 941: Acyclic Skeletal Replacement ('a') Nomenclature (IUPAC 2013 P-15.4, P-55).

Linear polyethers, polyamines, polythioethers containing 4 or more heteroatoms
are named using skeletal replacement ('a') nomenclature as Preferred IUPAC Names (PIN).
"""
from smiles2iupac import smiles_to_iupac


class TestPhase941AcyclicSkeletalReplacement:
    def test_tetraoxadodecane(self):
        # 12-atom chain: C-O-C-C-O-C-C-O-C-C-O-C
        assert smiles_to_iupac("COCCOCCOCCOC") == "2,5,8,11-tetraoxadodecane"

    def test_tetraoxatetradecane(self):
        # 14-atom chain: C-C-O-C-C-O-C-C-O-C-C-O-C-C
        assert smiles_to_iupac("CCOCCOCCOCCOCC") == "3,6,9,12-tetraoxatetradecane"

    def test_pentaoxapentadecane(self):
        # 15-atom chain
        assert smiles_to_iupac("COCCOCCOCCOCCOC") == "2,5,8,11,14-pentaoxapentadecane"

    def test_pentaazatridecane(self):
        # 13-atom chain: N-C-C-N-C-C-N-C-C-N-C-C-N
        assert smiles_to_iupac("NCCNCCNCCNCCN") == "1,4,7,10,13-pentaazatridecane"

    def test_tetrathiadodecane(self):
        # 12-atom chain: C-S-C-C-S-C-C-S-C-C-S-C
        assert smiles_to_iupac("CSCCSCCSCCSC") == "2,5,8,11-tetrathiadodecane"
