"""Phase 931: the 9 "no leading indicated-H" hydro-prefix rows left over
from Phase 930's 106 confirmed bugs -- these ring systems use
"x,y-dihydro"/"x,y,z,w-tetrahydro" prefixes instead of an "NH-" indicated-H
prefix, so Phase 930's regex couldn't derive their locants automatically.
Locants instead come from sibling N-substituted molecules already verified
elsewhere in the suite (test_phase619/627/628/629/726 established the ring
numbering via C-substituted analogues; this phase fills in the N-locant
these files never happened to test), each re-validated with the same
methyl-probe technique as Phase 930.
"""
from smiles2iupac import smiles_to_iupac


class TestPhase931HydroPrefixIndicatedNFix:
    def test_perimidine_n3_methyl(self):
        assert smiles_to_iupac('CN1C=Nc2cccc3cccc1c23') == '3-methylperimidine'

    def test_tetrahydroisoquinoline_n2_methyl(self):
        assert smiles_to_iupac('CN1CCc2ccccc2C1') == '2-methyl-1,2,3,4-tetrahydroisoquinoline'

    def test_tetrahydroquinoxaline_n1_methyl(self):
        assert smiles_to_iupac('CN1CCNc2ccccc21') == '1-methyl-1,2,3,4-tetrahydroquinoxaline'

    def test_benzoxazine_n4_methyl(self):
        assert smiles_to_iupac('CN1CCOc2ccccc21') == '4-methyl-3,4-dihydro-2H-1,4-benzoxazine'

    def test_benzothiazine_n4_methyl(self):
        assert smiles_to_iupac('CN1CCSc2ccccc21') == '4-methyl-3,4-dihydro-2H-1,4-benzothiazine'

    def test_dihydroquinoline_n1_methyl(self):
        assert smiles_to_iupac('CN1CC=Cc2ccccc21') == '1-methyl-1,2-dihydroquinoline'

    def test_dihydroquinoxalinone_n4_methyl(self):
        assert smiles_to_iupac('CN1CC(=O)Nc2ccccc21') == '4-methyl-3,4-dihydroquinoxalin-2(1H)-one'

    def test_dihydroquinoxalinone_n1_methyl(self):
        assert smiles_to_iupac('CN1C(=O)CNc2ccccc21') == '1-methyl-3,4-dihydroquinoxalin-2-one'
