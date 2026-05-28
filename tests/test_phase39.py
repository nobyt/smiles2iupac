"""
Phase 39 テスト: 炭酸エステル (carbonate esters)

対象 (IUPAC P-65.1.2.3):
  R-O-C(=O)-O-R' パターン（カルボニル炭素に直接 C が付かない）。
  命名規則: {alkyl1} {alkyl2} carbonate (アルファベット順)
  例: CH3-O-C(=O)-O-CH3 → dimethyl carbonate
"""
from smiles2iupac import smiles_to_iupac


class TestSymmetricCarbonates:

    def test_dimethyl_carbonate(self):
        assert smiles_to_iupac("COC(=O)OC") == "dimethyl carbonate"

    def test_diethyl_carbonate(self):
        assert smiles_to_iupac("CCOC(=O)OCC") == "diethyl carbonate"

    def test_dipropyl_carbonate(self):
        assert smiles_to_iupac("CCCOC(=O)OCCC") == "dipropyl carbonate"


class TestAsymmetricCarbonates:

    def test_ethyl_methyl_carbonate(self):
        # アルファベット順: ethyl before methyl
        assert smiles_to_iupac("COC(=O)OCC") == "ethyl methyl carbonate"

    def test_methyl_propyl_carbonate(self):
        assert smiles_to_iupac("COC(=O)OCCC") == "methyl propyl carbonate"


class TestCarbonateVsEster:

    def test_methyl_ethanoate_is_ester(self):
        # ester (C 直結カルボニル) は carbonate と区別される
        assert smiles_to_iupac("CC(=O)OC") == "methyl acetate"
