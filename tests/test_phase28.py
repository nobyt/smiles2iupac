"""
Phase 28 テスト: フェニルアルキルアルコール / フェニルアルキルアミン

対象 (IUPAC P-63.1.1):
  OH が芳香環外の炭素に付いている場合は非環状命名パスで命名する。
  芳香環は "phenyl" 置換基として扱う。

  phenylmethanol: Ph-CH₂-OH
  1-phenylethan-1-ol: Ph-CH(OH)-CH₃
  2-phenylethan-1-ol: Ph-CH₂-CH₂-OH
"""
from smiles2iupac import smiles_to_iupac


class TestPhenylAlkylAlcohols:

    def test_phenylmethanol(self):
        # PhCH₂OH (benzyl alcohol)
        assert smiles_to_iupac("OCc1ccccc1") == "phenylmethanol"

    def test_1_phenylethan_1_ol(self):
        # Ph-CH(OH)-CH₃ → 1-phenylethanol (ethanol 保留名ベース)
        assert smiles_to_iupac("OC(C)c1ccccc1") == "1-phenylethanol"

    def test_2_phenylethan_1_ol(self):
        # Ph-CH₂-CH₂-OH → 2-phenylethanol (ethanol 保留名ベース)
        assert smiles_to_iupac("OCCc1ccccc1") == "2-phenylethanol"

    def test_3_phenylpropan_1_ol(self):
        # Ph-CH₂-CH₂-CH₂-OH
        assert smiles_to_iupac("OCCCc1ccccc1") == "3-phenylpropan-1-ol"


class TestLocantOmissionForOneCarbonChain:

    def test_methanol(self):
        # CH₃OH: 1炭素鎖のためロカント省略
        assert smiles_to_iupac("CO") == "methanol"

    def test_methanamine(self):
        # CH₃NH₂: 1炭素鎖のためロカント省略
        assert smiles_to_iupac("CN") == "methanamine"

    def test_methanethiol(self):
        # CH₃SH: 1炭素鎖のためロカント省略
        assert smiles_to_iupac("CS") == "methanethiol"

    def test_phenylmethanol_no_locant(self):
        # PhCH₂OH: 1炭素主鎖 → ロカント省略
        assert smiles_to_iupac("OCc1ccccc1") == "phenylmethanol"
