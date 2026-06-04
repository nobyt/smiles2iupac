"""
name_assembler のユニットテスト (RDKit 不要)。
IUPAC 命名ロジックの核心部分を単独でテストする。
"""

from smiles2iupac.name_assembler import assemble_name, _build_prefix


class TestBuildPrefix:

    def test_single_substituent(self):
        assert _build_prefix([(2, "methyl")]) == "2-methyl"

    def test_two_different_substituents_alpha_order(self):
        # chloro (c) comes before methyl (m)
        result = _build_prefix([(4, "chloro"), (2, "methyl")])
        assert result == "4-chloro-2-methyl"

    def test_two_same_substituents(self):
        result = _build_prefix([(2, "methyl"), (3, "methyl")])
        assert result == "2,3-dimethyl"

    def test_same_locant_twice(self):
        result = _build_prefix([(2, "methyl"), (2, "methyl")])
        assert result == "2,2-dimethyl"

    def test_no_substituents(self):
        assert _build_prefix([]) == ""

    def test_di_prefix_ignored_in_sort(self):
        # di-methyl vs chloro: 'm' vs 'c' → chloro first
        result = _build_prefix([(4, "chloro"), (2, "methyl"), (3, "methyl")])
        assert result == "4-chloro-2,3-dimethyl"


class TestAssembleName:

    # ─── アルカン ──────────────────────────────────────────────────

    def test_methane(self):
        assert assemble_name(1, "alkane", {"ene": [], "yne": []}, [], []) == "methane"

    def test_ethane(self):
        assert assemble_name(2, "alkane", {"ene": [], "yne": []}, [], []) == "ethane"

    def test_hexane(self):
        assert assemble_name(6, "alkane", {"ene": [], "yne": []}, [], []) == "hexane"

    def test_2_methylpropane(self):
        result = assemble_name(3, "alkane", {"ene": [], "yne": []}, [(2, "methyl")], [])
        assert result == "2-methylpropane"

    def test_3_methylpentane(self):
        result = assemble_name(5, "alkane", {"ene": [], "yne": []}, [(3, "methyl")], [])
        assert result == "3-methylpentane"

    def test_4_chloro_2_methylhexane(self):
        result = assemble_name(
            6, "alkane", {"ene": [], "yne": []},
            [(4, "chloro"), (2, "methyl")], []
        )
        assert result == "4-chloro-2-methylhexane"

    # ─── アルコール ────────────────────────────────────────────────

    def test_ethan_1_ol(self):
        # IUPAC 2013 保留名 (PIN): ethanol
        result = assemble_name(2, "alcohol", {"ene": [], "yne": []}, [], [], suffix_locant=1)
        assert result == "ethanol"

    def test_propan_2_ol(self):
        result = assemble_name(3, "alcohol", {"ene": [], "yne": []}, [], [], suffix_locant=2)
        assert result == "propan-2-ol"

    def test_pentan_2_ol(self):
        result = assemble_name(5, "alcohol", {"ene": [], "yne": []}, [], [], suffix_locant=2)
        assert result == "pentan-2-ol"

    # ─── アルデヒド ────────────────────────────────────────────────

    def test_ethanal(self):
        result = assemble_name(2, "aldehyde", {"ene": [], "yne": []}, [], [])
        assert result == "ethanal"

    def test_propanal(self):
        result = assemble_name(3, "aldehyde", {"ene": [], "yne": []}, [], [])
        assert result == "propanal"

    def test_butanal(self):
        result = assemble_name(4, "aldehyde", {"ene": [], "yne": []}, [], [])
        assert result == "butanal"

    # ─── カルボン酸 ────────────────────────────────────────────────

    def test_ethanoic_acid(self):
        # IUPAC 2013 保留名 (PIN): acetic acid
        result = assemble_name(2, "carboxylic_acid", {"ene": [], "yne": []}, [], [])
        assert result == "acetic acid"

    def test_propanoic_acid(self):
        result = assemble_name(3, "carboxylic_acid", {"ene": [], "yne": []}, [], [])
        assert result == "propanoic acid"

    def test_butanoic_acid(self):
        result = assemble_name(4, "carboxylic_acid", {"ene": [], "yne": []}, [], [])
        assert result == "butanoic acid"

    # ─── ケトン ────────────────────────────────────────────────────

    def test_propan_2_one(self):
        result = assemble_name(3, "ketone", {"ene": [], "yne": []}, [], [], suffix_locant=2)
        assert result == "propan-2-one"

    def test_butan_2_one(self):
        result = assemble_name(4, "ketone", {"ene": [], "yne": []}, [], [], suffix_locant=2)
        assert result == "butan-2-one"

    def test_pentan_3_one(self):
        result = assemble_name(5, "ketone", {"ene": [], "yne": []}, [], [], suffix_locant=3)
        assert result == "pentan-3-one"

    # ─── アルケン ──────────────────────────────────────────────────

    def test_ethene(self):
        # 2C: ロカント省略
        result = assemble_name(2, "alkene", {"ene": [1], "yne": []}, [], [])
        assert result == "ethene"

    def test_propene(self):
        result = assemble_name(3, "alkene", {"ene": [1], "yne": []}, [], [])
        assert result == "prop-1-ene"

    def test_but_2_ene(self):
        result = assemble_name(4, "alkene", {"ene": [2], "yne": []}, [], [])
        assert result == "but-2-ene"

    def test_hex_1_ene(self):
        result = assemble_name(6, "alkene", {"ene": [1], "yne": []}, [], [])
        assert result == "hex-1-ene"

    # ─── アルキン ──────────────────────────────────────────────────

    def test_ethyne(self):
        # 2C: ロカント省略
        result = assemble_name(2, "alkyne", {"ene": [], "yne": [1]}, [], [])
        assert result == "ethyne"

    def test_but_2_yne(self):
        result = assemble_name(4, "alkyne", {"ene": [], "yne": [2]}, [], [])
        assert result == "but-2-yne"

    def test_propyne(self):
        result = assemble_name(3, "alkyne", {"ene": [], "yne": [1]}, [], [])
        assert result == "prop-1-yne"
