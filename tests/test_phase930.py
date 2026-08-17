"""Phase 930: fused-ring "indicated-H nitrogen" substituent-drop bug.

Found via systematic adversarial probing of every heteroatom `None` locant
in `_FUSED_LOCANT_MAP` (docs/locant_map_audit.md): for each heteroatom that
still carries an implicit H (the only atoms a substituent could plausibly
land on), attach a methyl and check whether it survives into the assembled
name. 106 of 108 such positions silently dropped the substituent entirely
(heterocycle_handler.py:5024's `if loc is None: continue`); one case even
collapsed to the unrelated name "benzene".

This phase fixes the 63 of those 106 where the correct locant is simply the
already-known indicated-H locant baked into the ring system's own retained
name (e.g. "1H-pyrazolo[3,4-b]pyridine" -> locant 1), and where the fix was
verified not to disturb substructure-match selection for that entry. The
remaining ~43 (symmetric/tautomeric ring systems where multiple automorphic
substructure matches compete, plus a `_FUSED_HETERO_RETAINED` naming bug
that gives one entry the wrong name entirely) are deferred to a follow-up
phase -- see docs/locant_map_audit.md.
"""
from smiles2iupac import smiles_to_iupac


class TestPhase930IndicatedHSubstituentFix:
    def test_1h_1_2_3_triazole_n_methyl(self):
        """1H-1,2,3-triazole -> 1-methyl-1,2,3-triazole"""
        assert smiles_to_iupac('Cn1ccnn1') == '1-methyl-1,2,3-triazole'

    def test_1h_1_2_4_triazole_n_methyl(self):
        """1H-1,2,4-triazole -> 1-methyl-1,2,4-triazole"""
        assert smiles_to_iupac('Cn1cncn1') == '1-methyl-1,2,4-triazole'

    def test_1h_tetrazole_n_methyl(self):
        """1H-tetrazole -> 1-methyltetrazole"""
        assert smiles_to_iupac('Cn1cnnn1') == '1-methyltetrazole'

    def test_1h_1_2_3_triazolo_4_5_b_pyridine_n_methyl(self):
        """1H-[1,2,3]triazolo[4,5-b]pyridine -> 1-methyl[1,2,3]triazolo[4,5-b]pyridine"""
        assert smiles_to_iupac('Cn1nnc2cccnc21') == '1-methyl[1,2,3]triazolo[4,5-b]pyridine'

    def test_1h_1_2_3_triazolo_4_5_c_pyridazine_n_methyl(self):
        """1H-[1,2,3]triazolo[4,5-c]pyridazine -> 1-methyl[1,2,3]triazolo[4,5-c]pyridazine"""
        assert smiles_to_iupac('Cn1nnc2nnccc21') == '1-methyl[1,2,3]triazolo[4,5-c]pyridazine'

    def test_1h_1_2_3_triazolo_4_5_d_pyrimidine_n_methyl(self):
        """1H-[1,2,3]triazolo[4,5-d]pyrimidine -> 1-methyl[1,2,3]triazolo[4,5-d]pyrimidine"""
        assert smiles_to_iupac('Cn1nnc2ncncc21') == '1-methyl[1,2,3]triazolo[4,5-d]pyrimidine'

    def test_1h_1_2_3_triazolo_4_5_e_1_2_4_triazine_n_methyl(self):
        """1H-[1,2,3]triazolo[4,5-e][1,2,4]triazine -> 1-methyl[1,2,3]triazolo[4,5-e][1,2,4]triazine"""
        assert smiles_to_iupac('Cn1nnc2ncnnc21') == '1-methyl[1,2,3]triazolo[4,5-e][1,2,4]triazine'

    def test_1h_1_2_3_triazolo_4_5_e_pyrazine_n_methyl(self):
        """1H-[1,2,3]triazolo[4,5-e]pyrazine -> 1-methyl[1,2,3]triazolo[4,5-e]pyrazine"""
        assert smiles_to_iupac('Cn1nnc2nccnc21') == '1-methyl[1,2,3]triazolo[4,5-e]pyrazine'

    def test_1h_1_2_3_triazolo_5_4_d_1_2_3_triazine_n_methyl(self):
        """1H-[1,2,3]triazolo[5,4-d][1,2,3]triazine -> 1-methyl[1,2,3]triazolo[5,4-d][1,2,3]triazine"""
        assert smiles_to_iupac('Cn1nnc2cnnnc21') == '1-methyl[1,2,3]triazolo[5,4-d][1,2,3]triazine'

    def test_1h_imidazo_4_5_b_pyridine_n_methyl(self):
        """1H-imidazo[4,5-b]pyridine -> 1-methylimidazo[4,5-b]pyridine"""
        assert smiles_to_iupac('Cn1cnc2cccnc21') == '1-methylimidazo[4,5-b]pyridine'

    def test_1h_imidazo_4_5_b_quinoline_n_methyl(self):
        """1H-imidazo[4,5-b]quinoline -> 1-methylimidazo[4,5-b]quinoline"""
        assert smiles_to_iupac('Cn1cnc2cc3ccccc3nc21') == '1-methylimidazo[4,5-b]quinoline'

    def test_5h_imidazo_4_5_c_pyridazine_n_methyl(self):
        """5H-imidazo[4,5-c]pyridazine -> 5-methylimidazo[4,5-c]pyridazine"""
        assert smiles_to_iupac('Cn1cnc2nnccc21') == '5-methylimidazo[4,5-c]pyridazine'

    def test_1h_imidazo_4_5_c_pyridine_n_methyl(self):
        """1H-imidazo[4,5-c]pyridine -> 1-methylimidazo[4,5-c]pyridine"""
        assert smiles_to_iupac('Cn1cnc2cnccc21') == '1-methylimidazo[4,5-c]pyridine'

    def test_1h_imidazo_4_5_d_pyridazine_n_methyl(self):
        """1H-imidazo[4,5-d]pyridazine -> 1-methylimidazo[4,5-d]pyridazine"""
        assert smiles_to_iupac('Cn1cnc2cnncc21') == '1-methylimidazo[4,5-d]pyridazine'

    def test_7h_imidazo_4_5_e_1_2_4_triazine_n_methyl(self):
        """7H-imidazo[4,5-e][1,2,4]triazine -> 7-methylimidazo[4,5-e][1,2,4]triazine"""
        assert smiles_to_iupac('Cn1cnc2ncnnc21') == '7-methylimidazo[4,5-e][1,2,4]triazine'

    def test_1h_pyrazolo_3_4_b_pyridine_n_methyl(self):
        """1H-pyrazolo[3,4-b]pyridine -> 1-methylpyrazolo[3,4-b]pyridine"""
        assert smiles_to_iupac('Cn1ncc2cccnc21') == '1-methylpyrazolo[3,4-b]pyridine'

    def test_1h_pyrazolo_3_4_b_quinoline_n_methyl(self):
        """1H-pyrazolo[3,4-b]quinoline -> 1-methylpyrazolo[3,4-b]quinoline"""
        assert smiles_to_iupac('Cn1ncc2cc3ccccc3nc21') == '1-methylpyrazolo[3,4-b]quinoline'

    def test_1h_pyrazolo_3_4_c_pyridazine_n_methyl(self):
        """1H-pyrazolo[3,4-c]pyridazine -> 1-methylpyrazolo[3,4-c]pyridazine"""
        assert smiles_to_iupac('Cn1cc2ccnnc2n1') == '1-methylpyrazolo[3,4-c]pyridazine'

    def test_1h_pyrazolo_3_4_c_pyridine_n_methyl(self):
        """1H-pyrazolo[3,4-c]pyridine -> 1-methylpyrazolo[3,4-c]pyridine"""
        assert smiles_to_iupac('Cn1cc2ccncc2n1') == '1-methylpyrazolo[3,4-c]pyridine'

    def test_1h_pyrazolo_3_4_d_pyrimidine_n_methyl(self):
        """1H-pyrazolo[3,4-d]pyrimidine -> 1-methylpyrazolo[3,4-d]pyrimidine"""
        assert smiles_to_iupac('Cn1cc2cncnc2n1') == '1-methylpyrazolo[3,4-d]pyrimidine'

    def test_1h_pyrazolo_3_4_e_1_2_4_triazine_n_methyl(self):
        """1H-pyrazolo[3,4-e][1,2,4]triazine -> 1-methylpyrazolo[3,4-e][1,2,4]triazine"""
        assert smiles_to_iupac('Cn1ncnc2nncc1-2') == '1-methylpyrazolo[3,4-e][1,2,4]triazine'

    def test_1h_pyrazolo_4_3_b_pyridine_n_methyl(self):
        """1H-pyrazolo[4,3-b]pyridine -> 1-methylpyrazolo[4,3-b]pyridine"""
        assert smiles_to_iupac('Cn1cc2ncccc2n1') == '1-methylpyrazolo[4,3-b]pyridine'

    def test_1h_pyrazolo_4_3_d_1_2_3_triazine_n_methyl(self):
        """1H-pyrazolo[4,3-d][1,2,3]triazine -> 1-methylpyrazolo[4,3-d][1,2,3]triazine"""
        assert smiles_to_iupac('Cn1nncc2nncc1-2') == '1-methylpyrazolo[4,3-d][1,2,3]triazine'

    def test_1h_pyrazolo_4_5_c_pyridazine_n_methyl(self):
        """1H-pyrazolo[4,5-c]pyridazine -> 1-methylpyrazolo[4,5-c]pyridazine"""
        assert smiles_to_iupac('Cn1ncc2nnccc21') == '1-methylpyrazolo[4,5-c]pyridazine'

    def test_1h_pyrazolo_4_5_c_pyridine_n_methyl(self):
        """1H-pyrazolo[4,5-c]pyridine -> 1-methylpyrazolo[4,5-c]pyridine"""
        assert smiles_to_iupac('Cn1ncc2cnccc21') == '1-methylpyrazolo[4,5-c]pyridine'

    def test_1h_pyrazolo_4_5_d_pyridazine_n_methyl(self):
        """1H-pyrazolo[4,5-d]pyridazine -> 1-methylpyrazolo[4,5-d]pyridazine"""
        assert smiles_to_iupac('Cn1ncc2cnncc21') == '1-methylpyrazolo[4,5-d]pyridazine'

    def test_1h_pyrazolo_4_5_d_pyrimidine_n_methyl(self):
        """1H-pyrazolo[4,5-d]pyrimidine -> 1-methylpyrazolo[4,5-d]pyrimidine"""
        assert smiles_to_iupac('Cn1ncc2ncncc21') == '1-methylpyrazolo[4,5-d]pyrimidine'

    def test_1h_pyrazolo_4_5_e_1_2_4_triazine_n_methyl(self):
        """1H-pyrazolo[4,5-e][1,2,4]triazine -> 1-methylpyrazolo[4,5-e][1,2,4]triazine"""
        assert smiles_to_iupac('Cn1ncc2ncnnc21') == '1-methylpyrazolo[4,5-e][1,2,4]triazine'

    def test_1h_pyrazolo_4_5_e_pyrazine_n_methyl(self):
        """1H-pyrazolo[4,5-e]pyrazine -> 1-methylpyrazolo[4,5-e]pyrazine"""
        assert smiles_to_iupac('Cn1ncc2nccnc21') == '1-methylpyrazolo[4,5-e]pyrazine'

    def test_1h_pyrazolo_5_4_d_1_2_3_triazine_n_methyl(self):
        """1H-pyrazolo[5,4-d][1,2,3]triazine -> 1-methylpyrazolo[5,4-d][1,2,3]triazine"""
        assert smiles_to_iupac('Cn1ncc2cnnnc21') == '1-methylpyrazolo[5,4-d][1,2,3]triazine'

    def test_1h_pyrrolo_2_3_b_pyridine_n_methyl(self):
        """1H-pyrrolo[2,3-b]pyridine -> 1-methylpyrrolo[2,3-b]pyridine"""
        assert smiles_to_iupac('Cn1ccc2cccnc21') == '1-methylpyrrolo[2,3-b]pyridine'

    def test_1h_pyrrolo_2_3_b_quinoline_n_methyl(self):
        """1H-pyrrolo[2,3-b]quinoline -> 1-methylpyrrolo[2,3-b]quinoline"""
        assert smiles_to_iupac('Cn1ccc2cc3ccccc3nc21') == '1-methylpyrrolo[2,3-b]quinoline'

    def test_7h_pyrrolo_2_3_c_pyridazine_n_methyl(self):
        """7H-pyrrolo[2,3-c]pyridazine -> 7-methylpyrrolo[2,3-c]pyridazine"""
        assert smiles_to_iupac('Cn1ccc2ccnnc21') == '7-methylpyrrolo[2,3-c]pyridazine'

    def test_1h_pyrrolo_2_3_c_pyridine_n_methyl(self):
        """1H-pyrrolo[2,3-c]pyridine -> 1-methylpyrrolo[2,3-c]pyridine"""
        assert smiles_to_iupac('Cn1ccc2ccncc21') == '1-methylpyrrolo[2,3-c]pyridine'

    def test_7h_pyrrolo_2_3_d_1_2_3_triazine_n_methyl(self):
        """7H-pyrrolo[2,3-d][1,2,3]triazine -> 7-methylpyrrolo[2,3-d][1,2,3]triazine"""
        assert smiles_to_iupac('Cn1ccc2cnnnc21') == '7-methylpyrrolo[2,3-d][1,2,3]triazine'

    def test_1h_pyrrolo_2_3_d_pyridazine_n_methyl(self):
        """1H-pyrrolo[2,3-d]pyridazine -> 1-methylpyrrolo[2,3-d]pyridazine"""
        assert smiles_to_iupac('Cn1ccc2cnncc21') == '1-methylpyrrolo[2,3-d]pyridazine'

    def test_7h_pyrrolo_2_3_d_pyrimidine_n_methyl(self):
        """7H-pyrrolo[2,3-d]pyrimidine -> 7-methylpyrrolo[2,3-d]pyrimidine"""
        assert smiles_to_iupac('Cn1ccc2cncnc21') == '7-methylpyrrolo[2,3-d]pyrimidine'

    def test_1h_pyrrolo_2_3_e_1_2_4_triazine_n_methyl(self):
        """1H-pyrrolo[2,3-e][1,2,4]triazine -> 1-methylpyrrolo[2,3-e][1,2,4]triazine"""
        assert smiles_to_iupac('Cn1ncnc2nccc1-2') == '1-methylpyrrolo[2,3-e][1,2,4]triazine'

    def test_1h_pyrrolo_2_3_e_pyrazine_n_methyl(self):
        """1H-pyrrolo[2,3-e]pyrazine -> 1-methylpyrrolo[2,3-e]pyrazine"""
        assert smiles_to_iupac('Cn1ccnc2nccc1-2') == '1-methylpyrrolo[2,3-e]pyrazine'

    def test_1h_pyrrolo_3_2_b_pyridine_n_methyl(self):
        """1H-pyrrolo[3,2-b]pyridine -> 1-methylpyrrolo[3,2-b]pyridine"""
        assert smiles_to_iupac('Cn1ccc2ncccc21') == '1-methylpyrrolo[3,2-b]pyridine'

    def test_1h_pyrrolo_3_2_b_quinoline_n_methyl(self):
        """1H-pyrrolo[3,2-b]quinoline -> 1-methylpyrrolo[3,2-b]quinoline"""
        assert smiles_to_iupac('Cn1ccc2nc3ccccc3cc21') == '1-methylpyrrolo[3,2-b]quinoline'

    def test_1h_pyrrolo_3_2_c_pyridazine_n_methyl(self):
        """1H-pyrrolo[3,2-c]pyridazine -> 1-methylpyrrolo[3,2-c]pyridazine"""
        assert smiles_to_iupac('Cn1nccc2nccc1-2') == '1-methylpyrrolo[3,2-c]pyridazine'

    def test_1h_pyrrolo_3_2_c_pyridine_n_methyl(self):
        """1H-pyrrolo[3,2-c]pyridine -> 1-methylpyrrolo[3,2-c]pyridine"""
        assert smiles_to_iupac('Cn1ccc2cnccc21') == '1-methylpyrrolo[3,2-c]pyridine'

    def test_1h_pyrrolo_3_2_d_1_2_3_triazine_n_methyl(self):
        """1H-pyrrolo[3,2-d][1,2,3]triazine -> 1-methylpyrrolo[3,2-d][1,2,3]triazine"""
        assert smiles_to_iupac('Cn1nncc2nccc1-2') == '1-methylpyrrolo[3,2-d][1,2,3]triazine'

    def test_1h_pyrrolo_3_2_d_pyrimidine_n_methyl(self):
        """1H-pyrrolo[3,2-d]pyrimidine -> 1-methylpyrrolo[3,2-d]pyrimidine"""
        assert smiles_to_iupac('Cn1cncc2nccc1-2') == '1-methylpyrrolo[3,2-d]pyrimidine'

    def test_7h_pyrrolo_3_2_e_1_2_4_triazine_n_methyl(self):
        """7H-pyrrolo[3,2-e][1,2,4]triazine -> 7-methylpyrrolo[3,2-e][1,2,4]triazine"""
        assert smiles_to_iupac('Cn1ccc2ncnnc21') == '7-methylpyrrolo[3,2-e][1,2,4]triazine'

    def test_1h_pyrrolo_3_4_b_pyridine_n_methyl(self):
        """1H-pyrrolo[3,4-b]pyridine -> 1-methylpyrrolo[3,4-b]pyridine"""
        assert smiles_to_iupac('Cn1cccc2cncc1-2') == '1-methylpyrrolo[3,4-b]pyridine'

    def test_1h_pyrrolo_3_4_c_pyridazine_n_methyl(self):
        """1H-pyrrolo[3,4-c]pyridazine -> 1-methylpyrrolo[3,4-c]pyridazine"""
        assert smiles_to_iupac('Cn1nccc2cncc1-2') == '1-methylpyrrolo[3,4-c]pyridazine'

    def test_1h_pyrrolo_3_4_d_1_2_3_triazine_n_methyl(self):
        """1H-pyrrolo[3,4-d][1,2,3]triazine -> 1-methylpyrrolo[3,4-d][1,2,3]triazine"""
        assert smiles_to_iupac('Cn1nncc2cncc1-2') == '1-methylpyrrolo[3,4-d][1,2,3]triazine'

    def test_1h_pyrrolo_3_4_d_pyrimidine_n_methyl(self):
        """1H-pyrrolo[3,4-d]pyrimidine -> 1-methylpyrrolo[3,4-d]pyrimidine"""
        assert smiles_to_iupac('Cn1cncc2cncc1-2') == '1-methylpyrrolo[3,4-d]pyrimidine'

    def test_1h_pyrrolo_3_4_e_1_2_4_triazine_n_methyl(self):
        """1H-pyrrolo[3,4-e][1,2,4]triazine -> 1-methylpyrrolo[3,4-e][1,2,4]triazine"""
        assert smiles_to_iupac('Cn1ncnc2cncc1-2') == '1-methylpyrrolo[3,4-e][1,2,4]triazine'

    def test_1h_pyrrolo_3_4_e_pyrazine_n_methyl(self):
        """1H-pyrrolo[3,4-e]pyrazine -> 1-methylpyrrolo[3,4-e]pyrazine"""
        assert smiles_to_iupac('Cn1ccnc2cncc1-2') == '1-methylpyrrolo[3,4-e]pyrazine'

    def test_2h_1_2_3_triazolo_4_5_c_pyridine_n_methyl(self):
        """2H-[1,2,3]triazolo[4,5-c]pyridine -> 2-methyl[1,2,3]triazolo[4,5-c]pyridine"""
        assert smiles_to_iupac('Cn1nc2ccncc2n1') == '2-methyl[1,2,3]triazolo[4,5-c]pyridine'

    def test_7h_purine_n_methyl(self):
        """7H-purine -> 7-methylpurine"""
        assert smiles_to_iupac('Cn1cnc2ncncc21') == '7-methylpurine'

    def test_1h_naphtho_2_3_d_imidazole_n_methyl(self):
        """1H-naphtho[2,3-d]imidazole -> 1-methylnaphtho[2,3-d]imidazole"""
        assert smiles_to_iupac('Cn1cnc2cc3ccccc3cc21') == '1-methylnaphtho[2,3-d]imidazole'

    def test_1h_naphtho_2_3_d_pyrazole_n_methyl(self):
        """1H-naphtho[2,3-d]pyrazole -> 1-methylnaphtho[2,3-d]pyrazole"""
        assert smiles_to_iupac('Cn1ncc2cc3ccccc3cc21') == '1-methylnaphtho[2,3-d]pyrazole'

    def test_1h_naphtho_2_1_d_imidazole_n_methyl(self):
        """1H-naphtho[2,1-d]imidazole -> 1-methylnaphtho[2,1-d]imidazole"""
        assert smiles_to_iupac('Cn1cnc2c3ccccc3ccc21') == '1-methylnaphtho[2,1-d]imidazole'

    def test_1h_naphtho_2_1_d_pyrazole_n_methyl(self):
        """1H-naphtho[2,1-d]pyrazole -> 1-methylnaphtho[2,1-d]pyrazole"""
        assert smiles_to_iupac('Cn1ncc2ccc3ccccc3c21') == '1-methylnaphtho[2,1-d]pyrazole'

    def test_1h_naphtho_1_2_d_1_2_3_triazole_n_methyl(self):
        """1H-naphtho[1,2-d][1,2,3]triazole -> 1-methylnaphtho[1,2-d][1,2,3]triazole"""
        assert smiles_to_iupac('Cn1nnc2ccc3ccccc3c21') == '1-methylnaphtho[1,2-d][1,2,3]triazole'

    def test_1h_benzo_f_indole_n_methyl(self):
        """1H-benzo[f]indole -> 1-methylbenzo[f]indole"""
        assert smiles_to_iupac('Cn1ccc2cc3ccccc3cc21') == '1-methylbenzo[f]indole'

    def test_1h_benzo_g_indole_n_methyl(self):
        """1H-benzo[g]indole -> 1-methylbenzo[g]indole"""
        assert smiles_to_iupac('Cn1ccc2ccc3ccccc3c21') == '1-methylbenzo[g]indole'

    def test_9h_pyrido_2_3_b_indole_n_methyl(self):
        """9H-pyrido[2,3-b]indole -> 9-methylpyrido[2,3-b]indole"""
        assert smiles_to_iupac('Cn1c2ccccc2c2cccnc21') == '9-methylpyrido[2,3-b]indole'

    def test_9h_pyrido_3_4_b_indole_n_methyl(self):
        """9H-pyrido[3,4-b]indole -> 9-methylpyrido[3,4-b]indole"""
        assert smiles_to_iupac('Cn1c2ccccc2c2ccncc21') == '9-methylpyrido[3,4-b]indole'
