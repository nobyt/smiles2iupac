"""Phase 934: fixed the remaining 34 fused-ring "indicated-H nitrogen"
substituent-drop bugs left over from Phase 930/931 -- these needed a
structurally different fix, not just a locant-table patch.

Root cause #1 (tautomer disambiguation): `_try_fused_hetero_retained()`'s
fallback for N-substituted fused heterocycles (heterocycle_handler.py
~4817-4874) used to extract the ring-only fragment SMILES text and, when
that failed to parse standalone (because the substituent replaced the
ring N-H), try converting each remaining bare aromatic 'n' into '[nH]' in
TEXT order, taking the first candidate that happened to match some
_FUSED_HETERO_RETAINED entry. When two DIFFERENT tautomeric SMILES
representations of the same underlying ring system (confirmed via
matching InChI) both resolve to the identical retained name, this could
pick the WRONG one -- structurally correct as a NAME, but keyed to a
_FUSED_LOCANT_MAP entry that doesn't correspond to where the substituent
actually sits, silently dropping it (or, in one case, producing a
completely unrelated name, "benzene").

Fixed by replacing the ambiguous text-position guessing with a direct,
unambiguous reconstruction: identify the exact ring nitrogen atom that
lost its H to the substituent (there's only one, in the cases this phase
covers), explicitly set an H back on THAT specific atom, remove the
substituent atoms, and canonicalize -- eliminating the tautomer-choice
ambiguity by construction rather than guessing among candidates.

Root cause #2 (indicated-H over-retention): with root cause #1 fixed, 14
of the 34 still failed OPSIN round-trip verification -- the *locant*
assigned to the newly-correctly-identified core was itself wrong for
about half of these, traced to a second, independent issue: the "leading
digit" of a retained name like "1H-tetrazole" reflects the IUPAC
convention of numbering the UNSUBSTITUTED parent to give indicated
hydrogen the lowest possible locant -- it does NOT necessarily correspond
to which physical atom bears the H in any one particular SMILES text
representation of that same parent (verified via InChI: multiple SMILES
texts, with H written at different ring positions, can be genuine
tautomers of the identical named parent). The correct locant for each of
these 34 substitution sites was instead derived empirically: for each
candidate locant N, check whether "N-methyl{ring name}" round-trips via
OPSIN to the same structure as the actual substituted molecule, and use
whichever N does. Once that produces a locant that legitimately differs
from the retained name's own leading digit, the "1H-"/"3H-" prefix
inherited unchanged from the retained-name lookup becomes actively WRONG
(OPSIN then tries to satisfy an H at position 1 AND a methyl at position
2/3 simultaneously, producing a different, non-aromatic structure) --
fixed by having _try_fused_hetero_retained() unconditionally drop the
leading indicated-H prefix whenever the direct-reconstruction path (root
cause #1's fix) identified the substituted atom, regardless of whether
its locant numerically matches the prefix's own digit -- the direct
reconstruction already proves that atom is where the parent's
tautomeric ambiguity lived, so once it is substituted the ambiguity is
gone.

All 34 fixed names verified via `java -jar opsin-cli-2.9.0-jar-with-dependencies.jar`
round-trip to the same canonical structure as the SMILES tested.
"""
from smiles2iupac import smiles_to_iupac


class TestPhase934TautomerDisambiguationFix:
    def test_2_methyltetrazole(self):
        assert smiles_to_iupac('Cn1ncnn1') == '2-methyltetrazole'

    def test_2_methyl_1_2_3_triazole(self):
        assert smiles_to_iupac('Cn1nccn1') == '2-methyl-1,2,3-triazole'

    def test_1_methyl_1_2_3_triazolo_4_5_b_pyridine(self):
        assert smiles_to_iupac('Cn1nnc2ncccc21') == '1-methyl[1,2,3]triazolo[4,5-b]pyridine'

    def test_1_methyl_1_2_3_triazolo_4_5_d_1_2_3_triazine(self):
        assert smiles_to_iupac('Cn1nnc2nnncc21') == '1-methyl[1,2,3]triazolo[4,5-d][1,2,3]triazine'

    def test_3_methyl_1_2_3_triazolo_5_4_c_pyridazine(self):
        assert smiles_to_iupac('Cn1nnc2ccnnc21') == '3-methyl[1,2,3]triazolo[5,4-c]pyridazine'

    def test_3_methyl_1_2_3_triazolo_5_4_c_pyridine(self):
        assert smiles_to_iupac('Cn1nnc2ccncc21') == '3-methyl[1,2,3]triazolo[5,4-c]pyridine'

    def test_3_methyl_1_2_3_triazolo_5_4_d_pyrimidine(self):
        assert smiles_to_iupac('Cn1nnc2cncnc21') == '3-methyl[1,2,3]triazolo[5,4-d]pyrimidine'

    def test_1_methyl_1_2_3_triazolo_5_4_e_1_2_4_triazine(self):
        assert smiles_to_iupac('Cn1nnc2nncnc21') == '1-methyl[1,2,3]triazolo[5,4-e][1,2,4]triazine'

    def test_3_methylimidazo_4_5_c_pyridine(self):
        assert smiles_to_iupac('Cn1cnc2ccncc21') == '3-methylimidazo[4,5-c]pyridine'

    def test_5_methylimidazo_4_5_d_1_2_3_triazine(self):
        assert smiles_to_iupac('Cn1cnc2nnncc21') == '5-methylimidazo[4,5-d][1,2,3]triazine'

    def test_1_methylimidazo_4_5_e_pyrazine(self):
        assert smiles_to_iupac('Cn1cnc2nccnc21') == '1-methylimidazo[4,5-e]pyrazine'

    def test_7_methylimidazo_4_5_d_1_2_3_triazine(self):
        assert smiles_to_iupac('Cn1cnc2cnnnc21') == '7-methylimidazo[4,5-d][1,2,3]triazine'

    def test_5_methylimidazo_5_4_e_1_2_4_triazine(self):
        assert smiles_to_iupac('Cn1cnc2nncnc21') == '5-methylimidazo[5,4-e][1,2,4]triazine'

    def test_2_methylpyrazolo_3_4_b_pyridine(self):
        assert smiles_to_iupac('Cn1cc2cccnc2n1') == '2-methylpyrazolo[3,4-b]pyridine'

    def test_6_methylpyrazolo_3_4_d_1_2_3_triazine(self):
        assert smiles_to_iupac('Cn1cc2cnnnc2n1') == '6-methylpyrazolo[3,4-d][1,2,3]triazine'

    def test_2_methylpyrazolo_3_4_d_pyridazine(self):
        assert smiles_to_iupac('Cn1cc2cnncc2n1') == '2-methylpyrazolo[3,4-d]pyridazine'

    def test_2_methylpyrazolo_3_4_e_pyrazine(self):
        assert smiles_to_iupac('Cn1cc2nccnc2n1') == '2-methylpyrazolo[3,4-e]pyrazine'

    def test_2_methylpyrazolo_4_3_c_pyridazine(self):
        assert smiles_to_iupac('Cn1cc2nnccc2n1') == '2-methylpyrazolo[4,3-c]pyridazine'

    def test_2_methylpyrazolo_4_3_d_pyrimidine(self):
        assert smiles_to_iupac('Cn1cc2ncncc2n1') == '2-methylpyrazolo[4,3-d]pyrimidine'

    def test_2_methylpyrazolo_4_3_e_1_2_4_triazine(self):
        assert smiles_to_iupac('Cn1cc2ncnnc2n1') == '2-methylpyrazolo[4,3-e][1,2,4]triazine'

    def test_1_methylpyrazolo_4_5_b_pyridine(self):
        assert smiles_to_iupac('Cn1ncc2ncccc21') == '1-methylpyrazolo[4,5-b]pyridine'

    def test_1_methylpyrazolo_5_4_c_pyridazine(self):
        assert smiles_to_iupac('Cn1ncc2ccnnc21') == '1-methylpyrazolo[5,4-c]pyridazine'

    def test_1_methylpyrazolo_5_4_c_pyridine(self):
        assert smiles_to_iupac('Cn1ncc2ccncc21') == '1-methylpyrazolo[5,4-c]pyridine'

    def test_1_methylpyrazolo_5_4_d_pyrimidine(self):
        assert smiles_to_iupac('Cn1ncc2cncnc21') == '1-methylpyrazolo[5,4-d]pyrimidine'

    def test_2_methyl_1_2_3_triazolo_4_5_b_pyridine(self):
        assert smiles_to_iupac('Cn1nc2cccnc2n1') == '2-methyl[1,2,3]triazolo[4,5-b]pyridine'

    def test_2_methyl_1_2_3_triazolo_4_5_c_pyridazine(self):
        assert smiles_to_iupac('Cn1nc2ccnnc2n1') == '2-methyl[1,2,3]triazolo[4,5-c]pyridazine'

    def test_2_methyl_1_2_3_triazolo_4_5_d_1_2_3_triazine(self):
        assert smiles_to_iupac('Cn1nc2cnnnc2n1') == '2-methyl[1,2,3]triazolo[4,5-d][1,2,3]triazine'

    def test_2_methyl_1_2_3_triazolo_4_5_d_pyrimidine(self):
        assert smiles_to_iupac('Cn1nc2cncnc2n1') == '2-methyl[1,2,3]triazolo[4,5-d]pyrimidine'

    def test_2_methyl_1_2_3_triazolo_4_5_e_1_2_4_triazine(self):
        assert smiles_to_iupac('Cn1nc2ncnnc2n1') == '2-methyl[1,2,3]triazolo[4,5-e][1,2,4]triazine'

    def test_2_methyl_1_2_3_triazolo_4_5_e_pyrazine(self):
        assert smiles_to_iupac('Cn1nc2nccnc2n1') == '2-methyl[1,2,3]triazolo[4,5-e]pyrazine'

    def test_1_methylimidazo_4_5_b_pyridine(self):
        assert smiles_to_iupac('Cn1cnc2ncccc21') == '1-methylimidazo[4,5-b]pyridine'

    def test_1_methylimidazo_4_5_b_quinoline(self):
        assert smiles_to_iupac('Cn1cnc2nc3ccccc3cc21') == '1-methylimidazo[4,5-b]quinoline'

    def test_9_methylpurine(self):
        assert smiles_to_iupac('Cn1cnc2cncnc21') == '9-methylpurine'

    def test_1_methylnaphtho_2_1_d_imidazole(self):
        assert smiles_to_iupac('Cn1cnc2ccc3ccccc3c21') == '1-methylnaphtho[2,1-d]imidazole'
