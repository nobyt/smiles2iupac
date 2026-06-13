"""Phase 467: thiazolo/oxazolo/imidazo[x,y-c]pyridine,
1H-[1,2,3]triazolo[4,5-b]pyridine, pyrazolo[3,4-c] and [4,5-c]
pyridine/pyridazine isomers (IUPAC 2013 P-31.1.3 fusion nomenclature).
Also fixes Phase 465 pyrazolo[3,4-c]pyridine → pyrazolo[4,5-c]pyridine.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # c-bond pyridine fusions with S/O/N heteroatom in ring5
    ("c1cc2ncsc2cn1",     "thiazolo[5,4-c]pyridine"),
    ("c1cc2ncoc2cn1",     "oxazolo[5,4-c]pyridine"),
    ("c1cc2[nH]cnc2cn1",  "1H-imidazo[4,5-c]pyridine"),
    # b-bond pyridine fusion: 1H-[1,2,3]triazolo
    ("c1cnc2nn[nH]c2c1",  "1H-[1,2,3]triazolo[4,5-b]pyridine"),
    # pyrazolo c-bond isomers: [4,5-c] (corrected from Phase 465) and [3,4-c] (new)
    ("c1cc2[nH]ncc2cn1",  "1H-pyrazolo[4,5-c]pyridine"),
    ("c1cc2c[nH]nc2cn1",  "1H-pyrazolo[3,4-c]pyridine"),
    # pyrazolo c-bond pyridazine isomers
    ("c1cc2[nH]ncc2nn1",  "1H-pyrazolo[4,5-c]pyridazine"),
    ("c1cc2c[nH]nc2nn1",  "1H-pyrazolo[3,4-c]pyridazine"),
])
def test_phase467_heteroaromatic_c_fusions(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
