"""Phase 821: tetrazolo-fused heteroaromatic α-ol/thiol → tautomers.

- tetrazolo[1,5-a]pyridine C5 → 5(1H)-one/thione
- tetrazolo[1,5-a]pyrimidine C5 → 5(1H)-one/thione; C7 → 7(1H)-one/thione
- tetrazolo[1,5-b]pyridazine C6 → 6(4H)-one/thione
- tetrazolo[1,5-b][1,2,4]triazine C6 → 6(4H)-one/thione; C7 → 7(1H)-one/thione
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # tetrazolo[1,5-a]pyridine C5-OH/SH
    ("Oc1cccc2nnnn12",   "tetrazolo[1,5-a]pyridin-5(1H)-one"),
    ("Sc1cccc2nnnn12",   "tetrazolo[1,5-a]pyridin-5(1H)-thione"),
    # tetrazolo[1,5-a]pyrimidine C5-OH/SH
    ("Oc1ccn2nnnc2n1",   "tetrazolo[1,5-a]pyrimidin-5(1H)-one"),
    ("Sc1ccn2nnnc2n1",   "tetrazolo[1,5-a]pyrimidin-5(1H)-thione"),
    # tetrazolo[1,5-a]pyrimidine C7-OH/SH
    ("Oc1ccnc2nnnn12",   "tetrazolo[1,5-a]pyrimidin-7(1H)-one"),
    ("Sc1ccnc2nnnn12",   "tetrazolo[1,5-a]pyrimidin-7(1H)-thione"),
    # tetrazolo[1,5-b]pyridazine C6-OH/SH
    ("Oc1ccc2nnnn2n1",   "tetrazolo[1,5-b]pyridazin-6(4H)-one"),
    ("Sc1ccc2nnnn2n1",   "tetrazolo[1,5-b]pyridazin-6(4H)-thione"),
    # tetrazolo[1,5-b][1,2,4]triazine C6-OH/SH
    ("Oc1cnc2nnnn2n1",   "tetrazolo[1,5-b][1,2,4]triazin-6(4H)-one"),
    ("Sc1cnc2nnnn2n1",   "tetrazolo[1,5-b][1,2,4]triazin-6(4H)-thione"),
    # tetrazolo[1,5-b][1,2,4]triazine C7-OH/SH
    ("Oc1cnn2nnnc2n1",   "tetrazolo[1,5-b][1,2,4]triazin-7(1H)-one"),
    ("Sc1cnn2nnnc2n1",   "tetrazolo[1,5-b][1,2,4]triazin-7(1H)-thione"),
    # Regressions: parent rings unchanged
    ("c1ccn2nnnc2c1",    "tetrazolo[1,5-a]pyridine"),
    ("c1cnc2nnnn2c1",    "tetrazolo[1,5-a]pyrimidine"),
    ("c1cnn2nnnc2c1",    "tetrazolo[1,5-b]pyridazine"),
    ("c1cnn2nnnc2n1",    "tetrazolo[1,5-b][1,2,4]triazine"),
])
def test_phase821_tetrazolo_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
