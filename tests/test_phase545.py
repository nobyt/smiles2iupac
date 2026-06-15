"""Phase 545: N-substituted aldehyde hydrazone on rings.

_name_substituted_hydrazone called find_principal_chain which gives a 1-carbon
chain for ring-adjacent anchors, producing 'methanal N-methylhydrazone'.
Added the same ring-detection block as _name_semicarbazone (phase 542): find
an aromatic ring neighbor, use _aryl_sulfonyl_prefix, build the carbaldehyde
parent name, then attach the N-alkyl prefix.
"""
import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # N-substituted aldhydrazone on heteroaromatic rings
    ("CNN=Cc1ccccn1",    "pyridine-2-carbaldehyde N-methylhydrazone"),
    ("CCNN=Cc1ccccn1",   "pyridine-2-carbaldehyde N-ethylhydrazone"),
    ("CNN=Cc1cccnc1",    "pyridine-3-carbaldehyde N-methylhydrazone"),
    ("CNN=Cc1cccs1",     "thiophene-2-carbaldehyde N-methylhydrazone"),
    ("CNN=Cc1ccc[nH]1",  "1H-pyrrole-2-carbaldehyde N-methylhydrazone"),
    # Benzene (now fixed too)
    ("CNN=Cc1ccccc1",    "benzaldehyde N-methylhydrazone"),
    # Regression: aliphatic and unsubstituted unchanged
    ("CNN=CC",           "ethanal N-methylhydrazone"),
    ("NN=Cc1ccccn1",     "pyridine-2-carbaldehyde hydrazone"),
    ("NN=Cc1ccccc1",     "benzaldehyde hydrazone"),
])
def test_phase545_n_sub_hydrazone_ring(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
