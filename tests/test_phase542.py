"""Phase 542: Carbaldehyde hydrazone/semicarbazone/thiosemicarbazone on rings.

Three root causes fixed:
1. _name_semicarbazone used find_principal_chain which gives 1-carbon chain for
   ring-adjacent anchors → "methanal semicarbazone". Now detects aromatic ring
   neighbor and returns "{ring}carbaldehyde {suffix}" directly.
2. name_heterocycle's _EXOCYCLIC_SUFFIX lacked "aldhydrazone" (unsubstituted
   hydrazone bypasses PGRP_DISPATCH and hits name_heterocycle).
3. constants.py aldhydrazone had benzene_name=None; set to "benzaldehyde hydrazone"
   so _name_cyclic generates the correct benzene name.
"""
import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Aldhydrazone on heteroaromatic rings
    ("NN=Cc1ccccn1",    "pyridine-2-carbaldehyde hydrazone"),
    ("NN=Cc1cccnc1",    "pyridine-3-carbaldehyde hydrazone"),
    ("NN=Cc1cccs1",     "thiophene-2-carbaldehyde hydrazone"),
    ("NN=Cc1ccc[nH]1",  "1H-pyrrole-2-carbaldehyde hydrazone"),
    # Aldsemicarbazone on heteroaromatic rings
    ("NC(=O)NN=Cc1ccccn1", "pyridine-2-carbaldehyde semicarbazone"),
    ("NC(=O)NN=Cc1cccnc1", "pyridine-3-carbaldehyde semicarbazone"),
    ("NC(=O)NN=Cc1cccs1",  "thiophene-2-carbaldehyde semicarbazone"),
    # Aldthiosemicarbazone on heteroaromatic rings
    ("NC(=S)NN=Cc1ccccn1", "pyridine-2-carbaldehyde thiosemicarbazone"),
    # Benzene (now fixed too)
    ("NN=Cc1ccccc1",        "benzaldehyde hydrazone"),
    ("NC(=O)NN=Cc1ccccc1",  "benzaldehyde semicarbazone"),
    ("NC(=S)NN=Cc1ccccc1",  "benzaldehyde thiosemicarbazone"),
    # Regression: aliphatic unchanged
    ("NN=CC",           "ethanal hydrazone"),
    ("NC(=O)NN=CC",     "ethanal semicarbazone"),
    ("NC(=S)NN=CC",     "ethanal thiosemicarbazone"),
])
def test_phase542_hydrazone_semicarbazone_ring(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
