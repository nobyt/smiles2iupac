"""Phase 873: ring (pseudo-asymmetric) stereochemistry via CIP r/s.

Ring stereocenters whose CIP descriptor is the lowercase r/s "reference"
form -- e.g. 1,4-disubstituted cyclohexanes, which are not chiral but whose
relative configuration (cis/trans) still distinguishes two compounds -- were
dropped entirely: cis- and trans-1,4-dimethylcyclohexane both gave the bare
"1,4-dimethylcyclohexane".

Root cause: the graph builder used the legacy Chem.AssignStereochemistry, which
does not assign _CIPCode to pseudo-asymmetric centers. Switching to
rdCIPLabeler.AssignCIPLabels adds the lowercase r/s labels while leaving every
existing uppercase R/S label identical. rdCIPLabeler moves bond E/Z into the
bond _CIPCode (converting the GetStereo enum to CIS/TRANS), so bond stereo is
now read from _CIPCode with a fallback to the legacy enum.

The codebase already uses CIP descriptors exclusively (no cis/trans words), and
uppercase-R/S ring stereo like (1S,2S)-1,2-dimethylcyclohexane already worked;
this extends the same convention to the lowercase r/s centers.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1,4-disubstituted cyclohexane: pseudo-asymmetric r/s (two distinct isomers)
    ("C[C@H]1CC[C@@H](C)CC1", "(1s,4s)-1,4-dimethylcyclohexane"),
    ("C[C@H]1CC[C@H](C)CC1",  "(1r,4r)-1,4-dimethylcyclohexane"),
    ("O[C@H]1CC[C@@H](O)CC1", "(1s,4s)-cyclohexane-1,4-diol"),
    ("O[C@H]1CC[C@H](O)CC1",  "(1r,4r)-cyclohexane-1,4-diol"),
    # regression: true chiral ring centers (uppercase R/S) unchanged
    ("C[C@H]1CCCC[C@@H]1C", "(1S,2S)-1,2-dimethylcyclohexane"),
    # regression: acyclic R/S unchanged
    ("C[C@H](O)CC",        "(2S)-butan-2-ol"),
    ("C[C@H](N)C(=O)O",    "(2S)-2-aminopropanoic acid"),
    ("F[C@H](Cl)Br",       "(R)-bromochlorofluoromethane"),
    ("C[C@H](Cl)[C@@H](Cl)C", "(2S,3S)-2,3-dichlorobutane"),
    # regression: E/Z double-bond stereo unchanged (read via CIP _CIPCode now)
    ("C/C=C/C",            "(2E)-but-2-ene"),
    (r"C/C=C\C",           "(2Z)-but-2-ene"),
    ("C/C=C/C=C/C",        "(2E,4E)-hexa-2,4-diene"),
    ("ClC=C",              "chloroethene"),
])
def test_phase873_ring_pseudoasymmetric_stereo(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phase873_cis_trans_distinguished():
    # The two 1,4-dimethylcyclohexane diastereomers must get distinct names.
    trans = smiles_to_iupac("C[C@H]1CC[C@@H](C)CC1")
    cis = smiles_to_iupac("C[C@H]1CC[C@H](C)CC1")
    assert trans != cis
