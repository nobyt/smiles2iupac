"""Phase 544: Two substituent-naming fixes for locant-bearing heteroaryl groups.

Fix A — _make_oxy_name: "pyridin-2-yl" → stripped "pyridin-2-" + "oxy"
= "pyridin-2-oxy" (wrong). Base ending with "-" signals a heteroaryl-yl
name; keep the full "-yl" suffix before appending "oxy": "pyridin-2-yloxy".

Fix B — _name_sulfonamide N-substituent: "pyridin-2-yl" lacks parens in
"N-pyridin-2-yl..." output. Extended the sub_str rule to also wrap when the
name contains a digit, matching the same rule as S-thioate (phase 543).
Benzene/methyl (no digit) are unaffected.
"""
import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Fix A: heteroaryloxy substituent uses "-yloxy" not "-oxy"
    ("OC(=O)COc1ccccn1",    "2-(pyridin-2-yloxy)acetic acid"),
    ("OC(=O)CCOc1ccccn1",   "3-(pyridin-2-yloxy)propanoic acid"),
    ("OC(=O)COc1cccnc1",    "2-(pyridin-3-yloxy)acetic acid"),
    ("OC(=O)COc1cccs1",     "2-(thiophen-2-yloxy)acetic acid"),
    # Fix B: N-sulfonamide substituent parentheses
    ("CS(=O)(=O)Nc1ccccn1", "N-(pyridin-2-yl)methanesulfonamide"),
    ("CS(=O)(=O)Nc1cccnc1", "N-(pyridin-3-yl)methanesulfonamide"),
    ("CS(=O)(=O)Nc1cccs1",  "N-(thiophen-2-yl)methanesulfonamide"),
    # Regression A: simple ether names unchanged
    ("OC(=O)COc1ccccc1",    "2-phenoxyacetic acid"),
    ("OC(=O)COCC",          "2-ethoxyacetic acid"),
    # Regression B: simple N-subs unchanged
    ("CS(=O)(=O)Nc1ccccc1", "N-phenylmethanesulfonamide"),
    ("CS(=O)(=O)NC",        "N-methylmethanesulfonamide"),
])
def test_phase544_heteroaryl_substituent_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
