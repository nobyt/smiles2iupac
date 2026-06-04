"""Phase 389: Guanidine amine-N labelling is SMILES-order independent (IUPAC 2013 P-66.4).

Previously the two amine N atoms in guanidine were labelled N / N' in atom-index
order (which depends on SMILES parsing direction), so CNC(=N)N gave
'N-methylguanidine' but NC(=N)NC gave 'N''-methylguanidine' for the same molecule.

Fix: sort amine nitrogens by descending substituent count first, then alphabetically
by substituent name, so the N with more / earlier-alphabetical substituents always
gets the unprimed 'N' locant.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # mono N-substituted (same compound, two SMILES directions)
    ("CNC(=N)N",            "N-methylguanidine"),
    ("NC(=N)NC",            "N-methylguanidine"),
    # symmetric N,N'-disubstituted
    ("CNC(=N)NC",           "N,N'-dimethylguanidine"),
    ("CCNC(=N)NCC",         "N,N'-diethylguanidine"),
    # asymmetric: alphabetical order (ethyl before methyl)
    ("CNC(=N)NCC",          "N-ethyl-N'-methylguanidine"),
    ("CCNC(=N)NC",          "N-ethyl-N'-methylguanidine"),
    # unsubstituted (regression)
    ("NC(=N)N",             "guanidine"),
    # N-methyl, N'-methyl with ethyl on imine N
    ("CNC(=NCC)NC",         "N''-ethyl-N,N'-dimethylguanidine"),
    # N,N-disubstituted on one amine N
    ("CN(C)C(=N)N",         "N,N-dimethylguanidine"),
    ("NC(=N)N(C)C",         "N,N-dimethylguanidine"),
])
def test_phase389_guanidine_n_labelling(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
