# OPSIN Verification Skip List

Entries in `_FUSED_HETERO_RETAINED` that cannot have `_FUSED_LOCANT_MAP` entries built via OPSIN.
Recorded during Phase 616–620 locant-map construction (2026-06-19).

---

## Category 1: OPSIN returns None (unsupported name)

These compound names are not implemented by OPSIN (v2.9.0).
Methyl-derivative locants cannot be verified, so locant maps cannot be built.

| Key SMILES | Registered name | Note |
|---|---|---|
| `c1cn2ncncc2n1` | `imidazo[1,2-f][1,2,4]triazine` | f-descriptor fused ring; OPSIN unsupported |
| `c1cc2nnnn2nn1` | `tetrazolo[1,5-f][1,2,3]triazine` | f-descriptor fused ring; OPSIN unsupported |
| `c1ncc2nnnn2n1` | `tetrazolo[1,5-f][1,2,4]triazine` | f-descriptor fused ring; OPSIN unsupported |

---

## Category 2: Wrong SMILES key (SMILES ≠ name)

The SMILES registered as a key does not represent the named compound
(discovered via InChI mismatch between key SMILES and OPSIN-generated SMILES).
The correct SMILES exists elsewhere in `_FUSED_HETERO_RETAINED` or needs a separate fix.

| Key SMILES (wrong) | Registered name | Correct canonical SMILES |
|---|---|---|
| `c1cc2n[nH]cc2cn1` | `1H-pyrazolo[3,4-c]pyridine` | `c1cc2c[nH]nc2cn1` (different isomer) |
| `c1c[nH]c2ncnc-2n1` | `1H-pyrazolo[3,4-d]pyrimidine` | `c1ncc2c[nH]nc2n1` (different isomer) |
| `c1ccc2c(c1)CCO2` | `1,3-dihydro-2-benzofuran` | `c1ccc2c(c1)COC2` (O position differs) |
| `c1ccc2c(c1)CCS2` | `1,3-dihydro-2-benzothiophene` | `c1ccc2c(c1)CSC2` (S position differs) |
| `C1=Nc2cccc3cccc1c23` | `perimidine` | `C1=Nc2cccc3cccc(c23)N1` (**fixed** 2026-06-19) |

---

## Category 3: Fully substituted retained names (no free C for methyl)

These entries are retained names for specific substituted compounds
(e.g. purines with =O or amino groups). They have no unsubstituted C atoms
available for further methyl substitution, so locant maps are unnecessary.

| Key SMILES | Registered name |
|---|---|
| `Cn1c(=O)c2c(ncn2C)n(C)c1=O` | `1,3,7-trimethyl-3,7-dihydro-1H-purine-2,6-dione` (caffeine) |
| `Cn1c(=O)c2[nH]cnc2n(C)c1=O` | `1,3-dimethyl-3,7-dihydro-1H-purine-2,6-dione` (theophylline) |
| `Cn1cnc2c1c(=O)[nH]c(=O)n2C` | `3,7-dimethyl-3,7-dihydro-1H-purine-2,6-dione` (theobromine) |
| `CN1C(=O)NC2=NC=NC21` | `3-methyl-3,7-dihydro-1H-purine-2,6-dione` |
| `CN1NC2N=CN=C2C1=O` | `1-methyl-3,7-dihydro-1H-purine-2,6-dione` |
| `O=c1[nH]cnc2nc[nH]c12` | `1,7-dihydro-6H-purin-6-one` (hypoxanthine) |
| `O=c1[nH]c(=O)c2[nH]c(=O)[nH]c2[nH]1` | `7,9-dihydro-1H-purine-2,6,8(3H)-trione` (uric acid) |
| `Nc1nc2[nH]cnc2c(=O)[nH]1` | `2-amino-3,7-dihydro-1H-purin-6-one` (guanine) |

---

## Summary

| Category | Count | Action |
|---|---|---|
| OPSIN unsupported (f-descriptor) | 3 | Skip permanently |
| Wrong SMILES key | 4 (+1 fixed) | Fix keys in `_FUSED_HETERO_RETAINED` |
| Fully substituted retained names | 8 | No locant map needed |
| **Total missing** | **15** | — |

`_FUSED_HETERO_RETAINED` total: 549 entries  
`_FUSED_LOCANT_MAP` coverage: 534 / 549 (97.3%)
