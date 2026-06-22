# smiles2iupac Handoff — Phase 633 Complete

**Date:** 2026-06-20  
**Project:** `/Users/nobuya/Projects/smiles2iupac`  
**Branch:** `main`  
**Test count:** 7390 passing, 2 skipped  
**Standing instruction:** "keep on" — systematically expand IUPAC 2013 naming coverage

---

## What Was Done This Session

### Phase 633: 4,5,6,7-Tetrahydrobenzo-Fused 5-Membered Rings

Added retained names + locant maps + 34 tests for:
- `4,5,6,7-tetrahydrobenzofuran` (SMILES key: `c1cc2c(o1)CCCC2`)
- `4,5,6,7-tetrahydrobenzothiophene` (key: `c1cc2c(s1)CCCC2`)
- `4,5,6,7-tetrahydro-1H-indole` (key: `c1cc2c([nH]1)CCCC2`)
- `4,5,6,7-tetrahydro-1H-benzimidazole` (key: `c1nc2c([nH]1)CCCC2`)
- `4,5,6,7-tetrahydro-1H-indazole` (key: `c1n[nH]c2c1CCCC2`)

All positions tested including N-1-methyl variants.

### Infrastructure Fixes (broadly applicable)

Three related bugs fixed in `_try_fused_hetero_retained` (line ~2590 of `heterocycle_handler.py`):

1. **`uniquify=False` in `GetSubstructMatches`** — was `uniquify=True` (default). Change at lines ~2733–2735. Now returns all symmetric matches, allowing `min(matches, key=_sub_locants)` to correctly pick the IUPAC-preferred (lowest) locant for C2-symmetric compounds.

2. **Penalty for substituents at None-locant atoms** — `_sub_locants` function (line ~2740) now adds `10000` instead of silently skipping when a substituent lands on a junction/unlabeled atom. Prevents symmetric matches that "lose" a substituent from winning via empty-list comparison.

3. **Fixed pyrene locant map** (line ~2113) — atoms in the same D2h equivalence class now all share the minimum locant:
   - Class {1,3,6,8} → all get locant `1`  
   - Class {2,7} → all get locant `2`  
   - Class {4,5,9,10} → all get locant `4`

### NH Locant Maps Fixed

Added proper locant `1` for NH atoms in:
- `c1cc2c([nH]1)CCCC2`: atom 4 (N-1) now locant `1` (was `None`)
- `c1nc2c([nH]1)CCCC2`: atom 4 (N-1) → `1`, atom 1 (N-3) → `3` (were `None`)
- `c1n[nH]c2c1CCCC2`: atom 2 (N-1) → `1` (was `None`)

### Tests Updated to Use IUPAC-Preferred Locants

Six tests corrected to use the lower (preferred) locant for symmetric compounds:
| File | Old expected | New expected | Reason |
|---|---|---|---|
| `test_phase551.py` line 16 | `6-hydroxyxanthen-9-one` | `3-hydroxyxanthen-9-one` | C2 sym, 3<6 |
| `test_phase552.py` line 13 | `3-methylquinoxaline` | `2-methylquinoxaline` | C2 sym, 2<3 |
| `test_phase552.py` line 22 | `4-methylphthalazine` | `1-methylphthalazine` | C2 sym, 1<4 |
| `test_phase556.py` line 19 | `6-methyl-1,3-benzodioxole` | `5-methyl-1,3-benzodioxole` | C2 sym, 5<6 |
| `test_phase624.py` line 40 | `8-methylxanthen-9-one` | `1-methylxanthen-9-one` | C2 sym, 1<8 |
| `test_phase624.py` line 41 | `7-methylxanthen-9-one` | `2-methylxanthen-9-one` | C2 sym, 2<7 |

---

## Key Technical Context

### Architecture
- **Main naming function:** `smiles_to_iupac(smiles)` in `src/smiles2iupac/__init__.py`
- **Fused heterocycle/retained-name logic:** `src/smiles2iupac/heterocycle_handler.py`
  - `_FUSED_HETERO_RETAINED`: dict `canonical_smiles → retained_name_string`
  - `_FUSED_LOCANT_MAP`: dict `canonical_smiles → {atom_idx: locant_or_None}`
  - `_try_fused_hetero_retained()` at line ~2590: looks up compound, finds substituents, assembles name

### How to Add a New Retained-Name Compound
1. Get canonical SMILES: `from rdkit.Chem import MolFromSmiles, MolToSmiles; MolToSmiles(MolFromSmiles(smi))`
2. Verify name with OPSIN: `echo "name" | java -jar /Users/nobuya/Projects/smiles2iupac/opsin-cli-2.9.0-jar-with-dependencies.jar -osmi`
3. Add to `_FUSED_HETERO_RETAINED` (around line 830–1600, organized by compound type)
4. Compute locant map: add methyl at each substitutable ring atom → compare canonical SMILES with OPSIN output for known named positions
5. Add to `_FUSED_LOCANT_MAP` (after line ~1800, organized similarly)
6. **For symmetric compounds:** all atoms in the same equivalence class MUST share the same (minimum) locant in the map

### Java/OPSIN
```
JAVA=/opt/homebrew/Cellar/openjdk/26.0.1/libexec/openjdk.jdk/Contents/Home/bin/java
JAR=/Users/nobuya/Projects/smiles2iupac/opsin-cli-2.9.0-jar-with-dependencies.jar
echo "compound-name" | $JAVA -jar $JAR -osmi
```

### Python venv
```
.venv/bin/python -m pytest tests/ -q
```

### Naming convention for retained-name substituents
- When parent name contains a locant (e.g., `fluoren-9-one`, `quinolizine`): no hyphen before substituent prefix → `1-methylfluoren-9-one`, `2-methylquinolizine`
- Normal fused heterocycles: hyphen → `2-methylquinoline`

---

## Recent Phase History (Phases 630–633)

| Phase | Content |
|---|---|
| 630 | Benzo-fused 7-membered saturated/partially unsaturated rings (benzazepine, benzoxepine, etc.) |
| 631 | Tricyclic partially saturated: 1,2,3,4-tetrahydrophenanthrene, 9,10-dihydrophenanthrene, etc. + fluoren-9-one fix |
| 632 | 5,6,7,8-tetrahydro N-heterocycles (tetrahydroquinoline, tetrahydroisoquinoline, etc.) + quinolizine + benz[a]anthracene, benzo[a/b/c]fluorene, dibenz[a,h]anthracene |
| 633 | 4,5,6,7-tetrahydrobenzo-fused 5-membered rings (this session) |

---

## Suggested Next Phases

Good candidates to continue with:
1. **4,5,6,7-tetrahydrobenzo-fused 5-membered rings with 2 N** — e.g., `4,5,6,7-tetrahydro-1H-benzotriazole`, `4,5,6,7-tetrahydro-1H-benzo[d][1,2,3]triazole`
2. **2,3-dihydro-1H-pyrrolo/imidazo benzo-fused** — partially saturated with one sp3 C in 5-membered ring
3. **Octahydro/perhydro fully saturated bridged systems** — decalin variants
4. **Additional partially saturated naphthyridines** — 1,2,3,4-tetrahydronaphthyridines
5. **Partially saturated acridine, phenanthridine** — similar to tetrahydroquinoline pattern

---

## Suggested Skills

- `/memory` or check memory files at `/Users/nobuya/.claude/projects/-Users-nobuya-Projects-smiles2iupac/memory/` for user preferences and project context
- When running tests: `.venv/bin/python -m pytest tests/ -q`
- When verifying IUPAC names: use OPSIN (see path above)
- When checking for regressions after adding new entries: always run full test suite before writing new test file
