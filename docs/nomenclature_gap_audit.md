# Nomenclature Gap Audit (Track B)

Generated 2026-08-18 per the "implement missing IUPAC features" plan. Unlike
`docs/locant_map_audit.md` (Track A: known-wrong answers from an existing
code path), this doc tracks **rule classes with no code path at all**.

## Method

1. **B1 — mine `tests/pubchem_cache.json`**: ran every `validated: False`
   cache entry (4,011 of 6,978) through `smiles_to_iupac`, diffed against
   PubChem's `IUPACName`, and classified mismatches with the existing
   `iupac_verdict.classify_mismatch()`.
2. **B2 — targeted probes**: directly exercised specific rule classes
   suspected to be unimplemented based on `README.md`'s coverage table and
   the small size of `stereochemistry.py` (64 lines).

## B1 results: real-world corpus shows almost no hard gaps

| Bucket | Count | Meaning |
|---|---|---|
| Exact match | 1 | (most cache entries are mismatches by construction — this cache tracks known diffs) |
| Crash (exception) | 1 | `CC=N(CC)CC` — `ValueError: Invalid SMILES`. Root cause: the cached SMILES itself is chemically invalid (neutral trivalent N can't carry a double bond + 2 single bonds without a formal + charge) — a **bad cache entry**, not an engine gap. |
| `pubchem_wrong` | 392 | Already-classified: PubChem violates IUPAC 2013 formatting (not our gap) |
| `pubchem_retained` / `our_retained` / `tautomer` | 47 | Already-classified retained-name/tautomer preference differences |
| `needs_review` (unclassified) | 1,947 | See below |

The 1,947 "needs_review" rows were spot-checked (30-row sample) and are
overwhelmingly **PubChem choosing a trivial/retained name where we produce
the systematic PIN** (`alanine` vs `2-aminopropanoic acid`, `malonic acid`
vs `propanedioic acid`, `cumene` vs `(propan-2-yl)benzene`, `citric acid`
vs `2-hydroxypropane-1,2,3-tricarboxylic acid`) or a **locant-omission
convention difference** (`nitroethane` vs PubChem's `1-nitroethane` on an
unambiguous monosubstituted chain). None of the sampled rows indicate a
missing rule class — they're the same "which name is PIN" judgment calls
this project's PubChem-verdict pipeline already triages routinely.

**Conclusion: the engine's functional-group/ring coverage is already close
to complete for real-world compounds.** B1 does not surface a hidden crash
backlog. The genuine gaps are structural (whole descriptor classes), found
instead by B2's targeted probing.

## B2 results: three confirmed structural gaps

### 1. Isotopic labeling — not implemented at all
```python
smiles_to_iupac("[2H]C([2H])([2H])C(=O)O")   # -> "acetic acid"   (deuterium silently dropped)
smiles_to_iupac("C[13CH2]C(=O)O")            # -> "propanoic acid" (13C silently dropped)
```
No isotope-related code exists anywhere in `src/` (confirmed by grep). IUPAC
2013 P-82 isotope nomenclature (e.g. `(2,2,2-²H₃)acetic acid`) has zero
support. Any isotope-labeled input silently loses the label with no error —
same silent-wrongness shape as the Track A locant bugs, but structural
rather than a table lookup gap.

### 2. Axial chirality (allenes, atropisomers) — not implemented
```python
smiles_to_iupac("CC(F)=C=C(Cl)C")              # -> "4-chloro-2-fluoropenta-2,3-diene" (no descriptor)
smiles_to_iupac(r"C(/C(=C=C(\C)F)Cl)")         # -> "2-chloro-4-fluoropenta-2,3-diene" (no descriptor, DIFFERENT locants for the same input intent)
```
`stereochemistry.py` (64 lines) only maps CIP tetrahedral/double-bond codes
to `(R)`/`(E)`-style descriptors; it has no path for axial (`Ra`/`Sa`) or
planar chirality. Worse than a simple omission: the two input
representations above encode the same substitution pattern with opposite
explicit stereo bonds, yet produce differently-numbered names with neither
descriptor — suggesting the numbering routine isn't even stably ignoring
the stereo bonds, it's arbitrarily picking different equivalent atoms.

**Investigated 2026-08-19: this is a hard blocker, not just missing code.**
This project delegates all stereo perception to RDKit (per `README.md`).
The installed RDKit (2026.03.2) does not preserve or perceive allene axial
chirality through any of its documented SMILES-parsing paths: `[C@]`/
`[C@@]` on an allene's central sp carbon (e.g.
`"ClC(Br)=[C@]=C(F)I"`) parses without error but the atom's
`GetChiralTag()` comes back `CHI_UNSPECIFIED` — confirmed with both legacy
and non-legacy stereo perception (`Chem.SetUseLegacyStereoPerception(False)`
+ `Chem.FindPotentialStereo`, which only reports the two cumulated bonds as
generic `Bond_Double` elements, never a dedicated axial/allene stereo
type). Both `@`/`@@` variants of the same input canonicalize to the
identical SMILES string — RDKit is silently discarding the distinction at
parse time, not merely failing to label it. Implementing axial chirality
naming on top of this would require either (a) a custom CIP/stereo-
perception implementation bypassing RDKit for this one case, or (b) an
RDKit version bump to one with proper allene/atropisomer support — both
materially bigger scope decisions than "add a naming feature," and (b) in
particular isn't something to decide unilaterally given this project's
large, mature test suite's dependency on RDKit's current behavior.
**Recommendation: do not attempt until either is explicitly decided**, not
merely deprioritized.

### 3. Complex (non-simple) polycyclic von Baeyer systems — confirmed wrong, not just unsupported
```python
smiles_to_iupac("C1CC2CCC1CC2")        # -> "bicyclo[2.2.2]octane"  (correct)
smiles_to_iupac("C1CC2CC1C1CCCCC21")   # -> "cyclohexane"           (WRONG — silently drops 5 of 11 ring atoms)
```
The second SMILES is a genuine, valid, connected tricyclic hydrocarbon (11
atoms, 3 SSSR rings of size 5/5/6 — confirmed via RDKit `GetRingInfo`).
Simple bridged bicyclics (`polycyclic_handler.py`) and retained names
(adamantane) work; this shows the fallback path for *irregular* tricyclic+
fused/bridged combinations doesn't fail loudly — it falls through to a
handler that names only a fragment as if it were the whole molecule. This
is the most severe of the three B2 findings (wrong answer with no error,
versus 1/2's silent-drop-but-otherwise-sane names) and the best B3
starting point given it's a correctness bug in existing code, not
new-nomenclature-from-scratch work like isotopes/axial chirality.

**Status: partially addressed in Phase 932.** Root cause was
`find_principal_ring()`'s multi-ring fallback (`ring_handler.py`, "複数環:
最大環を選択") — for any molecule not matching one of the specific shapes
handled earlier (spiro/bicyclo/cage-retained/naphthalene/anthracene/
phenanthrene/PAH/biphenyl), it picks whichever SSSR ring is largest and
names *only* that ring; any other ring sharing atoms with it (i.e. actually
fused/bridged, not just chain-connected) has its atoms silently discarded.
Phase 932 added a guard: when another ring shares atoms with the chosen
principal ring, raise `ValueError` instead of emitting the wrong name.
**This does NOT implement general fused/bridged 3+-ring von Baeyer
nomenclature** (IUPAC P-23.2.3-5 main-ring/main-bridge selection is still
unimplemented) — it only converts the silent-wrong-answer failure mode
into a loud, honest one. Actually naming these systems correctly (e.g.
`tricyclo[5.2.1.0²,⁶]decane`-style names) remains open B3 work. The guard
was initially too broad (fired for ANY other ring anywhere in the
molecule, including chain-connected separate rings like diphenylmethane's
two phenyls, which `collect_ring_substituents` already handles correctly)
— narrowed to only fire when a ring actually shares atoms with the chosen
one; caught by the full test suite (8 failures on the first attempt, 0 on
the corrected version) — see `tests/test_phase932.py`.

### 4. Multi-component (salt/hydrate) naming produced invalid or non-standard names

```python
smiles_to_iupac("O.O.CC(=O)[O-].[Ca+2].[O-]C(C)=O")   # was -> "calcium diacetate diwater" (non-standard; real term is "dihydrate")
smiles_to_iupac("CC(=O)O.CC(=O)O")                     # -> "diacetic acid" (confirmed OPSIN-unparseable — not a valid name at all)
```
`_name_multicomponent()` (`__init__.py`, Phase 147) was originally scoped
and tested only for charged salt components (cation/anion names getting a
multiplier prefix, e.g. `"calcium diacetate"`, `"dipotassium sulfate"` —
all in `tests/test_phase147.py`). Its generic `_compact()` multiplier
logic was also being applied unconditionally to *neutral* components,
which nobody had verified: water got `"di"`+`"water"` instead of the
IUPAC/CAS-standard `"dihydrate"` term, and two disconnected identical
neutral organic molecules (no salt/hydrate relationship at all) produced
`"diacetic acid"` — an invalid name (OPSIN can't parse it).

**Status: the water/hydrate case fixed in Phase 935** (the common,
well-established real-world case — pharma hydrates like "amoxicillin
trihydrate" are everyday terminology). Water fragments are now pulled out
of the generic neutral-component path and appended as a single
`"{multiplier}hydrate"` suffix (`"mono"` for count=1, since
`"monohydrate"`, not bare `"hydrate"`, is the standard term), verified via
OPSIN round-trip for 1/2/3-water cases. **The non-water case (two
identical disconnected neutral organic molecules) is deliberately NOT
fixed** — it's a synthetic edge case (0 hits in `tests/pubchem_cache.json`'s
6,978 entries) and the correct general IUPAC multi-component convention
for arbitrary neutral organic co-crystals isn't as unambiguous as the
water case; would need its own investigation if a real use case surfaces.

## Recommended B3 priority

1. **General fused/bridged 3+-ring von Baeyer naming** (finding 3,
   continued) — Phase 932 stopped the silent wrong-answer bleeding but did
   not implement the actual nomenclature; affected molecules now raise
   `ValueError` rather than get a name at all. **Checked 2026-08-20: 0 of
   6,978 `tests/pubchem_cache.json` entries trigger this path** — real-world
   frequency in this project's reference corpus is effectively zero.
   Implementing IUPAC P-23.2.3-5's main-ring/main-bridge selection
   algorithm for arbitrary N-ring fused/bridged systems is a substantial,
   high-regression-risk feature for near-zero measured payoff — **not
   recommended as a priority unless a concrete real-world need surfaces**;
   if picked up anyway, scope it as its own session with dedicated
   OPSIN/PubChem verification against a range of topologies, not a quick
   add-on.
2. **Isotopic labeling** (finding 1) — **implemented in Phase 933** for the
   dominant real-world case: isotope labels on principal-chain atoms via
   the acyclic naming path (`_name_acyclic`/`find_principal_chain`).
   `AtomInfo.isotope` (new field) + `format_isotope_descriptor()`
   (`name_assembler.py`) build the P-82 parenthetical descriptor, verified
   against OPSIN including a non-obvious prefix-ordering rule (isotope
   descriptor goes AFTER ordinary substituent prefixes, BEFORE the parent
   name — "2-chloro(2,2-2H2)acetic acid", not the other order). **Still
   open**: isotope labels on ring atoms, substituent branches, and the
   amide/amine/amidine/imine early-return special paths in `_name_acyclic`
   are out of scope — they continue to silently drop the label (no
   regression, just not yet extended to those paths).
3. **Multi-component hydrate naming** (finding 4) — **implemented in Phase
   935** for water; the non-water disconnected-duplicate-neutral-component
   case remains open but is low priority (synthetic edge case, 0 real-world
   hits).
4. **Axial/planar chirality** (finding 2) — highest implementation cost
   (needs new stereo-perception logic beyond what `stereochemistry.py`
   currently does, likely including RDKit's newer stereo API for
   atropisomers), and blocked pending an RDKit-version/custom-CIP decision
   (see investigation note above). Do not attempt without that decision.

Explicitly not pursued in this pass (per plan Track B4): the 3 OPSIN-
unsupported f-locant fused-ring names in `docs/opsin_skip_list.md`, and
anything that would require new work in the Rust port.

## Relative stereodescriptors (`rel-`, `R*`/`S*`) — investigated, not applicable

Briefly considered as a candidate gap (partial/relative-only stereo
configuration, as opposed to fully absolute or fully unspecified). Not
pursued: `rel-`/`R*`/`S*` nomenclature describes a molecule whose *relative*
configuration is known but *absolute* configuration is not (e.g. from a
racemic synthesis where only the diastereomer relationship is
established) — this is metadata that plain SMILES `@`/`@@` notation cannot
express at all (SMILES stereo bonds always encode a specific absolute
configuration once written). Supporting this would require a different
*input* format or an explicit "this is relative-only" flag alongside the
SMILES, not a gap in the naming logic itself — out of scope for a
SMILES-in, name-out pipeline. (Partial stereo — some centers defined,
others left completely unspecified — already works correctly today, e.g.
`smiles_to_iupac("OC(=O)[C@@H](O)C(O)C(=O)O")` → `"(2S)-2,3-dihydroxybutanedioic acid"`,
correctly omitting a descriptor for the undefined center.)
