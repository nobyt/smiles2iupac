# smiles2iupac

A library and CLI tool that converts SMILES strings to **IUPAC 2013 preferred names**.

SMILES parsing and molecular graph construction are delegated to [RDKit](https://www.rdkit.org/);
the naming logic is a from-scratch implementation verified by 4,086+ tests.

[日本語版 README](README.ja.md)

---

## Coverage

| Category | Examples |
|---|---|
| Alkanes, alkenes, alkynes | butane, but-2-ene, but-2-yne |
| Halides, alcohols | 2-chloropropane, propan-1-ol |
| Carboxylic acids, esters, amides | acetic acid, ethyl acetate, acetamide |
| Cyclic compounds | cyclohexane, benzene, naphthalene |
| Stereochemistry | (R)-alanine, (E)-but-2-ene |
| Heterocycles | pyridine, furan, imidazole |
| Nitrogen functional groups | amine, nitrile, amidine, hydrazide |
| Sulfur functional groups | thiol, sulfonic acid, sulfoxide |
| Phosphorus and silicon compounds | phosphate ester, trimethylsilanol |
| Other | ethers, peroxides, isocyanates, ylides, multiplicative, 'a' nomenclature, 3+-ring von Baeyer, isotopes, etc. |

---

## Limitations & Unsupported IUPAC 2013 Rules

While `smiles2iupac` achieves near-100% coverage of systematic Preferred IUPAC Names (PIN) for common organic chemistry, the following specialized categories from IUPAC 2013 Blue Book are intentionally out of scope or blocked by upstream libraries:

1. **Axial & Planar Chirality (IUPAC P-92.4, P-92.5)**
   - Substituted allenes (`(Ra)/(Sa)` or `(P)/(M)`), atropisomeric biaryls, helicenes, cyclophane planar chirality.
   - *Upstream constraint*: RDKit's SMILES parser discards allene chiral markers (e.g. `[C@]` in `C/C=[C@]=C/C`) as `CHI_UNSPECIFIED`. See [`docs/RDKIT_ALLENE_STEREO_LIMITATION.md`](docs/RDKIT_ALLENE_STEREO_LIMITATION.md) for details.
2. **Phane Nomenclature (IUPAC P-26)**
   - Macrocyclic phanes such as `[2.2]paracyclophane`.
3. **Fullerenes (IUPAC P-27)**
   - Closed cage carbon allotropes like `[60]fullerene`.
4. **Transition Metal Complexes & Metallocenes (IUPAC P-68)**
   - Ferrocene, coordination polymers, and transition metal organometallics. Main group organometallics (Si, Ge, Sn, Pb, B, P, As, Sb, Bi, Hg) are fully supported.
5. **Systematic Parent Derivation for Giant Natural Products (IUPAC P-28)**
   - Skeletal derivation of complex steroids (`cholestane`) or alkaloids (retained amino acid and simple natural product PINs are supported via lookup tables).
6. **Free Radicals & Carbenium Ions (IUPAC P-61)**
   - Isolated reactive intermediate radicals (`methyl radical`) or cations (`methylium`). Stable onium salts, carboxylates, and alkoxides are fully supported.

---

## Installation

### Requirements

- Python 3.11+
- [uv](https://docs.astral.sh/uv/) or pip
- RDKit 2023.9+

### uv (recommended)

```bash
git clone https://github.com/yourname/smiles2iupac.git
cd smiles2iupac

# Create virtual environment and install dependencies
uv sync

# Optional: install CLI globally
uv tool install .
```

After `uv sync` the CLI is available at `.venv/bin/smiles2iupac`.

### pip

```bash
git clone https://github.com/yourname/smiles2iupac.git
cd smiles2iupac

pip install .
```

> **Note on RDKit**: the `rdkit` PyPI package works on most platforms, but if you
> encounter issues the conda distribution is more stable:
>
> ```bash
> conda install -c conda-forge rdkit
> ```

---

## CLI usage

### Single SMILES

```bash
smiles2iupac "CC(=O)O"
# acetic acid

smiles2iupac "c1ccccc1"
# benzene

smiles2iupac "C[C@@H](N)C(=O)O"
# D-alanine
```

### Multiple inputs (stdin)

Pass one SMILES per line; one IUPAC name is printed per line.

```bash
printf 'CC\nCCC\nCCCC\n' | smiles2iupac -
# ethane
# propane
# butane

smiles2iupac - < smiles.txt
```

### Help

```bash
smiles2iupac --help
smiles2iupac --version
```

### More examples

```bash
smiles2iupac "CC(C)C"              # 2-methylpropane
smiles2iupac "c1ccc(O)cc1"         # phenol
smiles2iupac "CC(=O)OCC"           # ethyl acetate
smiles2iupac "c1ccc(cc1)C(=O)O"   # benzoic acid
smiles2iupac "ClC(Cl)Cl"           # trichloromethane
```

---

## Library usage

```python
from smiles2iupac import smiles_to_iupac

print(smiles_to_iupac("CC"))                    # ethane
print(smiles_to_iupac("c1ccccc1"))              # benzene
print(smiles_to_iupac("CC(=O)O"))               # acetic acid
print(smiles_to_iupac("C[C@@H](N)C(=O)O"))     # D-alanine
```

---

## Development

```bash
# Run all tests
uv run pytest

# Run tests for a specific phase
uv run pytest tests/test_phase1.py -v
```

### Test organization

Each implementation increment ("phase") historically got its own
`tests/test_phaseNNN.py` file. That history is load-bearing — phase numbers
are referenced from commit messages, code comments, and `docs/REFACTORING_PLAN.md`
— so **existing `test_phaseNNN.py` files are frozen**: don't merge, rename,
or move tests out of them.

New tests that aren't tied to a specific implementation phase (e.g. a
regression test discovered during refactoring, or a cross-cutting property
test) go in `tests/domains/test_<domain>.py` instead (e.g.
`tests/domains/test_fused_thiones.py`). Create the `tests/domains/` directory
the first time it's needed.

---

## License

MIT
