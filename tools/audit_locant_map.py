"""Refactoring plan Step 2: audit `_FUSED_LOCANT_MAP` for None locants on
heteroatoms.

A `None` locant means the atom's substituent is silently dropped from the
assembled name (see Phase 849 / 851 bugs, both caused by a heteroatom being
mapped to `None`). Carbon atoms with `None` are normal and expected (only
atoms that can carry IUPAC-numbered substituents need a locant). Heteroatoms
(N/O/S/Se/Te) mapped to `None` are the risky cases worth reviewing.

Usage:
    uv run python tools/audit_locant_map.py
"""
from __future__ import annotations

from rdkit import Chem

from smiles2iupac.heterocycle_handler import _FUSED_LOCANT_MAP

_HETERO_SYMBOLS = {"N", "O", "S", "Se", "Te"}


def audit() -> list[tuple[str, int, str]]:
    """Return (smiles, atom_idx, element) for every heteroatom mapped to None."""
    flagged: list[tuple[str, int, str]] = []
    for smiles, locant_map in _FUSED_LOCANT_MAP.items():
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            flagged.append((smiles, -1, "UNPARSEABLE"))
            continue
        for atom_idx, locant in locant_map.items():
            if locant is not None:
                continue
            if atom_idx >= mol.GetNumAtoms():
                continue
            atom = mol.GetAtomWithIdx(atom_idx)
            symbol = atom.GetSymbol()
            if symbol in _HETERO_SYMBOLS:
                flagged.append((smiles, atom_idx, symbol))
    return flagged


def main() -> None:
    flagged = audit()
    total_entries = len(_FUSED_LOCANT_MAP)
    total_none = sum(
        1 for m in _FUSED_LOCANT_MAP.values() for v in m.values() if v is None
    )
    print(f"_FUSED_LOCANT_MAP: {total_entries} ring SMILES entries, "
          f"{total_none} total None locants")
    print(f"Heteroatom (N/O/S/Se/Te) None locants: {len(flagged)}\n")
    for smiles, atom_idx, symbol in flagged:
        print(f"  {smiles:45s}  atom {atom_idx:>2}  {symbol}  -> None")


if __name__ == "__main__":
    main()
