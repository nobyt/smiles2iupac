"""
スピロ化合物・架橋二環式化合物（von Baeyer）の検出と命名。

IUPAC 2013 Blue Book:
  P-31.1.3.1  spiro[m.n]alkane
  P-31.1.3.2  bicyclo[l.m.n]alkane (von Baeyer)
"""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .molecule_analyzer import MoleculeGraph

from .constants import CHAIN_PREFIX, MULTIPLIER


def _ring_neighbors(graph: "MoleculeGraph") -> dict[int, list[int]]:
    """環原子インデックス → 環内隣接原子リスト を返す。"""
    from .molecule_analyzer import get_atom
    ring_set = {a.idx for a in graph.atoms if a.in_ring}
    result: dict[int, list[int]] = {}
    for idx in ring_set:
        result[idx] = [nb for nb in graph.adjacency[idx] if nb in ring_set]
    return result


def _connected_components(
    atom_set: set[int],
    ring_nbrs: dict[int, list[int]],
    exclude: set[int],
) -> list[set[int]]:
    """
    atom_set 内の原子を ring_nbrs で接続し、exclude を除いた
    連結成分のリストを返す。
    """
    remaining = atom_set - exclude
    visited: set[int] = set()
    components: list[set[int]] = []

    for start in remaining:
        if start in visited:
            continue
        component: set[int] = set()
        stack = [start]
        while stack:
            node = stack.pop()
            if node in visited:
                continue
            visited.add(node)
            component.add(node)
            for nb in ring_nbrs.get(node, []):
                if nb not in visited and nb not in exclude:
                    stack.append(nb)
        components.append(component)

    return components


def _alkane_name(n: int) -> str:
    """炭素数 n のアルカン名を返す（例: 7 → heptane）。"""
    stem = CHAIN_PREFIX.get(n, f"C{n}")
    return f"{stem}ane"


_HET_A_SYM = {"N": "aza", "O": "oxa", "S": "thia"}
_HET_A_PRI = {"O": 0, "S": 1, "N": 2}
_HET_A_MULT = {1: "", 2: "di", 3: "tri", 4: "tetra"}


def _het_a_prefix(graph: "MoleculeGraph", het_atoms: set[int], lmap: dict[int, int]) -> str:
    """ヘテロ原子置換接頭辞を生成: e.g. {N:1,N:4} → '1,4-diaza'."""
    from .molecule_analyzer import get_atom
    from collections import defaultdict
    by_sym: dict[str, list[int]] = defaultdict(list)
    for idx in het_atoms:
        by_sym[get_atom(graph, idx).symbol].append(lmap[idx])
    parts: list[str] = []
    for sym in sorted(by_sym, key=lambda s: _HET_A_PRI.get(s, 99)):
        locs = sorted(by_sym[sym])
        loc_str = ",".join(str(x) for x in locs)
        mult = _HET_A_MULT.get(len(locs), str(len(locs)))
        parts.append(f"{loc_str}-{mult}{_HET_A_SYM.get(sym, sym.lower())}")
    return "-".join(parts)


def _bicyclo_unsaturated_name(l: int, m: int, n: int, total: int, db_locs: list[int]) -> str:
    """bicyclo[l.m.n]alkane/ene/diene ベース名を返す。"""
    bracket = f"[{l}.{m}.{n}]"
    stem = CHAIN_PREFIX.get(total, f"C{total}")
    if not db_locs:
        return f"bicyclo{bracket}{stem}ane"
    if len(db_locs) == 1:
        return f"bicyclo{bracket}{stem}-{db_locs[0]}-ene"
    loc_str = ",".join(str(x) for x in db_locs)
    mult = {2: "diene", 3: "triene"}.get(len(db_locs), f"{len(db_locs)}ene")
    return f"bicyclo{bracket}{stem}a-{loc_str}-{mult}"


def _trace_bridge_path(
    start: int, end: int, bridge_atoms: set[int],
    ring_nbrs: dict[int, list[int]],
) -> list[int]:
    """
    start から end へ bridge_atoms を通る経路の原子リスト（start/end 除く）を返す。
    ゼロ橋（bridge_atoms={}）の場合は [] を返す。
    """
    if not bridge_atoms:
        return []
    # bridge_atoms 内で start に隣接する原子を探す
    first = next((nb for nb in ring_nbrs[start] if nb in bridge_atoms), None)
    if first is None:
        return []
    path = [first]
    prev = start
    curr = first
    while curr != end:
        nxt = next(
            (nb for nb in ring_nbrs[curr]
             if nb != prev and (nb in bridge_atoms or nb == end)),
            None,
        )
        if nxt is None:
            break
        prev, curr = curr, nxt
        if curr != end:
            path.append(curr)
    return path


def _collect_polycyclic_subs(
    graph: "MoleculeGraph",
    ring_atoms: set[int],
    locant_map: dict[int, int],
) -> list[tuple[int, str]]:
    """環外置換基を収集して (ロカント, 名前) のリストを返す。"""
    from .molecule_analyzer import get_atom
    from .substituent import name_substituent

    subs: list[tuple[int, str]] = []
    for ring_idx, locant in locant_map.items():
        for nb in graph.adjacency[ring_idx]:
            if nb in ring_atoms:
                continue
            if get_atom(graph, nb).symbol == "H":
                continue
            sub_name = name_substituent(graph, nb, ring_atoms)
            subs.append((locant, sub_name))
    return sorted(subs)


_PCG_TO_SUFFIX: dict[str, str] = {
    "hydroxy": "ol",
    "oxo":     "one",
    "amino":   "amine",
}
# PCGs that keep the terminal "e" before the suffix (no "ane"→"an-" elision)
_PCG_KEEPS_E: set[str] = {"carboxy"}
_PCG_TO_SUFFIX_KEEPE: dict[str, str] = {
    "carboxy": "carboxylic acid",
}
_ALL_PCG: set[str] = set(_PCG_TO_SUFFIX) | _PCG_KEEPS_E


def _format_polycyclic_name(base: str, subs: list[tuple[int, str]]) -> str:
    """置換基プレフィクスを付与した完全名を返す。

    単一の主要特性基 (hydroxy/oxo/amino/carboxy) は接尾辞形式に変換する:
      bicyclo[2.2.2]octane + (2, hydroxy) → bicyclo[2.2.2]octan-2-ol
      bicyclo[2.2.2]octane + (2, carboxy) → bicyclo[2.2.2]octane-2-carboxylic acid
    """
    if not subs:
        return base

    from collections import defaultdict

    # Separate PCG from regular substituents
    pcg_subs = [(loc, nm) for loc, nm in subs if nm in _ALL_PCG]
    reg_subs = [(loc, nm) for loc, nm in subs if nm not in _ALL_PCG]

    # Apply single PCG as a suffix (works for "-ane" and "-ene"/"-diene"/"-triene")
    if len(pcg_subs) == 1:
        loc_pcg, nm_pcg = pcg_subs[0]
        new_base: str | None = None
        if nm_pcg in _PCG_TO_SUFFIX:
            suffix = _PCG_TO_SUFFIX[nm_pcg]
            if base.endswith("ane"):
                stem = base[:-3]
                new_base = f"{stem}an-{loc_pcg}-{suffix}"
            elif re.search(r"(a-\d.*ene|ene)$", base):
                # e.g. "hept-2-ene" → "hept-2-en-{loc}-{suffix}"
                # strip trailing "e" of "ene" → add "-{loc}-{suffix}"
                new_base = base[:-1] + f"-{loc_pcg}-{suffix}"
        else:  # carboxy etc — keep terminal "e"
            suffix = _PCG_TO_SUFFIX_KEEPE[nm_pcg]
            if base.endswith("ane") or re.search(r"ene$", base):
                new_base = f"{base}-{loc_pcg}-{suffix}"
        if new_base is not None:
            return _format_polycyclic_name(new_base, reg_subs)

    grouped: dict[str, list[int]] = defaultdict(list)
    for locant, name in subs:
        grouped[name].append(locant)

    def _needs_parens(nm: str) -> bool:
        if re.search(r"[0-9]", nm):
            return True
        for pfx in ("di", "tri", "tetra", "penta",
                    "fluoro", "chloro", "bromo", "iodo",
                    "hydroxy", "amino", "oxo", "thio", "oxy"):
            if nm.startswith(pfx) and len(nm) > len(pfx):
                return True
        return False

    parts: list[str] = []
    for name in sorted(grouped):
        locs = sorted(grouped[name])
        loc_str = ",".join(str(l) for l in locs)
        n = len(locs)
        mult = MULTIPLIER.get(n, f"{n}") if n > 1 else ""
        parens = _needs_parens(name)
        parts.append(
            f"{loc_str}-{mult}({name})" if parens else f"{loc_str}-{mult}{name}"
        )

    def _alpha_key(s: str) -> str:
        return re.sub(r"^[\d,\-]+", "", s)

    parts.sort(key=_alpha_key)
    sep = "-" if base[:1].isdigit() else ""
    prefix = "-".join(parts[:-1]) + "-" + parts[-1] if len(parts) > 1 else parts[0]
    return f"{prefix}{sep}{base}"


# ─── スピロ化合物 ─────────────────────────────────────────────────────────────

def _try_spiro(graph: "MoleculeGraph") -> str | None:
    """
    スピロ化合物なら 'spiro[m.n]alkane' (ヘテロ原子・置換基付きも含む) を返す。
    スピロ原子 = 環結合数 4 の炭素原子がちょうど 1 つ存在する非芳香族環。
    """
    from .molecule_analyzer import get_atom

    ring_nbrs = _ring_neighbors(graph)
    ring_atoms = set(ring_nbrs.keys())

    # 非芳香族の C/N/O/S のみ対応
    _ALLOWED = frozenset({"C", "N", "O", "S"})
    for idx in ring_atoms:
        a = get_atom(graph, idx)
        if a.symbol not in _ALLOWED or a.is_aromatic:
            return None

    het_atoms = {idx for idx in ring_atoms if get_atom(graph, idx).symbol != "C"}

    deg4 = [idx for idx, nbs in ring_nbrs.items() if len(nbs) == 4]
    deg3 = [idx for idx, nbs in ring_nbrs.items() if len(nbs) == 3]

    if len(deg4) != 1 or deg3:
        return None

    spiro_atom = deg4[0]
    if get_atom(graph, spiro_atom).symbol != "C":
        return None

    components = _connected_components(ring_atoms, ring_nbrs, exclude={spiro_atom})

    if len(components) != 2:
        return None

    comp_list = list(components)
    comp_list.sort(key=len)  # smaller first
    sizes = [len(c) for c in comp_list]
    m, n = sizes[0], sizes[1]
    total = m + n + 1
    ring_name_base = f"spiro[{m}.{n}]{_alkane_name(total)}"

    heavy = {a.idx for a in graph.atoms if a.symbol not in ("H",)}
    if not (heavy - ring_atoms) and not het_atoms:
        return ring_name_base

    def _traverse_ring_from(start: int, comp: set[int]) -> list[int]:
        path = [start]
        prev = spiro_atom
        curr = start
        while True:
            nxt = next(
                (nb for nb in ring_nbrs[curr]
                 if nb != prev and (nb in comp or nb == spiro_atom)),
                None,
            )
            if nxt is None or nxt == spiro_atom:
                break
            path.append(nxt)
            prev, curr = curr, nxt
        return path

    spiro_nbs = ring_nbrs[spiro_atom]

    best_subs: list[tuple[int, str]] | None = None
    best_lmap: dict[int, int] | None = None
    best_score: tuple | None = None

    # 同サイズの場合は両コンポーネントを small/large 両方で試す
    assignments: list[tuple[set[int], set[int]]] = []
    if m == n:
        assignments = [(comp_list[0], comp_list[1]), (comp_list[1], comp_list[0])]
    else:
        assignments = [(comp_list[0], comp_list[1])]

    for small_comp, large_comp in assignments:
        small_starts = [nb for nb in spiro_nbs if nb in small_comp]
        large_starts = [nb for nb in spiro_nbs if nb in large_comp]

        for s_start in small_starts:
            for l_start in large_starts:
                lmap: dict[int, int] = {}
                pos = 1

                small_path = _traverse_ring_from(s_start, small_comp)
                large_path = _traverse_ring_from(l_start, large_comp)

                for a in small_path:
                    lmap[a] = pos; pos += 1
                lmap[spiro_atom] = pos; pos += 1
                for a in large_path:
                    lmap[a] = pos; pos += 1

                if len(lmap) != len(ring_atoms):
                    continue

                het_locs = sorted(lmap[idx] for idx in het_atoms)
                subs = _collect_polycyclic_subs(graph, ring_atoms, lmap)
                sub_locs = sorted(loc for loc, _ in subs)
                score = (het_locs, sub_locs)

                if best_score is None or score < best_score:
                    best_score = score
                    best_subs = subs
                    best_lmap = lmap

    if best_lmap is None:
        return ring_name_base if not het_atoms else None

    het_prefix = _het_a_prefix(graph, het_atoms, best_lmap) if het_atoms else ""
    final_base = het_prefix + ring_name_base
    return _format_polycyclic_name(final_base, best_subs or [])


# ─── 架橋二環式化合物 ─────────────────────────────────────────────────────────

def _try_bicyclo(graph: "MoleculeGraph") -> str | None:
    """
    von Baeyer 式二環化合物なら 'bicyclo[l.m.n]alkane' を返す。
    架橋頭炭素 = 環結合数 3 の原子がちょうど 2 つ存在する純炭素環。
    置換基付きも対応。
    """
    from .molecule_analyzer import get_atom
    from itertools import permutations

    ring_nbrs = _ring_neighbors(graph)
    ring_atoms = set(ring_nbrs.keys())

    # 非芳香族の C/N/O/S のみ対応
    _ALLOWED = frozenset({"C", "N", "O", "S"})
    for idx in ring_atoms:
        a = get_atom(graph, idx)
        if a.symbol not in _ALLOWED or a.is_aromatic:
            return None

    # ヘテロ原子セット
    het_atoms = {idx for idx in ring_atoms if get_atom(graph, idx).symbol != "C"}

    deg3 = [idx for idx, nbs in ring_nbrs.items() if len(nbs) == 3]
    deg4 = [idx for idx, nbs in ring_nbrs.items() if len(nbs) >= 4]

    if len(deg3) != 2 or deg4:
        return None

    b1, b2 = deg3
    components = _connected_components(ring_atoms, ring_nbrs, exclude={b1, b2})

    bridge_lengths = sorted((len(c) for c in components), reverse=True)

    has_zero = b2 in ring_nbrs[b1]
    if has_zero:
        bridge_lengths.append(0)
        bridge_lengths.sort(reverse=True)
    if len(bridge_lengths) != 3:
        return None

    l, m, n = bridge_lengths
    total = l + m + n + 2

    # Phase 276: 環内二重結合を検出
    from .molecule_analyzer import get_bond_order as _gbo
    ring_db_pairs: list[tuple[int, int]] = []
    for a in sorted(ring_atoms):
        for b in ring_nbrs[a]:
            if b > a and _gbo(graph, a, b) == 2.0:
                ring_db_pairs.append((a, b))

    # 環外置換基チェック
    heavy = {a.idx for a in graph.atoms if a.symbol not in ("H",)}
    has_subs = bool(heavy - ring_atoms)

    if not has_subs and not ring_db_pairs and not het_atoms:
        return f"bicyclo[{l}.{m}.{n}]{_alkane_name(total)}"

    # ─ ロカント付与 ─
    # IUPAC P-31.1.2.2: bh1=1, bridge1(l atoms)=2..l+1, bh2=l+2,
    #   bridge2(m atoms from bh2)=l+3..l+m+2,
    #   bridge3(n atoms from bh1)=l+m+3..l+m+n+2
    comp_sets: list[set[int]] = list(components)
    if has_zero:
        comp_sets.append(set())
    # Sort by size descending
    comp_sets.sort(key=len, reverse=True)

    best_subs: list[tuple[int, str]] | None = None
    best_lmap: dict[int, int] | None = None
    best_db_locs: list[int] = []
    best_score: tuple | None = None

    for bh1, bh2 in [(b1, b2), (b2, b1)]:
        for perm in permutations(comp_sets):
            sizes = [len(c) for c in perm]
            if sizes != sorted(sizes, reverse=True):
                continue  # 降順でない → 無効

            # build locant map
            lmap: dict[int, int] = {}
            lmap[bh1] = 1
            pos = 2

            # bridge1: from bh1 toward bh2
            path1 = _trace_bridge_path(bh1, bh2, perm[0], ring_nbrs)
            for a in path1:
                lmap[a] = pos; pos += 1
            lmap[bh2] = pos; pos += 1

            # bridge2: from bh2 toward bh1
            path2 = _trace_bridge_path(bh2, bh1, perm[1], ring_nbrs)
            for a in path2:
                lmap[a] = pos; pos += 1

            # bridge3: from bh1 toward bh2
            path3 = _trace_bridge_path(bh1, bh2, perm[2], ring_nbrs)
            for a in path3:
                lmap[a] = pos; pos += 1

            if len(lmap) != len(ring_atoms):
                continue

            het_locs = sorted(lmap[idx] for idx in het_atoms)
            db_locs = sorted(min(lmap[a], lmap[b]) for a, b in ring_db_pairs)
            subs = _collect_polycyclic_subs(graph, ring_atoms, lmap)
            # Phase 287/289: PCG ロカントを二重結合ロカントより優先
            pcg_locs = sorted(loc for loc, nm in subs if nm in _ALL_PCG)
            reg_locs = sorted(loc for loc, nm in subs if nm not in _ALL_PCG)
            # Priority: het > PCG > double bonds > other subs
            score = (het_locs, pcg_locs, db_locs, reg_locs)

            if best_score is None or score < best_score:
                best_score = score
                best_subs = subs
                best_lmap = lmap
                best_db_locs = db_locs

    if best_lmap is None:
        return f"bicyclo[{l}.{m}.{n}]{_alkane_name(total)}"

    het_prefix = _het_a_prefix(graph, het_atoms, best_lmap) if het_atoms else ""
    final_base = het_prefix + _bicyclo_unsaturated_name(l, m, n, total, best_db_locs)
    return _format_polycyclic_name(final_base, best_subs or [])


_CAGE_REF_SEED: dict[str, str] = {
    "cubane": "C12C3C4C1C5C2C3C45",
    "adamantane": "C1C2CC3CC1CC(C2)C3",
}
_cage_ref_canon: dict[str, str] = {}


def _cage_ref(name: str) -> str | None:
    """ケージ保留名の炭素骨格の標準 (canonical) SMILES を返す (キャッシュ)。"""
    if name not in _cage_ref_canon:
        from rdkit import Chem
        m = Chem.MolFromSmiles(_CAGE_REF_SEED[name])
        _cage_ref_canon[name] = Chem.MolToSmiles(m) if m is not None else None
    return _cage_ref_canon[name]


def _ring_skeleton_canon(graph: "MoleculeGraph", ring_atoms) -> str | None:
    """環原子だけからなる炭化水素骨格を組み立て、その標準 SMILES を返す。

    次数分布が同一でも位相が異なる縮退ケージ (twistane vs adamantane,
    cuneane vs cubane 等) を確実に区別するために使う。
    """
    from rdkit import Chem

    rw = Chem.RWMol()
    idx_map: dict[int, int] = {}
    for gi in ring_atoms:
        idx_map[gi] = rw.AddAtom(Chem.Atom(6))
    for gi in ring_atoms:
        for nb in graph.adjacency[gi]:
            if nb in idx_map and gi < nb:
                rw.AddBond(idx_map[gi], idx_map[nb], Chem.BondType.SINGLE)
    m = rw.GetMol()
    try:
        Chem.SanitizeMol(m)
    except Exception:
        return None
    return Chem.MolToSmiles(m)


def _try_cage_retained(graph: "MoleculeGraph") -> str | None:
    """
    adamantane / cubane などのケージ化合物保留名を返す。
    単置換も対応: 1-X (bridgehead) または 2-X (bridge) adamantane。
    """
    from .molecule_analyzer import get_atom

    ring_nbrs = _ring_neighbors(graph)
    ring_atoms = set(ring_nbrs.keys())

    # 純炭化水素・非芳香族のみ
    for idx in ring_atoms:
        a = get_atom(graph, idx)
        if a.symbol != "C" or a.is_aromatic:
            return None

    n = len(ring_atoms)
    deg_counts: dict[int, int] = {}
    for nbs in ring_nbrs.values():
        d = len(nbs)
        deg_counts[d] = deg_counts.get(d, 0) + 1

    # Cubane: 8 atoms, all deg-3 (次数だけでは cuneane と区別できないため骨格照合)
    if n == 8 and deg_counts.get(3) == 8 and len(deg_counts) == 1:
        heavy = {a.idx for a in graph.atoms if a.symbol not in ("H",)}
        if not (heavy - ring_atoms) and _ring_skeleton_canon(graph, ring_atoms) == _cage_ref("cubane"):
            return "cubane"

    # Adamantane: 10 atoms, 4×deg-3 (bridgehead) + 6×deg-2 (bridge)
    if not (n == 10 and deg_counts.get(3) == 4 and deg_counts.get(2) == 6):
        return None

    # Phase 947: 次数分布が同一でも位相の異なる縮退 C10 ケージ (twistane 等) を
    # adamantane と誤認しないよう、環骨格の標準 SMILES を実際の adamantane と照合する。
    if _ring_skeleton_canon(graph, ring_atoms) != _cage_ref("adamantane"):
        return None

    heavy = {a.idx for a in graph.atoms if a.symbol not in ("H",)}
    sub_anchors = heavy - ring_atoms
    if not sub_anchors:
        return "adamantane"

    # 置換アダマンタン: 1種類の置換基 + 1箇所のみ対応
    bridgeheads = {idx for idx, nbs in ring_nbrs.items() if len(nbs) == 3}
    bridges = {idx for idx, nbs in ring_nbrs.items() if len(nbs) == 2}

    subs: list[tuple[int, str]] = []
    from .substituent import name_substituent
    for ring_idx in ring_atoms:
        for nb in graph.adjacency[ring_idx]:
            if nb in ring_atoms or get_atom(graph, nb).symbol == "H":
                continue
            sub_name = name_substituent(graph, nb, ring_atoms)
            loc = 1 if ring_idx in bridgeheads else 2
            subs.append((loc, sub_name))
    if not subs:
        return "adamantane"

    from .name_assembler import _build_prefix
    from collections import defaultdict
    grouped: dict[str, list[int]] = defaultdict(list)
    for loc, nm in subs:
        grouped[nm].append(loc)

    def _needs_p(nm: str) -> bool:
        if re.search(r"[0-9]", nm):
            return True
        for pfx in ("di", "tri", "tetra", "penta",
                    "fluoro", "chloro", "bromo", "iodo",
                    "hydroxy", "amino", "oxo", "thio", "oxy"):
            if nm.startswith(pfx) and len(nm) > len(pfx):
                return True
        return False

    parts: list[str] = []
    for nm in sorted(grouped):
        locs = sorted(set(grouped[nm]))
        loc_str = ",".join(str(l) for l in locs)
        cnt = len(subs)
        mult = MULTIPLIER.get(cnt, str(cnt)) if cnt > 1 else ""
        parts.append(
            f"{loc_str}-{mult}({nm})" if _needs_p(nm) else f"{loc_str}-{mult}{nm}"
        )

    def _alpha_key(s: str) -> str:
        return re.sub(r"^[\d,\-]+", "", s)

    parts.sort(key=_alpha_key)
    prefix = "-".join(parts)
    return f"{prefix}adamantane"


_SUP = {
    "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴",
    "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹",
    ",": ",",
}


def _to_sup(s: str) -> str:
    return "".join(_SUP.get(c, c) for c in s)


_POLY_RING_PREFIX = {
    3: "tricyclo",
    4: "tetracyclo",
    5: "pentacyclo",
    6: "hexacyclo",
    7: "heptacyclo",
    8: "octacyclo",
}


def _ring_subgraph_has_cut_vertex(ring_nbrs: dict[int, list[int]]) -> bool:
    """環部分グラフに関節点 (cut vertex / articulation point) が存在するか判定する。

    純粋な von Baeyer 縮合/架橋多環系は 2-連結 (biconnected) で関節点を持たない。
    関節点があれば、その原子はスピロ結合点 (両側の環系を1原子だけで繋ぐ) であり、
    von Baeyer nomenclature では正しく表現できない。
    """
    nodes = list(ring_nbrs.keys())
    if not nodes:
        return False
    disc: dict[int, int] = {}
    low: dict[int, int] = {}
    timer = [0]
    cut = [False]

    import sys
    old_limit = sys.getrecursionlimit()
    sys.setrecursionlimit(max(old_limit, len(nodes) * 4 + 100))

    def _dfs(u: int, parent: int) -> None:
        disc[u] = low[u] = timer[0]
        timer[0] += 1
        children = 0
        for v in ring_nbrs[u]:
            if v not in disc:
                children += 1
                _dfs(v, u)
                low[u] = min(low[u], low[v])
                if parent != -1 and low[v] >= disc[u]:
                    cut[0] = True
                if parent == -1 and children > 1:
                    cut[0] = True
            elif v != parent:
                low[u] = min(low[u], disc[v])

    try:
        _dfs(nodes[0], -1)
    finally:
        sys.setrecursionlimit(old_limit)
    return cut[0]


def _try_polycyclic_von_baeyer(graph: "MoleculeGraph") -> str | None:
    """
    3環以上の架橋多環式炭化水素 (von Baeyer: tricyclo[...], tetracyclo[...]) を命名する (IUPAC 2013 P-23.2.3 - P-23.2.5)。
    """
    from .molecule_analyzer import get_atom
    import itertools

    ring_nbrs = _ring_neighbors(graph)
    ring_atoms = set(ring_nbrs.keys())
    if not ring_atoms:
        return None

    # 全て炭素かつ非芳香族
    for idx in ring_atoms:
        a = get_atom(graph, idx)
        if a.symbol != "C" or a.is_aromatic:
            return None

    # 環外置換基のない純炭化水素を対象とする
    heavy = {a.idx for a in graph.atoms if a.symbol != "H"}
    if heavy != ring_atoms:
        return None

    n_atoms = len(ring_atoms)
    # 環内の辺の数をカウント
    n_bonds = sum(len(nbs) for nbs in ring_nbrs.values()) // 2
    n_rings = n_bonds - n_atoms + 1

    if n_rings < 3:
        return None

    # Phase 948: スピロ原子 (環部分グラフの関節点) を含む系は 2-連結でなく、
    # 純粋な von Baeyer 縮合/架橋多環ではない。誤った tricyclo 記述子
    # (橋の本数が環数と矛盾する) を生成しないようここで撤退し、上位の
    # 「未対応」ガード / スピロハンドラに委ねる。
    if _ring_subgraph_has_cut_vertex(ring_nbrs):
        return None

    ring_prefix = _POLY_RING_PREFIX.get(n_rings, f"{n_rings}-cyclo")

    # 全ての単純閉路を探索
    def _find_all_cycles() -> list[list[int]]:
        cycles = []
        visited_cycles = set()
        nodes = sorted(ring_atoms)

        def _dfs(start, curr, path, visited):
            for nb in ring_nbrs[curr]:
                if nb == start and len(path) >= 3:
                    cycle = list(path)
                    # 正規化
                    min_idx = cycle.index(min(cycle))
                    rot = cycle[min_idx:] + cycle[:min_idx]
                    if rot[1] > rot[-1]:
                        rot = [rot[0]] + list(reversed(rot[1:]))
                    t_rot = tuple(rot)
                    if t_rot not in visited_cycles:
                        visited_cycles.add(t_rot)
                        cycles.append(rot)
                elif nb not in visited and nb > start:
                    visited.add(nb)
                    _dfs(start, nb, path + [nb], visited)
                    visited.remove(nb)

        for n in nodes:
            _dfs(n, n, [n], {n})
        return cycles

    all_cycles = _find_all_cycles()
    if not all_cycles:
        return None

    # 主環（最大長）の候補
    max_c_len = max(len(c) for c in all_cycles)
    candidate_main_rings = [c for c in all_cycles if len(c) == max_c_len]

    best_candidate: tuple | None = None

    for main_ring in candidate_main_rings:
        main_ring_set = set(main_ring)
        mr_len = len(main_ring)
        # 主環内の辺集合
        main_edges = set()
        for i in range(mr_len):
            u, v = main_ring[i], main_ring[(i + 1) % mr_len]
            main_edges.add((min(u, v), max(u, v)))

        # 主環上のペア (b1, b2) を探索
        for i in range(mr_len):
            for j in range(i + 1, mr_len):
                b1, b2 = main_ring[i], main_ring[j]
                # 主環を 2 つの経路に分割
                # 経路 1: i -> j
                p1_nodes = main_ring[i:j + 1]
                # 経路 2: j -> i (環を回る)
                p2_nodes = main_ring[j:] + main_ring[:i + 1]
                l = len(p1_nodes) - 2
                m = len(p2_nodes) - 2
                if m > l:
                    l, m = m, l
                    p1_nodes, p2_nodes = p2_nodes, p1_nodes

                # 主環外を通る b1-b2 経路 (主橋) を探索
                def _find_bridges_between(start, end, forbidden_edges):
                    paths = []
                    def _bfs():
                        queue = [(start, [start], {start})]
                        while queue:
                            curr, path, vis = queue.pop(0)
                            if curr == end and len(path) >= 2:
                                paths.append(path)
                                continue
                            for nb in ring_nbrs[curr]:
                                edge = (min(curr, nb), max(curr, nb))
                                if edge in forbidden_edges:
                                    continue
                                if nb not in vis:
                                    queue.append((nb, path + [nb], vis | {nb}))
                    _bfs()
                    return paths

                bridge_paths = _find_bridges_between(b1, b2, main_edges)
                if not bridge_paths:
                    continue

                for mb_path in bridge_paths:
                    n = len(mb_path) - 2
                    if n > m:  # l >= m >= n の順序
                        continue
                    # 主二環系に含まれる辺
                    mb_edges = set()
                    for k in range(len(mb_path) - 1):
                        u, v = mb_path[k], mb_path[k + 1]
                        mb_edges.add((min(u, v), max(u, v)))

                    bicycle_edges = main_edges | mb_edges
                    bicycle_nodes = set(main_ring) | set(mb_path)

                    # 残りの辺 (副橋)
                    all_edges = {
                        (min(u, v), max(u, v))
                        for u in ring_atoms for v in ring_nbrs[u]
                    }
                    rem_edges = all_edges - bicycle_edges

                    # 2方向での番号付け
                    for bh1, bh2 in [(b1, b2), (b2, b1)]:
                        # bh1 から出発して最長経路 l を通り bh2 へ、次に m を通り bh1 へ、主橋 n を通る
                        path_l = p1_nodes if p1_nodes[0] == bh1 else list(reversed(p1_nodes))
                        path_m = p2_nodes if p2_nodes[0] == bh2 else list(reversed(p2_nodes))
                        path_n = mb_path if mb_path[0] == bh1 else list(reversed(mb_path))

                        lmap: dict[int, int] = {}
                        lmap[bh1] = 1
                        pos = 2
                        for a in path_l[1:-1]:
                            lmap[a] = pos; pos += 1
                        lmap[bh2] = pos; pos += 1
                        for a in path_m[1:-1]:
                            lmap[a] = pos; pos += 1
                        for a in path_n[1:-1]:
                            lmap[a] = pos; pos += 1

                        # 残りの原子（副橋の内部原子）に番号付け
                        rem_nodes = ring_atoms - set(lmap.keys())
                        for a in sorted(rem_nodes):
                            lmap[a] = pos; pos += 1

                        # 副橋の同定
                        # rem_edges から連結成分または単一辺を抽出
                        sec_bridges: list[tuple[int, int, int]] = []
                        visited_rem = set()

                        for edge in sorted(rem_edges):
                            if edge in visited_rem:
                                continue
                            visited_rem.add(edge)
                            u, v = edge
                            # 直接結合の場合
                            if u in lmap and v in lmap and abs(lmap[u] - lmap[v]) > 0:
                                x, y = sorted([lmap[u], lmap[v]])
                                sec_bridges.append((0, x, y))

                        # Phase 950: この候補が要求される副橋の本数
                        # (n_rings - 2; tricyclo=1本, tetracyclo=2本, ...) と
                        # 一致しない場合は不正な記述子 (例: tricyclo なのに
                        # 副橋が2本) になるため、この候補を棄却する。
                        # 以前はここが no-op (`pass`) で、不一致を検出しながら
                        # 何もせず不正な記述子をそのまま採用していた。
                        # rem_edges に複数原子からなる副橋 (直接辺でない) が
                        # 含まれる場合も現状 sec_bridges に反映されないため
                        # 同様に棄却対象となる。
                        if len(sec_bridges) != n_rings - 2:
                            continue

                        # ソート: 長さ降順、x昇順、y昇順
                        sec_bridges.sort(key=lambda item: (-item[0], item[1], item[2]))

                        sec_locs_key = tuple((item[0], item[1], item[2]) for item in sec_bridges)
                        score = (-mr_len, -l, -m, -n, sec_locs_key)

                        if best_candidate is None or score < best_candidate[0]:
                            best_candidate = (score, l, m, n, sec_bridges, n_atoms)

    if best_candidate is None:
        return None

    _, l, m, n, sec_bridges, total_c = best_candidate
    alkane_name = _alkane_name(total_c)

    sec_parts = []
    for d, x, y in sec_bridges:
        sup_str = _to_sup(f"{x},{y}")
        sec_parts.append(f".{d}{sup_str}")

    sec_str = "".join(sec_parts)
    return f"{ring_prefix}[{l}.{m}.{n}{sec_str}]{alkane_name}"


def name_polycyclic(graph: "MoleculeGraph") -> str | None:
    """
    スピロ・架橋二環式・ケージ化合物・多環式 von Baeyer なら IUPAC 名を返す。
    対象外の場合は None を返す。
    """
    cage = _try_cage_retained(graph)
    if cage is not None:
        return cage
    result = _try_spiro(graph)
    if result is not None:
        return result
    bicyclic = _try_bicyclo(graph)
    if bicyclic is not None:
        return bicyclic
    return _try_polycyclic_von_baeyer(graph)

