"""
ヘテロ環保留名の検出と命名。

IUPAC 2013 Blue Book P-31.1.3.4 / P-31.1.3.5 / P-31.1.3.6

対応範囲:
  - 芳香族ヘテロ環: pyridine, furan, thiophene, 1H-pyrrole, 1H-imidazole,
                    pyrimidine, pyrazine, pyridazine
  - 飽和ヘテロ環: piperidine, morpholine
  - 置換ヘテロ環: 2-methylpyridine, 3-chloropyridine 等
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .molecule_analyzer import MoleculeGraph


# ─── 優先度付きヘテロ原子記号 ─────────────────────────────────────────
# NH (芳香族 NH) は別扱い: 他の N より優先し "NH" として記録する

_PRIORITY: dict[str, int] = {
    "NH": 0,   # 最高優先 (1H-pyrrole, 1H-imidazole)
    "O":  1,
    "S":  2,
    "N":  3,
    "P":  4,
    "Se": 5,
    "Te": 6,
}


def _atom_sig(graph: "MoleculeGraph", idx: int) -> str:
    """原子インデックスをシグネチャ文字列に変換する。"""
    from .molecule_analyzer import get_atom
    a = get_atom(graph, idx)
    if a.symbol == "N" and a.is_aromatic:
        # build_molecule_graph は AddHs() 後にグラフを構築するため
        # num_hs ではなく隣接 H 原子の有無で NH を判定する
        has_h = any(get_atom(graph, nb).symbol == "H" for nb in graph.adjacency[idx])
        if has_h:
            return "NH"
    if a.symbol == "C":
        return "C"
    return a.symbol


def _canonical_sig(ordered: list[int], graph: "MoleculeGraph") -> tuple[str, ...]:
    """
    環の原子リスト（始端=ヘテロ原子）から正準シグネチャを生成する。

    2つの方向（時計回り・反時計回り）の両方を試し、
    辞書順で小さい方を選択する。
    """
    n = len(ordered)
    fwd = tuple(_atom_sig(graph, ordered[i % n]) for i in range(n))
    rev = tuple(_atom_sig(graph, ordered[(-i) % n]) for i in range(n))
    return min(fwd, rev)


def _hetero_priority(sig: str) -> int:
    return _PRIORITY.get(sig, 99)


def _find_best_start(ring: list[int], graph: "MoleculeGraph") -> list[int]:
    """
    ヘテロ原子の中から最高優先度のものを始端にした環の回転を返す。
    複数ある場合はロカント集合が最小になる回転を選択する。
    同点の場合は外環 C=O ロカント集合が最小の回転を選択する（IUPAC P-14.5）。
    """
    from .molecule_analyzer import get_atom, get_bond_order

    n = len(ring)
    ring_set = set(ring)
    # ヘテロ原子の位置を特定
    hetero_positions = [
        i for i, idx in enumerate(ring)
        if _atom_sig(graph, idx) != "C"
    ]
    if not hetero_positions:
        return ring

    # 最高優先ヘテロ原子の種類を特定
    best_priority = min(_hetero_priority(_atom_sig(graph, ring[i])) for i in hetero_positions)
    best_positions = [
        i for i in hetero_positions
        if _hetero_priority(_atom_sig(graph, ring[i])) == best_priority
    ]

    # 候補の回転を生成: 各最高優先ヘテロ原子を位置1に置く (順/逆の両方向を試す)
    candidates: list[tuple[tuple[str, ...], list[int]]] = []
    for start in best_positions:
        fwd = ring[start:] + ring[:start]
        rev = [ring[(start - i) % n] for i in range(n)]
        for rotated in (fwd, rev):
            sig = _canonical_sig(rotated, graph)
            candidates.append((sig, rotated))

    # ロカント集合が最小の回転を選択
    candidates.sort(key=lambda x: x[0])
    best_sig = candidates[0][0]
    best_candidates = [rot for sig, rot in candidates if sig == best_sig]

    def _hetero_locant_key(rot: list[int]) -> list[int]:
        return sorted(i + 1 for i, a in enumerate(rot) if _atom_sig(graph, a) != "C")

    def _exo_co_locant_key(rot: list[int]) -> list[int]:
        locs = []
        for i, idx in enumerate(rot):
            if get_atom(graph, idx).symbol != "C":
                continue
            has_exo_co = any(
                get_atom(graph, nb).symbol == "O"
                and get_bond_order(graph, idx, nb) == 2.0
                for nb in graph.adjacency[idx]
                if nb not in ring_set
            )
            if has_exo_co:
                locs.append(i + 1)
        return sorted(locs)

    def _composite_key(rot: list[int]) -> tuple:
        return (_hetero_locant_key(rot), _exo_co_locant_key(rot))

    best_candidates.sort(key=_composite_key)
    return best_candidates[0]


# ─── 保留名テーブル ────────────────────────────────────────────────────
# key: (is_aromatic, canonical_signature_tuple)
# value: (retained_name, nh_indicator)
#   nh_indicator: True → "1H-" プレフィクスを付ける
#   注: canonical_sig は始端=最高優先ヘテロ原子、辞書順最小方向

# ─── Hantzsch-Widman 名テーブル ───────────────────────────────────────
# key: (ring_size, heteroatom_sig)  — 単一ヘテロ原子・非芳香族のみ
_HW_NAMES: dict[tuple[int, str], str] = {
    (3, "O"): "oxirane",
    (3, "S"): "thiirane",
    (3, "N"): "aziridine",
    (3, "P"): "phosphirane",
    (4, "O"): "oxetane",
    (4, "S"): "thietane",
    (4, "N"): "azetidine",
    (4, "P"): "phosphetane",
    (5, "O"): "oxolane",
    (5, "S"): "thiolane",
    (5, "N"): "pyrrolidine",
    (5, "P"): "phospholane",
    (6, "O"): "oxane",
    (6, "S"): "thiane",
    (6, "P"): "phosphinane",
    # (6, "N") → piperidine は保留名テーブルで処理済み
    (7, "N"): "azepane",
    (7, "O"): "oxepane",
    (7, "S"): "thiepane",
    (7, "P"): "phosphepane",
    (8, "N"): "azocane",
    (8, "O"): "oxocane",
    (8, "S"): "thiocane",
    (8, "P"): "phosphocane",
}


def _match_hantzsch_widman(ring: list[int], graph: "MoleculeGraph") -> str | None:
    """
    単一ヘテロ原子の飽和環を Hantzsch-Widman 名で返す。
    一致しない場合は None。
    """
    if _is_aromatic_ring(ring, graph):
        return None
    n = len(ring)
    if n not in (3, 4, 5, 6, 7, 8):
        return None
    # ヘテロ原子を収集
    hetero_sigs: list[str] = []
    for idx in ring:
        sig = _atom_sig(graph, idx)
        if sig != "C":
            hetero_sigs.append(sig)
    # 単一ヘテロ原子のみ対応
    if len(hetero_sigs) != 1:
        return None
    return _HW_NAMES.get((n, hetero_sigs[0]))


# 多ヘテロ原子 HW 名: 7-10員環 (ring_size → suffix without heteroatom prefix)
_HW_MULTI_SUFFIX: dict[int, str] = {
    7: "epane", 8: "ocane", 9: "onane", 10: "ecane",
}
# ヘテロ原子ベースプレフィックス (末尾 'a' を除いた形)
_HW_HET_BASE: dict[str, str] = {"O": "ox", "N": "az", "S": "thi"}
_HW_HET_PRIO: dict[str, int] = {"O": 0, "S": 1, "N": 2}
_HW_MULT: dict[int, str] = {1: "", 2: "di", 3: "tri", 4: "tetra"}


def _match_multi_het_ring(ring: list[int], graph: "MoleculeGraph") -> str | None:
    """
    複数ヘテロ原子を持つ非芳香族環 (7-10員) を a-命名法で返す。
    例: 7員環 2O(1,4) → '1,4-dioxepane', 7員環 N+O → '1-oxa-4-azepane'
    """
    from collections import defaultdict

    if _is_aromatic_ring(ring, graph):
        return None
    n = len(ring)
    suffix = _HW_MULTI_SUFFIX.get(n)
    if suffix is None:
        return None

    het_indices = [idx for idx in ring if _atom_sig(graph, idx) != "C"]
    if len(het_indices) < 2:
        return None

    syms = [_atom_sig(graph, idx) for idx in het_indices]
    if not all(s in _HW_HET_BASE for s in syms):
        return None

    rotation = _find_best_start(ring, graph)
    lmap = {idx: i + 1 for i, idx in enumerate(rotation)}

    by_sym: dict[str, list[int]] = defaultdict(list)
    for idx in het_indices:
        sym = _atom_sig(graph, idx)
        by_sym[sym].append(lmap[idx])

    sorted_syms = sorted(by_sym.keys(), key=lambda s: _HW_HET_PRIO.get(s, 99))
    parts: list[str] = []
    for i, sym in enumerate(sorted_syms):
        locs = sorted(by_sym[sym])
        loc_str = ",".join(str(x) for x in locs)
        mult = _HW_MULT.get(len(locs), str(len(locs)))
        base = _HW_HET_BASE[sym]
        if i == len(sorted_syms) - 1:
            parts.append(f"{loc_str}-{mult}{base}{suffix}")
        else:
            parts.append(f"{loc_str}-{mult}{base}a")

    return "-".join(parts)


# ─── Phase 272: 部分不飽和単環ヘテロ環 ─────────────────────────────────
# 非芳香族で環内二重結合を持つ 5/6員環 → "X,Y-dihydro[parent]" 形式

_PARTIAL_UNSAT_PARENT: dict[tuple[int, str], str] = {
    # 5-membered
    (5, "O"): "furan",
    (5, "S"): "thiophene",
    (5, "N"): "pyrrole",
    # 6-membered
    (6, "N"): "pyridine",
    (6, "O"): "pyran",
    (6, "S"): "thiopyran",
}

# Phase 508: 2-ヘテロ原子5員部分不飽和環の親芳香族名
# key = 最低ロカント方向で起算したシグネチャタプル
# value = 親芳香族名
_DIHETERO_PARENT_BY_LOWLOC_SIG: dict[tuple[str, ...], str] = {
    ("S", "C", "N", "C", "C"): "1,3-thiazole",
    ("O", "C", "N", "C", "C"): "1,3-oxazole",
    ("S", "N", "C", "C", "C"): "isothiazole",
    ("O", "N", "C", "C", "C"): "isoxazole",
    ("N", "C", "N", "C", "C"): "1H-imidazole",
    ("N", "N", "C", "C", "C"): "1H-pyrazole",
}

# Phase 509: 2-ヘテロ原子6員部分不飽和環の親名
# (sig) → (parent_name, is_aromatic_parent, indicated_h_pos)
# is_aromatic_parent=True → 親が芳香族 (pyrimidine 等), dihydro ロカントをH数差分で計算
# is_aromatic_parent=False → 非芳香族親 (4H-1,3-oxazine 等), indicated_h_pos は親の固定sp3位
_DIHETERO_PARENT_6MEM: dict[tuple[str, ...], tuple[str, bool, int | None]] = {
    # N+N → 芳香族親
    ("N", "C", "N", "C", "C", "C"): ("pyrimidine", True, None),
    ("N", "N", "C", "C", "C", "C"): ("pyridazine", True, None),
    ("N", "C", "C", "N", "C", "C"): ("pyrazine",   True, None),
    # N+O → 非芳香族 4H-1,3-oxazine (O が最高優先 → pos1, N が pos3)
    ("O", "C", "N", "C", "C", "C"): ("4H-1,3-oxazine",  False, 4),
    # N+S → 非芳香族 4H-1,3-thiazine (S が最高優先 → pos1, N が pos3)
    ("S", "C", "N", "C", "C", "C"): ("4H-1,3-thiazine", False, 4),
}


def _match_partial_unsat(ring: list[int], graph: "MoleculeGraph") -> str | None:
    """
    部分不飽和単一ヘテロ原子環の dihydro 名を返す。
    例: 2,3-dihydrofuran, 1,2-dihydropyridine
    """
    from .molecule_analyzer import get_atom, get_bond_order

    if _is_aromatic_ring(ring, graph):
        return None
    n = len(ring)
    if n not in (5, 6):
        return None

    # 環内二重結合を収集
    ring_set = set(ring)
    seen: set[frozenset] = set()
    db_pairs: list[tuple[int, int]] = []
    for idx in ring:
        for nb in graph.adjacency[idx]:
            if nb in ring_set:
                key = frozenset((idx, nb))
                if key not in seen:
                    seen.add(key)
                    if get_bond_order(graph, idx, nb) == 2.0:
                        db_pairs.append((idx, nb))

    if not db_pairs:
        return None  # 完全飽和 → HW 名で処理

    # 単一ヘテロ原子のみ対応
    hetero_atoms = [idx for idx in ring if _atom_sig(graph, idx) != "C"]
    if len(hetero_atoms) != 1:
        return None
    hetero_idx = hetero_atoms[0]
    hetero_sym = get_atom(graph, hetero_idx).symbol
    has_nh = any(get_atom(graph, nb).symbol == "H" for nb in graph.adjacency[hetero_idx])

    parent = _PARTIAL_UNSAT_PARENT.get((n, hetero_sym))
    if parent is None:
        return None

    # sp3 原子 = 環内二重結合に関与しない炭素 (ヘテロ原子除く)
    db_atoms = {a for pair in db_pairs for a in pair}
    sp3_atoms = [idx for idx in ring if idx not in db_atoms and idx != hetero_idx]

    # 5員環N (1H-pyrrole親): N-H は親の indicated hydrogen → dihydro ロカントに含めない
    # 6員環N (pyridine親): N が sp3 になる場合は dihydro 位置として扱う
    pyrrole_parent = (n == 5 and hetero_sym == "N")
    sp3_with_hetero = sp3_atoms[:]
    if has_nh and hetero_idx not in db_atoms and not pyrrole_parent:
        sp3_with_hetero.append(hetero_idx)

    if not sp3_with_hetero:
        return None

    # ヘテロ原子を 1 位として両方向で sp3 ロカントを求め最小を選ぶ
    hi = ring.index(hetero_idx)
    best_locs: list[int] | None = None
    for direction in (1, -1):
        ordered = [ring[(hi + direction * k) % n] for k in range(n)]
        loc_map = {ordered[i]: i + 1 for i in range(n)}
        locs = sorted(loc_map[a] for a in sp3_with_hetero)
        if best_locs is None or locs < best_locs:
            best_locs = locs

    if best_locs is None:
        return None

    # N-H を持つ pyrrole 親 (1H-pyrrole): indicated hydrogen "-1H-" を接頭辞に付ける
    # pyridine 親で N が sp3: N は pos 1 に来るため 1 in best_locs → indicated_h 不要
    indicated_h = ""
    if hetero_sym == "N" and has_nh and 1 not in best_locs:
        indicated_h = "-1H-"

    mult = {2: "di", 3: "tri", 4: "tetra", 5: "penta"}.get(len(best_locs), "")
    locs_str = ",".join(str(l) for l in best_locs)
    return f"{locs_str}-{mult}hydro{indicated_h}{parent}"


def _ring_has_endocyclic_db(ring: list[int], graph: "MoleculeGraph") -> bool:
    """環内に非芳香族二重結合があれば True を返す。"""
    from .molecule_analyzer import get_bond_order
    ring_set = set(ring)
    seen: set[frozenset] = set()
    for idx in ring:
        for nb in graph.adjacency[idx]:
            if nb in ring_set:
                key = frozenset((idx, nb))
                if key not in seen:
                    seen.add(key)
                    if get_bond_order(graph, idx, nb) == 2.0:
                        return True
    return False


def _match_partial_unsat_2het(ring: list[int], graph: "MoleculeGraph") -> str | None:
    """
    Phase 508/509: 2ヘテロ原子5/6員部分不飽和環の dihydro 名を返す。
    例: C1=NCCS1  → '4,5-dihydro-1,3-thiazole'  (5員)
        C1=NCCCN1 → '1,4,5,6-tetrahydropyrimidine' (6員)
        C1=NCCCO1 → '5,6-dihydro-4H-1,3-oxazine'  (6員)
    """
    from .molecule_analyzer import get_bond_order, get_atom as _ga2

    n = len(ring)
    if n not in (5, 6):
        return None
    if _is_aromatic_ring(ring, graph):
        return None

    ring_set = set(ring)
    seen: set[frozenset] = set()
    db_pairs: list[tuple[int, int]] = []
    for idx in ring:
        for nb in graph.adjacency[idx]:
            if nb in ring_set:
                key = frozenset((idx, nb))
                if key not in seen:
                    seen.add(key)
                    if get_bond_order(graph, idx, nb) == 2.0:
                        db_pairs.append((idx, nb))

    if not db_pairs:
        return None

    hetero_atoms = [idx for idx in ring if _atom_sig(graph, idx) != "C"]
    if len(hetero_atoms) != 2:
        return None

    best_prio = min(_hetero_priority(_atom_sig(graph, idx)) for idx in hetero_atoms)
    primary_hets = [idx for idx in hetero_atoms
                    if _hetero_priority(_atom_sig(graph, idx)) == best_prio]

    best_key: tuple | None = None
    best_name: str | None = None

    if n == 5:
        db_atoms = {a for pair in db_pairs for a in pair}
        hetero_set = set(hetero_atoms)
        sp3_atoms = [idx for idx in ring if idx not in db_atoms and idx not in hetero_set]
        if not sp3_atoms:
            return None

        for start_idx in primary_hets:
            si = ring.index(start_idx)
            fwd_rot = ring[si:] + ring[:si]
            rev_rot = [ring[(si - k) % n] for k in range(n)]

            for rot in (fwd_rot, rev_rot):
                sig = tuple(_atom_sig(graph, x) for x in rot)
                parent = _DIHETERO_PARENT_BY_LOWLOC_SIG.get(sig)
                if parent is None:
                    continue

                loc_map = {rot[i]: i + 1 for i in range(n)}
                locs = sorted(loc_map[a] for a in sp3_atoms)
                hlocs = sorted(k + 1 for k, x in enumerate(rot) if _atom_sig(graph, x) != "C")
                start_has_h = int(not any(
                    _ga2(graph, nb).symbol == "H"
                    for nb in graph.adjacency[rot[0]]
                ))

                key = (hlocs, start_has_h, locs)
                if best_key is None or key < best_key:
                    best_key = key
                    mult = {1: "", 2: "di", 3: "tri", 4: "tetra"}.get(len(locs), "")
                    locs_str = ",".join(str(l) for l in locs)
                    sep = "-" if parent[0].isdigit() else ""
                    best_name = f"{locs_str}-{mult}hydro{sep}{parent}"

    else:  # n == 6
        db_set = set(frozenset(p) for p in db_pairs)

        for start_idx in primary_hets:
            si = ring.index(start_idx)
            fwd_rot = ring[si:] + ring[:si]
            rev_rot = [ring[(si - k) % n] for k in range(n)]

            for rot in (fwd_rot, rev_rot):
                sig = tuple(_atom_sig(graph, x) for x in rot)
                entry = _DIHETERO_PARENT_6MEM.get(sig)
                if entry is None:
                    continue

                parent_name, is_aromatic, indicated_h_pos = entry
                loc_map = {rot[i]: i + 1 for i in range(n)}
                hlocs = sorted(k + 1 for k, x in enumerate(rot) if _atom_sig(graph, x) != "C")

                if is_aromatic:
                    # 親が芳香族: 各位置のH数が芳香族親より多ければ dihydro ロカント
                    # 芳香族親: C→1H, N→0H
                    dihydro_locs: list[int] = []
                    for k, atom_idx in enumerate(rot):
                        h_count = sum(1 for nb in graph.adjacency[atom_idx]
                                      if _ga2(graph, nb).symbol == "H")
                        parent_h = 0 if _atom_sig(graph, atom_idx) != "C" else 1
                        if h_count > parent_h:
                            dihydro_locs.append(k + 1)
                else:
                    # 非芳香族親: db が pos2-3 にあることを確認
                    if frozenset((rot[1], rot[2])) not in db_set:
                        continue
                    db_atoms2 = {a for pair in db_pairs for a in pair}
                    sp3_c_locs = sorted(
                        loc_map[idx] for idx in ring
                        if idx not in db_atoms2 and _atom_sig(graph, idx) == "C"
                    )
                    if indicated_h_pos is not None:
                        dihydro_locs = [l for l in sp3_c_locs if l != indicated_h_pos]
                    else:
                        dihydro_locs = sp3_c_locs

                if not dihydro_locs:
                    continue

                start_h = sum(1 for nb in graph.adjacency[rot[0]]
                              if _ga2(graph, nb).symbol == "H")
                start_has_h = int(not (start_h > 0))

                key = (hlocs, start_has_h, sorted(dihydro_locs))
                if best_key is None or key < best_key:
                    best_key = key
                    locs = sorted(dihydro_locs)
                    mult_map = {1: "", 2: "di", 3: "tri", 4: "tetra", 5: "penta", 6: "hexa"}
                    mult = mult_map.get(len(locs), "")
                    locs_str = ",".join(str(l) for l in locs)
                    sep = "-" if parent_name[0].isdigit() else ""
                    best_name = f"{locs_str}-{mult}hydro{sep}{parent_name}"

    return best_name


# ─── 保留名テーブル ────────────────────────────────────────────────────
# key: (is_aromatic, canonical_signature_tuple)
# value: (retained_name, nh_indicator)
#   nh_indicator: True → "1H-" プレフィクスを付ける
#   注: canonical_sig は始端=最高優先ヘテロ原子、辞書順最小方向

_RETAINED_NAMES: dict[tuple[bool, tuple[str, ...]], tuple[str, bool]] = {
    # 6員芳香族ヘテロ環
    # pyridine: N-C-C-C-C-C (canonical min of fwd/rev)
    (True,  ("N", "C", "C", "C", "C", "C")): ("pyridine",  False),
    # pyridinone tautomer: NH in aromatic ring → pyridin-2(1H)-one etc.
    (True,  ("NH", "C", "C", "C", "C", "C")): ("pyridine",  True),
    # pyrazine: N at 1,4 → min("N","C","C","N","C","C") vs ("N","C","C","N","C","C") same
    (True,  ("N", "C", "C", "N", "C", "C")): ("pyrazine",  False),
    # pyrimidine: N at 1,3 → best_start=[N,C,N,C,C,C], rev=[N,C,C,C,N,C], min=("N","C","C","C","N","C")
    (True,  ("N", "C", "C", "C", "N", "C")): ("pyrimidine", False),
    # pyrimidinone tautomer: NH in aromatic ring → pyrimidin-{loc}({nh}H)-one (Phase 397)
    (True,  ("NH", "C", "C", "C", "N", "C")): ("pyrimidine", True),
    # pyrazinone tautomer: N at 1,4 → pyrazin-{loc}({nh}H)-one (Phase 398)
    (True,  ("NH", "C", "C", "N", "C", "C")): ("pyrazine",   True),
    # pyridazinone tautomer: N at 1,2 → pyridazin-{loc}({nh}H)-one (Phase 398)
    (True,  ("NH", "C", "C", "C", "C", "N")): ("pyridazine", True),
    # pyridazine: N at 1,2 → best_start=[N,N,C,C,C,C], rev=[N,C,C,C,C,N], min=("N","C","C","C","C","N")
    (True,  ("N", "C", "C", "C", "C", "N")): ("pyridazine", False),
    # 5員芳香族ヘテロ環
    (True,  ("O", "C", "C", "C", "C")): ("furan",      False),
    (True,  ("S", "C", "C", "C", "C")): ("thiophene",  False),
    # pyrrole: 1 NH → ("NH","C","C","C","C")
    (True,  ("NH", "C", "C", "C", "C")): ("pyrrole",   True),   # 1H-pyrrole
    # N-substituted pyrrole: N has no H (replaced by alkyl/acyl)
    (True,  ("N",  "C", "C", "C", "C")): ("pyrrole",   False),  # 1-substituted pyrrole
    # imidazole (1H): NH at 1, N at 3 → best_start=[NH,C,N,C,C], rev=[NH,C,C,N,C], min=("NH","C","C","N","C")
    (True,  ("NH", "C", "C", "N", "C")): ("imidazole", True),   # 1H-imidazole
    # N-substituted imidazole
    (True,  ("N",  "C", "C", "N", "C")): ("imidazole", False),  # 1-substituted imidazole
    # 飽和6員ヘテロ環
    (False, ("N", "C", "C", "C", "C", "C")): ("piperidine", False),
    (False, ("O", "C", "C", "N", "C", "C")): ("morpholine",  False),
    (False, ("N", "C", "C", "O", "C", "C")): ("morpholine",  False),
    # ピペラジン: N at 1,4 (飽和)
    (False, ("N", "C", "C", "N", "C", "C")): ("piperazine",  False),
    # NH ピペラジン (N-H あり)
    (False, ("NH", "C", "C", "N",  "C", "C")): ("piperazine",  False),
    (False, ("NH", "C", "C", "NH", "C", "C")): ("piperazine",  False),
    # Phase 154: 5員 O,N 飽和環 (oxazolidine = 1,3-oxazolidine; isoxazolidine = 1,2-oxazolidine)
    (False, ("O", "C", "C", "N", "C")):           ("oxazolidine",   False),
    (False, ("O", "C", "C", "C", "N")):           ("isoxazolidine", False),
    # Phase 154: 5員 S,N 飽和環 (thiazolidine = 1,3-; isothiazolidine = 1,2-)
    (False, ("S", "C", "C", "N", "C")):           ("thiazolidine",   False),
    (False, ("S", "C", "C", "C", "N")):           ("isothiazolidine", False),
    # Phase 154: 6員 O,N 飽和環 (1,3-oxazinane, 1,4-oxazinane 等)
    (False, ("O", "C", "C", "C", "N", "C")):      ("1,3-oxazinane",   False),
    (False, ("O", "C", "N", "C", "C", "C")):      ("1,3-oxazinane",   False),  # alt sig
    # Phase 515: 6員 1,2-O,N / 1,2-S,N 飽和環
    (False, ("O", "C", "C", "C", "C", "N")):      ("1,2-oxazinane",   False),
    (False, ("S", "C", "C", "C", "C", "N")):      ("1,2-thiazinane",  False),
    # Phase 154: 6員 S,N 飽和環
    (False, ("S", "C", "C", "N", "C", "C")):      ("1,4-thiazinane",  False),
    (False, ("S", "C", "C", "C", "N", "C")):      ("1,3-thiazinane",  False),
    # Phase 152: imidazolidine / pyrazolidine (飽和5員 N,N 環)
    (False, ("N", "C", "C", "N", "C")):           ("imidazolidine",   False),
    (False, ("N", "C", "C", "C", "N")):           ("pyrazolidine",    False),
    # Phase 153: テトラゾール (5員 aromatic 4-N 環)
    (True,  ("NH", "C", "N", "N", "N")):          ("tetrazole",       True),
    (True,  ("N",  "C", "N", "N", "N")):          ("tetrazole",       False),  # N-substituted
    # Phase 153: 1,3,5-triazinane (6員 tri-N 環; kekulization 失敗で aromatic 扱い)
    (True,  ("NH", "C", "NH", "C", "NH", "C")):   ("1,3,5-triazinane", False),
    (False, ("N",  "C", "N",  "C", "N",  "C")):   ("1,3,5-triazinane", False),
    # Phase 153: hexahydropyrimidine 型 (非芳香族6員 N1,N3 環)
    (False, ("N", "C", "C", "C", "N", "C")):      ("hexahydropyrimidine", False),
    # Phase 269: 5員二ヘテロ芳香族
    (True,  ("NH", "C", "C", "C", "N")): ("pyrazole",    True),   # 1H-pyrazole
    (True,  ("N",  "C", "C", "C", "N")): ("pyrazole",    False),  # 1-substituted
    (True,  ("O",  "C", "C", "C", "N")): ("isoxazole",   False),
    (True,  ("O",  "C", "C", "N", "C")): ("oxazole",     False),
    (True,  ("S",  "C", "C", "N", "C")): ("thiazole",    False),
    (True,  ("S",  "C", "C", "C", "N")): ("isothiazole", False),
    # Phase 401: pyrimidine di-NH tautomer → pyrimidine-2,4(1H,3H)-dione (uracil/thymine)
    (True,  ("NH", "C", "C", "C", "NH", "C")): ("pyrimidine", True),
    # Phase 400: 5員 oxazolone/thiazolone tautomers (NH-first sig, O/S in ring)
    (True,  ("NH", "C", "C", "O", "C")): ("1,3-oxazole",  True),
    (True,  ("NH", "C", "C", "S", "C")): ("1,3-thiazole", True),
    # Phase 269: 6員 triazine/tetrazine
    (True,  ("N",  "C", "C", "C", "N", "N")): ("1,2,3-triazine",    False),
    (True,  ("N",  "C", "C", "N", "C", "N")): ("1,2,4-triazine",    False),
    (True,  ("N",  "C", "N", "C", "N", "C")): ("1,3,5-triazine",    False),
    (True,  ("N",  "C", "N", "N", "C", "N")): ("1,2,4,5-tetrazine", False),
    # Phase 151: 二ヘテロ原子飽和環 (1,3-dioxolane, 1,4-dioxane 等)
    (False, ("O", "C", "C", "O", "C")):           ("1,3-dioxolane",   False),
    (False, ("O", "C", "C", "O", "C", "C")):      ("1,4-dioxane",     False),
    (False, ("O", "C", "C", "C", "O", "C")):      ("1,3-dioxane",     False),
    (False, ("S", "C", "C", "S", "C")):           ("1,3-dithiolane",  False),
    (False, ("S", "C", "C", "C", "S")):           ("1,2-dithiolane",  False),
    (False, ("S", "C", "C", "S", "C", "C")):      ("1,4-dithiane",    False),
    (False, ("S", "C", "C", "C", "C", "S")):      ("1,2-dithiane",    False),
    (False, ("O", "C", "C", "S", "C")):           ("1,3-oxathiolane", False),
    (False, ("O", "C", "C", "S", "C", "C")):      ("1,4-oxathiane",   False),
    # Phase 380: 環状スルファート/スルフィット (O,O,S 環)
    (False, ("O", "C", "C", "O", "S")):           ("1,3,2-dioxathiolane", False),
    (False, ("O", "C", "C", "C", "O", "S")):      ("1,3,2-dioxathiane",   False),
    # Phase 232: 1,3,5-trioxane (6員 tri-O 環)
    (False, ("O", "C", "O", "C", "O", "C")): ("1,3,5-trioxane", False),
    # Phase 548: 1,3,5-trithiane (6員 tri-S 環)
    (False, ("S", "C", "S", "C", "S", "C")): ("1,3,5-trithiane", False),
    # Phase 196: 4員環 二ヘテロ原子 (1,3-dioxetane 等)
    (False, ("O", "C", "O", "C")):                ("1,3-dioxetane",   False),
    (False, ("O", "C", "N", "C")):                ("1,3-oxazetidine", False),
    (False, ("N", "C", "O", "C")):                ("1,3-oxazetidine", False),
    (False, ("O", "C", "C", "N")):                ("1,2-oxazetidine", False),
    (False, ("S", "C", "S", "C")):                ("1,3-dithietane",  False),
    (False, ("O", "C", "S", "C")):                ("1,3-oxathietane", False),
}


def _is_aromatic_ring(ring: list[int], graph: "MoleculeGraph") -> bool:
    from .molecule_analyzer import get_atom
    return all(get_atom(graph, idx).is_aromatic for idx in ring)


def _match_retained(ring: list[int], graph: "MoleculeGraph") -> tuple[str, bool, list[int]] | None:
    """
    環が保留名に一致するか確認する。

    Returns:
        (name, is_nh, best_rotation) または None
        - is_nh: True → "1H-" プレフィクス
        - best_rotation: ロカント1が始端のリスト
    """
    is_arom = _is_aromatic_ring(ring, graph)
    best_rotation = _find_best_start(ring, graph)
    sig = _canonical_sig(best_rotation, graph)

    entry = _RETAINED_NAMES.get((is_arom, sig))
    if entry is not None:
        # Phase 508: 飽和保留名テーブルに一致しても環内二重結合があれば別系
        if not is_arom and _ring_has_endocyclic_db(ring, graph):
            return None
        return entry[0], entry[1], best_rotation

    return None


# ─── ロカントマップ構築 ────────────────────────────────────────────────

def _build_locant_map(rotation: list[int]) -> dict[int, int]:
    """始端=1 のロカントマップを返す。"""
    return {idx: i + 1 for i, idx in enumerate(rotation)}


# ─── 置換基収集（ヘテロ環用） ──────────────────────────────────────────

def _collect_hetero_substituents(
    graph: "MoleculeGraph",
    ring_atoms: list[int],
    locant_map: dict[int, int],
    excluded_atoms: "set[int] | None" = None,
) -> list[tuple[int | str, str]]:
    """
    ヘテロ環の置換基を収集する。
    excluded_atoms: Phase 22 外環官能基などを除外するためのインデックス集合。
    Non-aromatic N atoms get string locants ('N', "N'", ...) per IUPAC 2013.
    Returns: [(locant, substituent_name), ...]  ソート済み
    """
    from .molecule_analyzer import get_atom, get_bond_order
    from .substituent import name_substituent

    ring_set = set(ring_atoms)
    # Phase 404: pre-compute rings that contain at least one ring atom (for fused-ring check)
    _other_ring_sets: list[frozenset[int]] = [
        frozenset(r) for r in graph.ring_atom_sets if frozenset(r) != ring_set
    ]

    # Build N-locant map: non-aromatic N atoms sorted by position → "N", "N'", "N''", ...
    _N_LOCANTS = ["N", "N'", "N''", "N'''"]
    non_arom_ns = sorted(
        [idx for idx in ring_atoms
         if get_atom(graph, idx).symbol == "N" and not get_atom(graph, idx).is_aromatic],
        key=lambda idx: locant_map[idx],
    )
    n_locant_map: dict[int, str] = {
        idx: _N_LOCANTS[i] for i, idx in enumerate(non_arom_ns) if i < len(_N_LOCANTS)
    }

    result: list[tuple[int | str, str]] = []

    for ring_idx in ring_atoms:
        ring_atom = get_atom(graph, ring_idx)
        locant: int | str = n_locant_map.get(ring_idx, locant_map[ring_idx])

        for nb_idx in graph.adjacency[ring_idx]:
            if nb_idx in ring_set:
                continue
            if excluded_atoms and nb_idx in excluded_atoms:
                continue
            nb = get_atom(graph, nb_idx)
            if nb.symbol == "H":
                continue
            if ring_atom.symbol in ("N", "O", "S") and nb.symbol == "H":
                continue
            # Phase 404: skip fused ring junction bonds (nb is in another ring with ring_idx)
            if any(ring_idx in rs and nb_idx in rs for rs in _other_ring_sets):
                continue
            sub_name = name_substituent(graph, nb_idx, ring_set)
            result.append((locant, sub_name))

    # Sort: numeric locants first (ascending), then string locants (N < N' < N'')
    def _sort_key(t: tuple[int | str, str]) -> tuple:
        loc, nm = t
        if isinstance(loc, int):
            return (0, loc, nm)
        return (1, loc, nm)

    result.sort(key=_sort_key)
    return result


# ─── 名前組み立て ─────────────────────────────────────────────────────

def _format_substituents(
    base: str,
    substituents: list[tuple[int | str, str]],
) -> str:
    """置換基プレフィクスを付与した名前を組み立てる。"""
    if not substituents:
        return base

    import re
    from collections import defaultdict
    from .constants import MULTIPLIER

    # 同一置換基をまとめる: name -> [locants]  (locant may be int or str like "N")
    grouped: dict[str, list[int | str]] = defaultdict(list)
    for locant, name in substituents:
        grouped[name].append(locant)

    def _cpd_parens(nm: str) -> bool:
        if re.search(r"[0-9]", nm):
            return True
        for pfx in ("di", "tri", "tetra", "penta", "fluoro", "chloro", "bromo", "iodo"):
            if nm.startswith(pfx) and len(nm) > len(pfx):
                return True
        return False

    def _sort_locs(locs: list[int | str]) -> list[int | str]:
        return sorted(locs, key=lambda x: (isinstance(x, str), x))

    parts: list[str] = []
    for name in sorted(grouped):
        locs = _sort_locs(grouped[name])
        parens = _cpd_parens(name)
        n_locs = len(locs)
        if n_locs == 1:
            parts.append(f"{locs[0]}-({name})" if parens else f"{locs[0]}-{name}")
        else:
            loc_str = ",".join(str(l) for l in locs)
            mult = MULTIPLIER.get(n_locs, f"{n_locs}")
            parts.append(f"{loc_str}-{mult}({name})" if parens else f"{loc_str}-{mult}{name}")

    # アルファベット順: strip leading locant+hyphen or N'/N-prefix; numeric before N-locant
    def alpha_key(s: str) -> tuple:
        is_n = bool(re.match(r"^N[',]*-", s))
        stripped = re.sub(r"^[\d,\-]+", "", s)
        stripped = re.sub(r"^N[',]*-", "", stripped)
        stripped = re.sub(r"^(di|tri|tetra|penta|hexa|hepta|octa|nona|deca|bis|tris)", "", stripped)
        return (stripped.lower(), is_n)

    sorted_parts = sorted(parts, key=alpha_key)
    sep = "-" if base[:1].isdigit() else ""
    if len(sorted_parts) == 1:
        return f"{sorted_parts[0]}{sep}{base}"
    prefix = "-".join(sorted_parts[:-1]) + "-" + sorted_parts[-1]
    return f"{prefix}{sep}{base}"


def _apply_hetero_suffixes(
    full_base: str,
    substituents: list[tuple[int, str]],
) -> str:
    """
    ヘテロ芳香族の principal group を -amine/-ol/-thiol サフィックス形式に変換して返す。
    いずれも含まない場合は _format_substituents の結果を返す。
    ただし、保留名自体が PCG を内包している場合（coumarin, chromone, 等）は
    サフィックス変換をスキップしてプレフィックス形式を用いる。
    """
    if not substituents:
        return full_base

    from .constants import MULTIPLIER

    # 保留名が PCG (ケトン・ラクトン等) を内包する場合、-ol/-amine/-thiol への変換をスキップ
    _base_has_pcg = (
        full_base.endswith(("one", "dione"))
        or full_base in {"coumarin", "isocoumarin"}
    )
    if not _base_has_pcg:
        for sub_nm, suffix, elide_e in (
            ("amino", "amine", True),
            ("hydroxy", "ol", True),
            ("sulfanyl", "thiol", False),
        ):
            entries = [(loc, nm) for loc, nm in substituents if nm == sub_nm]
            if not entries:
                continue
            other = [(loc, nm) for loc, nm in substituents if nm != sub_nm]
            locs = sorted(loc for loc, _ in entries)
            loc_str = ",".join(str(l) for l in locs)
            mult = MULTIPLIER.get(len(locs), "") if len(locs) > 1 else ""
            stem = full_base[:-1] if (elide_e and full_base.endswith("e")) else full_base
            base_with_suffix = f"{stem}-{loc_str}-{mult}{suffix}"
            if not other:
                return base_with_suffix
            return _format_substituents(base_with_suffix, other)

    return _format_substituents(full_base, substituents)


# ─── 縮合ヘテロ芳香族保留名テーブル (Phase 17) ───────────────────────────
# key: RDKit canonical SMILES (AddHs 前)
# value: IUPAC 保留名
_FUSED_HETERO_RETAINED: dict[str, str] = {
    "c1ccc2ncccc2c1":   "quinoline",
    "c1ccc2cnccc2c1":   "isoquinoline",
    "c1ccc2[nH]ccc2c1": "1H-indole",
    "c1ccc2[nH]cnc2c1": "1H-benzimidazole",
    "c1ccc2nccnc2c1":   "quinoxaline",        # benzo[g]pyrazine (Phase 130)
    "c1ccc2occc2c1":    "benzofuran",
    "c1ccc2sccc2c1":    "benzo[b]thiophene",
    "c1ncc2[nH]cnc2n1": "7H-purine",
    "c1ccc2nc3ccccc3cc2c1": "acridine",
    # Phase 130: 追加縮合ヘテロ芳香族
    "c1ccc2ocnc2c1":    "1,3-benzoxazole",
    "c1ccc2scnc2c1":    "1,3-benzothiazole",
    "c1ccc2ncncc2c1":   "quinazoline",
    "c1ccc2nnccc2c1":   "cinnoline",
    "c1ccc2cnncc2c1":   "phthalazine",
    "c1ccc2nc3ccccc3nc2c1": "phenazine",
    "c1cnc2ncncc2n1":   "pteridine",
    # Phase 131: 三環式縮合ヘテロ芳香族
    "c1ccc2c(c1)[nH]c1ccccc12": "9H-carbazole",
    "c1ccc2c(c1)sc1ccccc12":    "dibenzothiophene",
    "c1ccc2c(c1)oc1ccccc12":    "dibenzofuran",
    # Phase 132: ナフチリジン全6異性体 (canonical SMILES = RDKit canonical)
    "c1cnc2cccnc2c1":   "1,5-naphthyridine",
    "c1cnc2ccncc2c1":   "1,6-naphthyridine",
    "c1cnc2cnccc2c1":   "1,7-naphthyridine",
    "c1cnc2ncccc2c1":   "1,8-naphthyridine",
    "c1cc2cnccc2cn1":   "2,6-naphthyridine",
    "c1cc2ccncc2cn1":   "2,7-naphthyridine",
    # Phase 132: 縮合ヘテロ芳香族追加 (IUPAC 2013 P-31.1.3)
    "c1ccn2cccc2c1":    "indolizine",
    "c1ccc2[nH]ncc2c1": "1H-indazole",
    "c1ccc2n[nH]cc2c1": "2H-indazole",
    "c1ccc2[nH]nnc2c1": "1H-benzotriazole",
    "c1ccc2oncc2c1":    "1,2-benzisoxazole",
    # Phase 133: 部分飽和縮合環 保留名 (IUPAC 2013 P-31.1.2, P-31.1.6)
    "c1ccc2c(c1)CCC2":  "indane",
    "c1ccc2c(c1)CCCC2": "1,2,3,4-tetrahydronaphthalene",
    "C1=Cc2ccccc2CC1":  "1,2-dihydronaphthalene",
    # Phase 633: 4,5,6,7-tetrahydrobenzo-fused 5-membered aromatic heterocycles
    "c1cc2c(o1)CCCC2":    "4,5,6,7-tetrahydrobenzofuran",
    "c1cc2c(s1)CCCC2":    "4,5,6,7-tetrahydrobenzothiophene",
    "c1cc2c([nH]1)CCCC2": "4,5,6,7-tetrahydro-1H-indole",
    "c1nc2c([nH]1)CCCC2": "4,5,6,7-tetrahydro-1H-benzimidazole",
    "c1n[nH]c2c1CCCC2":   "4,5,6,7-tetrahydro-1H-indazole",
    # Phase 634: 4,5,6,7-tetrahydrobenzo[d] 5-membered heterocycles (oxazole/thiazole/isoxazole/isothiazole/triazole)
    "c1nc2c(o1)CCCC2":   "4,5,6,7-tetrahydrobenzo[d]oxazole",
    "c1nc2c(s1)CCCC2":   "4,5,6,7-tetrahydrobenzo[d]thiazole",
    "c1noc2c1CCCC2":     "4,5,6,7-tetrahydrobenzo[d]isoxazole",
    "c1nsc2c1CCCC2":     "4,5,6,7-tetrahydrobenzo[d]isothiazole",
    "C1CCc2[nH]nnc2C1":  "4,5,6,7-tetrahydro-1H-benzo[d][1,2,3]triazole",
    "C1CCc2n[nH]nc2C1":  "4,5,6,7-tetrahydro-2H-benzo[d][1,2,3]triazole",
    # Phase 632: 5,6,7,8-tetrahydro N-heterocycles (aromatic ring intact, benzo ring saturated)
    "c1cnc2c(c1)CCCC2":  "5,6,7,8-tetrahydroquinoline",
    "c1cc2c(cn1)CCCC2":  "5,6,7,8-tetrahydroisoquinoline",
    "c1cnc2c(n1)CCCC2":  "5,6,7,8-tetrahydroquinoxaline",
    "c1ncc2c(n1)CCCC2":  "5,6,7,8-tetrahydroquinazoline",
    "c1cc2c(nn1)CCCC2":  "5,6,7,8-tetrahydrocinnoline",
    "c1nncc2c1CCCC2":    "5,6,7,8-tetrahydrophthalazine",
    "c1ccc2c(c1)CNCN2":  "1,2,3,4-tetrahydroquinazoline",
    "C1=CCN2C=CC=CC2=C1": "quinolizine",
    # Phase 632: additional polycyclic aromatic hydrocarbons
    "c1ccc2cc3c(ccc4ccccc43)cc2c1": "benz[a]anthracene",
    "c1ccc2c(c1)Cc1c-2ccc2ccccc12": "benzo[a]fluorene",
    "c1ccc2c(c1)Cc1cc3ccccc3cc1-2": "benzo[b]fluorene",
    "c1ccc2c(c1)Cc1ccc3ccccc3c1-2": "benzo[c]fluorene",
    "c1ccc2c(c1)ccc1cc3c(ccc4ccccc43)cc12": "dibenz[a,h]anthracene",
    # Phase 631: tricyclic partially saturated ring systems
    "c1ccc2c3c(ccc2c1)CCCC3":  "1,2,3,4-tetrahydrophenanthrene",
    "c1ccc2c(c1)CCc1ccccc1-2": "9,10-dihydrophenanthrene",
    "C1=Cc2c(ccc3ccccc23)CC1": "1,2-dihydrophenanthrene",
    "c1ccc2cc3c(cc2c1)CCCC3":  "1,2,3,4-tetrahydroanthracene",
    "c1ccc2c(c1)Cc1ccccc1C2":  "9,10-dihydroanthracene",
    "C1=Cc2cc3ccccc3cc2CC1":   "1,2-dihydroanthracene",
    "c1cc2c3c(cccc3c1)CCC2":   "2,3-dihydro-1H-phenalene",
    "c1ccc2c(c1)CCN2":  "indoline",
    "c1ccc2c(c1)CCO2":  "2,3-dihydrobenzofuran",
    "c1ccc2c(c1)COC2":  "1,3-dihydro-2-benzofuran",
    "c1ccc2c(c1)CCS2":  "2,3-dihydrobenzothiophene",
    "c1ccc2c(c1)CSC2":  "1,3-dihydro-2-benzothiophene",
    # Phase 635: 2,3-dihydrobenzo[d] heterocycles (5-membered ring with N/O/S, aromatic benzo ring)
    "c1ccc2c(c1)NCO2":  "2,3-dihydrobenzo[d]oxazole",
    "c1ccc2c(c1)NCS2":  "2,3-dihydrobenzo[d]thiazole",
    "c1ccc2c(c1)CNO2":  "2,3-dihydrobenzo[d]isoxazole",
    "c1ccc2c(c1)CNS2":  "2,3-dihydrobenzo[d]isothiazole",
    "c1ccc2c(c1)NCN2":  "2,3-dihydro-1H-benzimidazole",
    "c1ccc2c(c1)CNN2":  "2,3-dihydro-1H-indazole",
    # Phase 636: 5,6,7,8-tetrahydronaphthyridines and 1,2,3,4-tetrahydronaphthyridines
    "c1cnc2c(c1)NCCC2":  "1,2,3,4-tetrahydro-1,5-naphthyridine",
    "c1cnc2c(c1)CNCC2":  "5,6,7,8-tetrahydro-1,6-naphthyridine",
    "c1cnc2c(c1)CCNC2":  "5,6,7,8-tetrahydro-1,7-naphthyridine",
    "c1cnc2c(c1)CCCN2":  "1,2,3,4-tetrahydro-1,8-naphthyridine",
    "c1cc2c(cn1)CCCN2":  "1,2,3,4-tetrahydro-1,6-naphthyridine",
    "c1cc2c(cn1)NCCC2":  "1,2,3,4-tetrahydro-1,7-naphthyridine",
    # Phase 637: 1,2,3,4-tetrahydronaphthyridines (2,x series, symmetric naphthyridines)
    "c1cc2c(cn1)CCNC2":  "1,2,3,4-tetrahydro-2,6-naphthyridine",
    "c1cc2c(cn1)CNCC2":  "1,2,3,4-tetrahydro-2,7-naphthyridine",
    # Phase 638: tetrahydroacridine, tetrahydrophenanthridine, dihydroacridine, octahydroacridine
    "c1ccc2nc3c(cc2c1)CCCC3": "1,2,3,4-tetrahydroacridine",
    "c1ccc2c3c(ncc2c1)CCCC3": "1,2,3,4-tetrahydrophenanthridine",
    "c1ccc2c(c1)Cc1ccccc1N2": "9,10-dihydroacridine",
    "c1c2c(nc3c1CCCC3)CCCC2": "1,2,3,4,5,6,7,8-octahydroacridine",
    # Phase 639: 1,2,3,4-tetrahydrobenzo[g/f/h]quinoline
    "c1ccc2cc3c(cc2c1)CCCN3": "1,2,3,4-tetrahydrobenzo[g]quinoline",
    "c1ccc2c3c(ccc2c1)NCCC3": "1,2,3,4-tetrahydrobenzo[f]quinoline",
    "c1ccc2c3c(ccc2c1)CCCN3": "1,2,3,4-tetrahydrobenzo[h]quinoline",
    # Phase 640: 9,10-dihydrophenanthridine and 5,6-dihydrophenanthridine
    "C1=Cc2cnc3ccccc3c2CC1": "9,10-dihydrophenanthridine",
    "c1ccc2c(c1)CNc1ccccc1-2": "5,6-dihydrophenanthridine",
    # Phase 641: 1,2,3,4-tetrahydrobenzo[f/g/h]isoquinoline
    "c1ccc2c3c(ccc2c1)CNCC3": "1,2,3,4-tetrahydrobenzo[f]isoquinoline",
    "c1ccc2cc3c(cc2c1)CCNC3": "1,2,3,4-tetrahydrobenzo[g]isoquinoline",
    "c1ccc2c3c(ccc2c1)CCNC3": "1,2,3,4-tetrahydrobenzo[h]isoquinoline",
    # Phase 642: 1,2-dihydroacridine
    "C1=Cc2nc3ccccc3cc2CC1": "1,2-dihydroacridine",
    # Phase 643: 3,4-dihydroacridine and 1,2-dihydrophenanthridine
    "C1=Cc2cc3ccccc3nc2CC1": "3,4-dihydroacridine",
    "C1=Cc2ncc3ccccc3c2CC1": "1,2-dihydrophenanthridine",
    # Phase 644: 3,4-dihydrophenanthridine
    "C1=Cc2c(ncc3ccccc23)CC1": "3,4-dihydrophenanthridine",
    "c1ccc2c(c1)CCCO2": "chromane",
    "c1ccc2c(c1)CCOC2": "isochromane",
    "c1ccc2c(c1)CCCN2": "1,2,3,4-tetrahydroquinoline",
    "c1ccc2c(c1)CCNC2": "1,2,3,4-tetrahydroisoquinoline",
    "c1ccc2c(c1)NCCN2": "1,2,3,4-tetrahydroquinoxaline",
    "c1ccc2c(c1)CCCS2": "thiochroman",
    # Phase 630: benzo-fused 7-membered saturated rings (benzazepines, benzoxepines, etc.)
    "c1ccc2c(c1)CCCCC2": "6,7,8,9-tetrahydro-5H-benzo[7]annulene",
    "c1ccc2c(c1)CCCCN2": "2,3,4,5-tetrahydro-1-benzazepine",
    "c1ccc2c(c1)CCCNC2": "2,3,4,5-tetrahydro-1H-2-benzazepine",
    "c1ccc2c(c1)CCNCC2": "2,3,4,5-tetrahydro-1H-3-benzazepine",
    "c1ccc2c(c1)CCCCO2": "2,3,4,5-tetrahydro-1-benzoxepine",
    "c1ccc2c(c1)CCCOC2": "2,3,4,5-tetrahydro-1H-2-benzoxepine",
    "c1ccc2c(c1)CCOCC2": "2,3,4,5-tetrahydro-1H-3-benzoxepine",
    "c1ccc2c(c1)CCCCS2": "2,3,4,5-tetrahydro-1-benzothiepine",
    "c1ccc2c(c1)NCCCN2": "2,3,4,5-tetrahydro-1H-1,5-benzodiazepine",
    "C1=Nc2ccccc2NCC1":  "2,3-dihydro-1H-1,5-benzodiazepine",
    "O=C1CCCc2ccccc2N1": "3,4-dihydro-1H-1-benzazepin-2(5H)-one",
    "O=C1NCCNc2ccccc21": "2,3-dihydro-1H-1,4-benzodiazepin-5(4H)-one",
    # Phase 628: benzo-fused di-heteroatom 6-membered saturated rings
    "c1ccc2c(c1)OCCO2": "2,3-dihydro-1,4-benzodioxine",
    "c1ccc2c(c1)NCCO2": "3,4-dihydro-2H-1,4-benzoxazine",
    "c1ccc2c(c1)NCCS2": "3,4-dihydro-2H-1,4-benzothiazine",
    # Phase 628: benzo-fused partially unsaturated N-heterocycles (one C=N or C=C remains)
    "C1=NCCc2ccccc21":  "3,4-dihydroisoquinoline",
    "C1=Nc2ccccc2CC1":  "3,4-dihydroquinoline",
    "C1=Nc2ccccc2NC1":  "3,4-dihydroquinoxaline",
    "C1=Cc2ccccc2NC1":  "1,2-dihydroquinoline",
    # Phase 629: benzo-fused lactams with dihydro prefix
    "O=C1CNc2ccccc2N1": "3,4-dihydroquinoxalin-2(1H)-one",
    "O=C1CCc2ccccc21":  "indan-1-one",
    "O=C1CC(=O)c2ccccc21": "indane-1,3-dione",
    "O=C1Cc2ccccc2N1":  "indolin-2-one",
    "O=C1NCc2ccccc21":  "isoindolin-1-one",
    # Phase 409: 6-membered benzo-fused lactams
    "O=C1CCc2ccccc2N1":  "3,4-dihydroquinolin-2(1H)-one",
    "O=C1NCCc2ccccc21":  "3,4-dihydroisoquinolin-1(2H)-one",
    # Phase 264: 縮合無水物・マレイン酸無水物
    "O=C1OC(=O)c2ccccc21": "isobenzofuran-1,3-dione",
    "O=C1C=CC(=O)O1":      "furan-2,5-dione",
    # Phase 267: 1,3-ベンゾジオキソール (メチレンジオキシベンゼン骨格)
    "c1ccc2c(c1)OCO2":  "1,3-benzodioxole",
    # Phase 134: 追加縮合環保留名 (IUPAC 2013 P-31.1.3, fluorene/xanthene 系)
    "c1ccc2c(c1)Cc1ccccc1-2":  "fluorene",
    "c1ccc2c(c1)Cc1ccccc1O2":  "xanthene",
    "c1ccc2c(c1)Cc1ccccc1S2":  "thioxanthene",
    "O=c1ccc2ccccc2o1":        "coumarin",
    "O=c1occc2ccccc12":        "isocoumarin",
    "C1=COc2ccccc2C1":         "4H-chromene",
    "C1=Cc2ccccc2OC1":         "2H-chromene",
    "c1ccc2c(c1)Nc1ccccc1O2":  "phenoxazine",
    "c1ccc2c(c1)Nc1ccccc1S2":  "phenothiazine",
    "O=c1c2ccccc2oc2ccccc12":  "xanthen-9-one",
    "O=C1c2ccccc2-c2ccccc21":  "fluoren-9-one",
    "O=c1c2ccccc2sc2ccccc12":  "thioxanthen-9-one",
    # Phase 410: benzo-fused saturated ketones and lactones
    "O=C1CCCc2ccccc21":  "3,4-dihydronaphthalen-1(2H)-one",
    "O=C1CCc2ccccc2C1":  "3,4-dihydronaphthalen-2(1H)-one",
    "O=C1OCCc2ccccc21":  "isochroman-1-one",
    "O=C1CCc2ccccc2O1":  "chroman-2-one",
    "O=C1CCOc2ccccc21":  "chroman-4-one",
    # Phase 411: chromone and indan-2-one
    "O=c1ccoc2ccccc12":  "chromone",
    "O=C1Cc2ccccc2C1":   "indan-2-one",
    # Phase 412: benzo-fused azolone C2=O derivatives
    "O=c1[nH]c2ccccc2[nH]1":  "1H-benzimidazol-2(3H)-one",
    "O=c1[nH]c2ccccc2s1":     "1,3-benzothiazol-2(3H)-one",
    "O=c1[nH]c2ccccc2o1":     "1,3-benzoxazol-2(3H)-one",
    # Phase 413: acridinone and quinoxalinone
    "O=c1c2ccccc2[nH]c2ccccc12":  "acridin-9(10H)-one",
    "O=c1cnc2ccccc2[nH]1":        "quinoxalin-2(1H)-one",
    # Phase 414: benzo-fused diazinones
    "O=c1[nH]ncc2ccccc12":  "phthalazin-1(2H)-one",
    "O=c1cn[nH]c2ccccc12":  "cinnolin-4(1H)-one",
    "O=c1[nH]cnc2ccccc12":  "quinazolin-4(3H)-one",
    # Phase 415: anthrone, quinolinone, quinazolinedione
    "O=C1c2ccccc2Cc2ccccc21":      "anthracen-9(10H)-one",
    "O=c1ccc2ccccc2[nH]1":         "quinolin-2(1H)-one",
    "O=c1[nH]c(=O)c2ccccc2[nH]1": "quinazoline-2,4(1H,3H)-dione",
    # Phase 416: quinolinones and isoquinolinones
    "O=c1cc[nH]c2ccccc12":  "quinolin-4(1H)-one",
    "O=c1[nH]ccc2ccccc12":  "isoquinolin-1(2H)-one",
    "O=c1cc2ccccc2c[nH]1":  "isoquinolin-3(2H)-one",
    # Phase 417: indazolone, phenanthridinone, phenanthridine
    "O=c1[nH][nH]c2ccccc12":         "1H-indazol-3(2H)-one",
    "O=c1[nH]c2ccccc2c2ccccc12":     "phenanthridin-6(5H)-one",
    "c1ccc2c(c1)cnc1ccccc12":        "phenanthridine",
    # Phase 418: isatin, benzofuranone, benzothiophenone
    "O=C1Nc2ccccc2C1=O":  "1H-indole-2,3-dione",
    "O=C1COc2ccccc21":    "benzofuran-2(3H)-one",
    "O=C1CSc2ccccc21":    "benzo[b]thiophen-2(3H)-one",
    # Phase 419: naphthalenones and 1,3-benzodioxol-2-one
    "O=C1CC=Cc2ccccc21":  "naphthalen-1(2H)-one",
    "O=C1C=Cc2ccccc2C1":  "naphthalen-2(1H)-one",
    "O=c1oc2ccccc2o1":    "1,3-benzodioxol-2-one",
    # Phase 420: monocyclic pyranones (alpha- and gamma-pyrone)
    "O=c1cccco1":  "2H-pyran-2-one",
    "O=c1ccocc1":  "4H-pyran-4-one",
    # Phase 421: acenaphthylene and acenaphthene (tricyclic 5+6+6 PAH) — duplicate of Phase 138 keys
    "C1=Cc2cccc3cccc1c23":  "acenaphthylene",
    "c1cc2c3c(cccc3c1)CC2":  "acenaphthene",
    # Phase 422: fluoranthene (C16H10, correct); pyrene already in Phase 138
    "c1cc2ccc3cccc4ccc(c1)c2c34":    "pyrene",
    "c1ccc2c(c1)-c1cccc3cccc-2c13": "fluoranthene",
    # Phase 423: phenoxathiin (new); triphenylene already in Phase 138
    "c1ccc2c(c1)Oc1ccccc1S2":         "phenoxathiin",
    "c1ccc2c(c1)c1ccccc1c1ccccc21":   "triphenylene",
    # Phase 424: tetracene (naphthacene) — chrysene already in Phase 138
    "c1ccc2cc3cc4ccccc4cc3cc2c1":    "tetracene",
    # Phase 426: thianthrene (S,S-bridged dibenzene) and correct benzo[a]pyrene (C20H12)
    "c1ccc2c(c1)Sc1ccccc1S2":             "thianthrene",
    "c1ccc2c(c1)cc1ccc3cccc4ccc2c1c34":  "benzo[a]pyrene",
    # Phase 427: 1H-phenalene (3×6 peri-fused, C13H10, sp3 CH2 at position 1)
    "C1=Cc2cccc3cccc(c23)C1":  "1H-phenalene",
    # Phase 428: phenanthroline retained names (IUPAC 2013 P-31.1.3.4)
    "c1cnc2c(c1)ccc1cccnc12":  "1,10-phenanthroline",
    "c1cnc2ccc3ncccc3c2c1":    "4,7-phenanthroline",
    "c1cnc2c(c1)ccc1ncccc12":  "1,7-phenanthroline",
    # Phase 429: benzo[f]quinoline and benzo[h]quinoline (IUPAC 2013 P-31.1.3)
    "c1ccc2c(c1)ccc1ncccc12":  "benzo[f]quinoline",
    "c1ccc2c(c1)ccc1cccnc12":  "benzo[h]quinoline",
    # Phase 430: benzo[f]isoquinoline and benzo[g]isoquinoline (IUPAC 2013 P-31.1.3)
    "c1ccc2cc3cnccc3cc2c1":    "benzo[g]isoquinoline",
    "c1ccc2c(c1)ccc1cnccc12":  "benzo[f]isoquinoline",
    # Phase 431: benz[a]acridine and benz[c]acridine (IUPAC 2013 P-31.1.3)
    "c1ccc2nc3ccc4ccccc4c3cc2c1":  "benz[a]acridine",
    "c1ccc2cc3nc4ccccc4cc3cc2c1":  "benz[c]acridine",
    # Phase 432: pyrido[b/c]indole carbolines (IUPAC 2013 P-31.1.3)
    # Note: [4,3-b] and [3,4-c] use old SMILES (not yet OPSIN-corrected);
    # [3,4-b] and [2,3-b] use OPSIN canonical SMILES and must appear LAST to win.
    "c1ccc2c(c1)[nH]c1ncccc12":   "9H-pyrido[4,3-b]indole",
    "c1ccc2c(c1)[nH]c1cnccc12":   "9H-pyrido[3,4-c]indole",
    "c1ccc2c(c1)[nH]c1cnccc12":   "9H-pyrido[3,4-b]indole",
    "c1ccc2c(c1)[nH]c1ncccc12":   "9H-pyrido[2,3-b]indole",
    # Phase 433: 1H-benzo[e]indole and 1H-benzo[f]indole (IUPAC 2013 P-31.1.3)
    # Note: 1H-benzo[e]indole OPSIN gives non-aromatic SMILES; skip that entry.
    # c1ccc2cc3[nH]ccc3cc2c1 is the OPSIN canonical for 1H-benzo[f]indole.
    "c1ccc2cc3[nH]ccc3cc2c1":     "1H-benzo[f]indole",
    # Phase 584: 1H-benzo[g]indole (IUPAC 2013 P-31.1.3)
    "c1ccc2c(c1)ccc1cc[nH]c12":   "1H-benzo[g]indole",
    # Phase 434: naphtho[2,1-b]furan and naphtho[2,1-b]thiophene (IUPAC 2013 P-31.1.3)
    "c1ccc2c(c1)ccc1occc12":      "naphtho[2,1-b]furan",
    "c1ccc2c(c1)ccc1sccc12":      "naphtho[2,1-b]thiophene",
    # Phase 435: drug-scaffold bicyclic heterocycles (IUPAC 2013 P-31.1.3)
    # corrected Phase 473: c1cnc2 ring6 is pyrazine (N at 1,4), not pyrimidine
    "c1c[nH]c2ccnc-2n1":  "1H-pyrrolo[2,3-e]pyrazine",
    "c1cnc2[nH]ncc2n1":  "1H-pyrazolo[4,5-e]pyrazine",
    "c1cc2nc[nH]c2cn1":  "1H-imidazo[4,5-c]pyridine",
    "c1cnc2sccc2n1":     "thieno[2,3-e]pyrazine",
    # Phase 436: more drug-scaffold bicyclics
    # corrected Phase 473: c1cnc2 ring6 is pyrazine
    "c1cnc2occc2n1":     "furo[2,3-e]pyrazine",
    "c1cnc2cocc2n1":     "furo[3,4-e]pyrazine",
    "c1ccn2nccc2c1":     "pyrazolo[1,5-a]pyridine",
    "c1ccn2cnnc2c1":     "[1,2,4]triazolo[4,3-a]pyridine",
    "c1cnc2occc2c1":     "furo[2,3-b]pyridine",
    "c1cnc2[nH]ncc2c1":  "1H-pyrazolo[3,4-b]pyridine",
    # Phase 437: imidazo/thieno-pyrimidine and triazolo-pyrimidine
    # corrected Phase 473: c1cnc2 ring6 is pyrazine
    "c1cnc2nccn2c1":     "imidazo[1,2-a]pyrimidine",
    "c1cnc2[nH]cnc2n1":  "1H-imidazo[4,5-e]pyrazine",
    "c1cnc2cscc2n1":     "thieno[3,4-e]pyrazine",
    "c1cnc2cnnn2c1":     "[1,2,3]triazolo[1,5-a]pyrimidine",
    # Phase 473: remaining pyrazine e-bond fusions
    # pyrazine C2 symmetry collapses [2,3]/[3,2], [4,5]/[5,4], and imidazo[4,5]/[5,4] each to one compound
    "c1c[nH]c2cncc-2n1":  "1H-pyrrolo[3,4-e]pyrazine",
    # Phase 474: isoxazolo/pyrazolo/[1,2,3]triazolo[x,y-e]pyrazine (IUPAC 2013 P-31.1.3)
    # pyrazine C2: [x,y] and [y,x] same compound; use lower locants; [4,5]<[5,4], [3,4]<[4,3]
    "c1cnc2nocc2n1":     "isoxazolo[3,4-e]pyrazine",
    "c1cnc2oncc2n1":     "isoxazolo[4,5-e]pyrazine",
    "c1cnc2n[nH]cc2n1":  "1H-pyrazolo[3,4-e]pyrazine",
    # triazolo and pyrazolo[4,5-e]: canonical SMILES also in Phase 469 section (corrected below)
    # Phase 475: isothiazolo[x,y-d]pyrimidine and [x,y-e]pyrazine (IUPAC 2013 P-31.1.3)
    "c1ncc2sncc2n1":     "isothiazolo[4,5-d]pyrimidine",
    "c1ncc2nscc2n1":     "isothiazolo[4,3-d]pyrimidine",
    "c1ncc2cnsc2n1":     "isothiazolo[5,4-d]pyrimidine",
    "c1ncc2csnc2n1":     "isothiazolo[3,4-d]pyrimidine",
    "c1cnc2sncc2n1":     "isothiazolo[4,5-e]pyrazine",
    "c1cnc2nscc2n1":     "isothiazolo[3,4-e]pyrazine",
    # Phase 476: isothiazolo[x,y-b]pyridine (4 isomers, IUPAC 2013 P-31.1.3)
    "c1cnc2nscc2c1":     "isothiazolo[3,4-b]pyridine",
    "c1cnc2csnc2c1":     "isothiazolo[4,3-b]pyridine",
    "c1cnc2cnsc2c1":     "isothiazolo[4,5-b]pyridine",
    "c1cnc2sncc2c1":     "isothiazolo[5,4-b]pyridine",
    # Phase 477: isothiazolo[x,y-c]pyridine (4 isomers, IUPAC 2013 P-31.1.3)
    "c1cc2csnc2cn1":     "isothiazolo[3,4-c]pyridine",
    "c1cc2nscc2cn1":     "isothiazolo[4,3-c]pyridine",
    "c1cc2sncc2cn1":     "isothiazolo[4,5-c]pyridine",
    "c1cc2cnsc2cn1":     "isothiazolo[5,4-c]pyridine",
    # Phase 478: isothiazolo[x,y-c]pyridazine (4 isomers, IUPAC 2013 P-31.1.3)
    "c1cc2csnc2nn1":     "isothiazolo[3,4-c]pyridazine",
    "c1cc2nscc2nn1":     "isothiazolo[4,3-c]pyridazine",
    "c1cc2sncc2nn1":     "isothiazolo[4,5-c]pyridazine",
    "c1cc2cnsc2nn1":     "isothiazolo[5,4-c]pyridazine",
    # Phase 438: pyrido-pyrimidine, pyrimido-pyrimidine, thieno-pyridine
    "c1cnc2ncncc2c1":    "pyrido[2,3-d]pyrimidine",
    "c1cc2cncnc2cn1":    "pyrido[3,4-d]pyrimidine",
    "c1cnc2ncncc2n1":    "pteridine",
    "c1cnc2cscc2c1":     "thieno[3,4-b]pyridine",
    # Phase 439: oxazolo/thiazolo-pyridine and -pyrazine (fixed Phase 470)
    "c1cnc2ocnc2c1":     "oxazolo[5,4-b]pyridine",
    "c1cnc2ncoc2c1":     "oxazolo[4,5-b]pyridine",
    "c1cnc2scnc2c1":     "thiazolo[5,4-b]pyridine",
    "c1cnc2ncsc2c1":     "thiazolo[4,5-b]pyridine",
    "c1cnc2ocnc2n1":     "oxazolo[4,5-e]pyrazine",
    "c1cnc2scnc2n1":     "thiazolo[4,5-e]pyrazine",
    # Phase 440: more thieno-pyridine and triazolo-pyridine
    "c1cc2cscc2cn1":     "thieno[3,4-c]pyridine",
    "c1ccn2nncc2c1":     "[1,2,3]triazolo[1,5-a]pyridine",
    "c1cnc2[nH]nnc2c1":  "1H-[1,2,3]triazolo[4,5-b]pyridine",
    "C1=Nc2cccc3cccc(c23)N1":  "perimidine",
    # Phase 142: 追加ヘテロ芳香族 (セレノフェン、縮合二環式)
    "c1cc[se]c1":       "selenophene",
    # Phase 255: テルロフェン
    "c1cc[te]c1":       "tellurophene",
    "c1ccn2ccnc2c1":    "imidazo[1,2-a]pyridine",
    # Phase 480: [1,2,4]triazolo[1,5-a]pyridine (IUPAC 2013 P-31.1.3)
    "c1ccn2ncnc2c1":    "[1,2,4]triazolo[1,5-a]pyridine",
    "c1cnc2ccnn2c1":    "pyrazolo[1,5-a]pyrimidine",
    # Phase 481: imidazo[1,2-b]pyridazine (IUPAC 2013 P-31.1.3)
    "c1cnn2ccnc2c1":    "imidazo[1,2-b]pyridazine",
    # Phase 482: missing imidazo/pyrazolo/triazolo fused bicyclics (IUPAC 2013 P-31.1.3)
    "c1cc2nccn2cn1":    "imidazo[1,2-c]pyrimidine",
    "c1cn2ccnc2cn1":    "imidazo[1,2-a]pyrazine",
    "c1cnc2cncn2c1":    "imidazo[1,5-a]pyrimidine",
    "c1cn2ncncc2n1":    "imidazo[1,2-f][1,2,4]triazine",
    "c1cn2nccc2cn1":    "pyrazolo[1,5-a]pyrazine",
    "c1cnn2cnnc2n1":    "[1,2,4]triazolo[4,3-b][1,2,4]triazine",
    # Phase 483: pyrazolo/triazolo fused with pyridazine and pyrazine (IUPAC 2013 P-31.1.3)
    "c1cn2ncnc2cn1":    "[1,2,4]triazolo[1,5-a]pyrazine",
    "c1cnn2nccc2c1":    "pyrazolo[1,5-b]pyridazine",
    "c1cnn2ncnc2c1":    "[1,2,4]triazolo[1,5-b]pyridazine",
    # Phase 484: imidazo/triazolo fused with [1,2,4]triazine and pyridazine (IUPAC 2013 P-31.1.3)
    "c1cnn2cncc2n1":    "imidazo[1,5-b][1,2,4]triazine",
    "c1cnn2ccnc2n1":    "imidazo[3,2-b][1,2,4]triazine",
    "c1cnn2nncc2c1":    "[1,2,3]triazolo[1,5-b]pyridazine",
    # Phase 485: pyrrolo/imidazo/triazolo/pyrazolo fused with pyrazine/pyridazine (IUPAC 2013 P-31.1.3)
    "c1cc2cnccn2c1":    "pyrrolo[1,2-a]pyrazine",
    "c1cnn2cccc2c1":    "pyrrolo[1,2-b]pyridazine",
    "c1cn2cncc2cn1":    "imidazo[1,5-a]pyrazine",
    "c1cnn2cncc2c1":    "imidazo[1,5-b]pyridazine",
    "c1cn2cnnc2cn1":    "[1,2,4]triazolo[4,3-a]pyrazine",
    "c1cnn2nccc2n1":    "pyrazolo[1,5-b][1,2,4]triazine",
    # Phase 486: imidazo[4,5-c]pyridazine, thieno/furo-pyridazine, pyrazolo-pyridine isomers (IUPAC 2013 P-31.1.3)
    "C1=NC2=CCN=NC2=N1": "3H-imidazo[4,5-c]pyridazine",
    "c1cc2ccsc2nn1":    "thieno[2,3-c]pyridazine",
    "c1cc2ccoc2nn1":    "furo[2,3-c]pyridazine",
    "c1cc2n[nH]cc2cn1": "1H-pyrazolo[3,4-c]pyridine",
    # Phase 487: thieno/furo pyridazine [3,4-c] and [x,y-d] isomers (IUPAC 2013 P-31.1.3)
    "c1cc2cscc2nn1":    "thieno[3,4-c]pyridazine",
    "c1cc2cocc2nn1":    "furo[3,4-c]pyridazine",
    "c1cc2cnncc2s1":    "thieno[3,2-d]pyridazine",
    "c1nncc2cscc12":    "thieno[3,4-d]pyridazine",
    "c1cc2cnncc2o1":    "furo[3,2-d]pyridazine",
    "c1nncc2cocc12":    "furo[3,4-d]pyridazine",
    # Phase 488: thieno/furo fused with [1,2,4]triazine at e-bond (IUPAC 2013 P-31.1.3)
    "c1nnc2ccsc2n1":    "thieno[2,3-e][1,2,4]triazine",
    "c1nnc2sccc2n1":    "thieno[3,2-e][1,2,4]triazine",
    "c1nnc2cscc2n1":    "thieno[3,4-e][1,2,4]triazine",
    "c1nnc2ccoc2n1":    "furo[2,3-e][1,2,4]triazine",
    "c1nnc2occc2n1":    "furo[3,2-e][1,2,4]triazine",
    "c1nnc2cocc2n1":    "furo[3,4-e][1,2,4]triazine",
    # Phase 489: isothiazolo/isoxazolo fused with [1,2,4]triazine at e-bond (IUPAC 2013 P-31.1.3)
    "c1nnc2csnc2n1":    "isothiazolo[3,4-e][1,2,4]triazine",
    "c1nnc2nscc2n1":    "isothiazolo[4,3-e][1,2,4]triazine",
    "c1nnc2cnsc2n1":    "isothiazolo[5,4-e][1,2,4]triazine",
    "c1nnc2conc2n1":    "isoxazolo[3,4-e][1,2,4]triazine",
    "c1nnc2nocc2n1":    "isoxazolo[4,3-e][1,2,4]triazine",
    "c1nnc2cnoc2n1":    "isoxazolo[5,4-e][1,2,4]triazine",
    # Phase 501: missing isoxazolo[4,5-e][1,2,4]triazine (IUPAC 2013 P-31.1.3)
    "c1nnc2oncc2n1":    "isoxazolo[4,5-e][1,2,4]triazine",
    # Phase 490: thieno/furo/isothiazolo/isoxazolo fused with [1,2,3]triazine at d-bond (IUPAC 2013 P-31.1.3)
    "c1cc2cnnnc2s1":    "thieno[2,3-d][1,2,3]triazine",
    "c1cc2nnncc2s1":    "thieno[3,2-d][1,2,3]triazine",
    "c1nnnc2cscc12":    "thieno[3,4-d][1,2,3]triazine",
    "c1cc2cnnnc2o1":    "furo[2,3-d][1,2,3]triazine",
    "c1cc2nnncc2o1":    "furo[3,2-d][1,2,3]triazine",
    "c1nnnc2cocc12":    "furo[3,4-d][1,2,3]triazine",
    "c1snc2cnnnc12":    "isothiazolo[4,3-d][1,2,3]triazine",
    "c1nnnc2nscc12":    "isothiazolo[3,4-d][1,2,3]triazine",
    "c1onc2cnnnc12":    "isoxazolo[4,3-d][1,2,3]triazine",
    "c1nnnc2nocc12":    "isoxazolo[3,4-d][1,2,3]triazine",
    # Phase 491: thiazolo/oxazolo fused with [1,2,3]triazine at d-bond (IUPAC 2013 P-31.1.3)
    "c1nc2cnnnc2s1":    "thiazolo[5,4-d][1,2,3]triazine",
    "c1nc2nnncc2s1":    "thiazolo[4,5-d][1,2,3]triazine",
    "c1nc2cnnnc2o1":    "oxazolo[5,4-d][1,2,3]triazine",
    "c1nc2nnncc2o1":    "oxazolo[4,5-d][1,2,3]triazine",
    # Phase 492: 1H-imidazo and 1H-pyrazolo fused with [1,2,3]triazine at d-bond (IUPAC 2013 P-31.1.3)
    "c1nc2nnncc2[nH]1": "1H-imidazo[4,5-d][1,2,3]triazine",
    "c1nc2cnnnc2[nH]1": "1H-imidazo[5,4-d][1,2,3]triazine",
    "c1nnnc2[nH]ncc12": "1H-pyrazolo[5,4-d][1,2,3]triazine",
    "c1nnnc2n[nH]cc12": "1H-pyrazolo[3,4-d][1,2,3]triazine",
    "c1nn[nH]c2cnnc1-2": "1H-pyrazolo[4,5-d][1,2,3]triazine",
    "c1nn[nH]c2cnnc1-2": "1H-pyrazolo[4,3-d][1,2,3]triazine",
    # Phase 493: 1H-imidazo and 1H-pyrazolo fused with [1,2,4]triazine at e-bond (IUPAC 2013 P-31.1.3)
    "c1nnc2nc[nH]c2n1":  "1H-imidazo[5,4-e][1,2,4]triazine",
    "c1nnc2[nH]cnc2n1":  "1H-imidazo[4,5-e][1,2,4]triazine",
    "c1nnc2[nH]ncc2n1":  "1H-pyrazolo[4,5-e][1,2,4]triazine",
    "c1nnc2n[nH]cc2n1":  "1H-pyrazolo[4,3-e][1,2,4]triazine",
    "c1n[nH]c2cnnc-2n1":  "1H-pyrazolo[5,4-e][1,2,4]triazine",
    "c1n[nH]c2cnnc-2n1":  "1H-pyrazolo[3,4-e][1,2,4]triazine",
    # Phase 494: 1H-pyrrolo/pyrazolo/imidazo fused with pyridazine (IUPAC 2013 P-31.1.3)
    "c1cc2cn[nH]c2nn1":  "1H-pyrazolo[5,4-c]pyridazine",
    "c1cc2n[nH]cc2nn1":  "1H-pyrazolo[4,3-c]pyridazine",
    "c1nc2cnncc2[nH]1":  "1H-imidazo[4,5-d]pyridazine",
    "c1nncc2[nH]ncc12":  "1H-pyrazolo[4,5-d]pyridazine",
    "c1nncc2n[nH]cc12":  "1H-pyrazolo[3,4-d]pyridazine",
    "c1cc2cnncc2[nH]1":  "1H-pyrrolo[2,3-d]pyridazine",
    "c1cc2cncc-2[nH]n1":  "1H-pyrrolo[3,4-c]pyridazine",
    "c1cc2cc[nH]c2nn1":  "1H-pyrrolo[2,3-c]pyridazine",
    # Phase 495: 1H-pyrazolo and 1H-pyrrolo fused with pyridine at b- and c-bonds (IUPAC 2013 P-31.1.3)
    "c1cnc2c[nH]nc2c1":  "1H-pyrazolo[4,3-b]pyridine",
    "c1c[nH]c2cncc-2c1":  "1H-pyrrolo[3,4-b]pyridine",
    "c1cc2cn[nH]c2cn1":  "1H-pyrazolo[5,4-c]pyridine",
    "C1=NCc2ccncc21":  "1H-pyrrolo[3,4-c]pyridine",
    # Phase 496: 1H-pyrrolo/[1,2,3]triazolo fused with [1,2,4]triazine at e-bond,
    #            and 1H-pyrrolo[3,4-d]pyridazine (IUPAC 2013 P-31.1.3)
    "c1cc2[nH]ncnc-2n1":  "1H-pyrrolo[2,3-e][1,2,4]triazine",
    "c1n[nH]c2cncc-2n1":  "1H-pyrrolo[3,4-e][1,2,4]triazine",
    "c1nnc2[nH]ccc2n1":  "1H-pyrrolo[3,2-e][1,2,4]triazine",
    "c1nnc2nn[nH]c2n1":  "1H-[1,2,3]triazolo[5,4-e][1,2,4]triazine",
    "c1nnc2[nH]nnc2n1":  "1H-[1,2,3]triazolo[4,5-e][1,2,4]triazine",
    "c1nnc2n[nH]nc2n1":  "2H-[1,2,3]triazolo[4,5-e][1,2,4]triazine",
    "C1=NC=C2CN=NC=C12":  "1H-pyrrolo[3,4-d]pyridazine",
    # Phase 497: 1H-pyrrolo/[1,2,3]triazolo fused with [1,2,3]triazine at d-bond (IUPAC 2013 P-31.1.3)
    "c1cc2[nH]nncc-2n1":  "1H-pyrrolo[3,2-d][1,2,3]triazine",
    "c1ncc2[nH]nncc1-2":  "1H-pyrrolo[3,4-d][1,2,3]triazine",
    "c1cc2cnnnc2[nH]1":  "1H-pyrrolo[2,3-d][1,2,3]triazine",
    "c1nnnc2nn[nH]c12":  "1H-[1,2,3]triazolo[4,5-d][1,2,3]triazine",
    "c1nnnc2[nH]nnc12":  "1H-[1,2,3]triazolo[5,4-d][1,2,3]triazine",
    "c1nnnc2n[nH]nc12":  "2H-[1,2,3]triazolo[4,5-d][1,2,3]triazine",
    # Phase 498: missing oxazolo/thiazolo/isothiazolo/isoxazolo at [1,2,4]triazine e-bond
    #            and isothiazolo/isoxazolo at [1,2,3]triazine d-bond (IUPAC 2013 P-31.1.3)
    "c1nnc2sncc2n1":   "isothiazolo[4,5-e][1,2,4]triazine",
    "c1nnc2ocnc2n1":   "oxazolo[4,5-e][1,2,4]triazine",
    "c1nnc2ncoc2n1":   "oxazolo[5,4-e][1,2,4]triazine",
    "c1nnc2scnc2n1":   "thiazolo[4,5-e][1,2,4]triazine",
    "c1nnc2ncsc2n1":   "thiazolo[5,4-e][1,2,4]triazine",
    "c1nnnc2sncc12":   "isothiazolo[5,4-d][1,2,3]triazine",
    "c1noc2cnnnc12":   "isoxazolo[4,5-d][1,2,3]triazine",
    # Phase 500: remaining isothiazolo/isoxazolo orientations at [1,2,3]triazine d-bond (IUPAC 2013 P-31.1.3)
    "c1nsc2cnnnc12":   "isothiazolo[4,5-d][1,2,3]triazine",
    "c1nnnc2oncc12":   "isoxazolo[5,4-d][1,2,3]triazine",
    # Phase 502: [1,2,3]oxadiazolo fused bicyclics (IUPAC 2013 P-31.1.3)
    "c1cnc2nnoc2c1":   "[1,2,3]oxadiazolo[4,5-b]pyridine",
    "c1cc2onnc2cn1":   "[1,2,3]oxadiazolo[4,5-c]pyridine",
    "c1cnc2onnc2n1":   "[1,2,3]oxadiazolo[4,5-e]pyrazine",
    "c1nncc2onnc12":   "[1,2,3]oxadiazolo[5,4-d]pyridazine",
    "c1nnc2onnc2n1":   "[1,2,3]oxadiazolo[4,5-e][1,2,4]triazine",
    "c1nnnc2onnc12":   "[1,2,3]oxadiazolo[5,4-d][1,2,3]triazine",
    "c1nnnc2nnoc12":   "[1,2,3]oxadiazolo[4,5-d][1,2,3]triazine",
    # Phase 503: [1,2,3]thiadiazolo fused bicyclics (IUPAC 2013 P-31.1.3)
    "c1cnc2nnsc2c1":   "[1,2,3]thiadiazolo[4,5-b]pyridine",
    "c1cc2snnc2cn1":   "[1,2,3]thiadiazolo[4,5-c]pyridine",
    "c1cnc2snnc2n1":   "[1,2,3]thiadiazolo[4,5-e]pyrazine",
    "c1nncc2snnc12":   "[1,2,3]thiadiazolo[5,4-d]pyridazine",
    "c1nnc2snnc2n1":   "[1,2,3]thiadiazolo[4,5-e][1,2,4]triazine",
    "c1nnnc2snnc12":   "[1,2,3]thiadiazolo[5,4-d][1,2,3]triazine",
    "c1nnnc2nnsc12":   "[1,2,3]thiadiazolo[4,5-d][1,2,3]triazine",
    # Phase 504: [1,2,5]oxadiazolo (furazano) fused bicyclics (IUPAC 2013 P-31.1.3)
    "c1cnc2nonc2c1":   "[1,2,5]oxadiazolo[3,4-b]pyridine",
    "c1cnc2nonc2n1":   "[1,2,5]oxadiazolo[3,4-e]pyrazine",
    "c1nncc2nonc12":   "[1,2,5]oxadiazolo[3,4-d]pyridazine",
    "c1nnc2nonc2n1":   "[1,2,5]oxadiazolo[3,4-e][1,2,4]triazine",
    # Phase 505: [1,2,5]thiadiazolo fused bicyclics (IUPAC 2013 P-31.1.3)
    "c1cnc2nsnc2c1":   "[1,2,5]thiadiazolo[3,4-b]pyridine",
    "c1cnc2nsnc2n1":   "[1,2,5]thiadiazolo[3,4-e]pyrazine",
    "c1nncc2nsnc12":   "[1,2,5]thiadiazolo[3,4-d]pyridazine",
    "c1nnc2nsnc2n1":   "[1,2,5]thiadiazolo[3,4-e][1,2,4]triazine",
    # Phase 507: purine-2,6-dione derivatives (IUPAC 2013 P-31.1.7)
    "Cn1c(=O)c2c(ncn2C)n(C)c1=O":                  "1,3,7-trimethyl-3,7-dihydro-1H-purine-2,6-dione",
    "Cn1c(=O)c2[nH]cnc2n(C)c1=O":                  "1,3-dimethyl-3,7-dihydro-1H-purine-2,6-dione",
    "Cn1cnc2c1c(=O)[nH]c(=O)n2C":                  "3,7-dimethyl-3,7-dihydro-1H-purine-2,6-dione",
    "CN1C(=O)NC2=NC=NC21":                          "3-methyl-3,7-dihydro-1H-purine-2,6-dione",
    "CN1NC2N=CN=C2C1=O":                            "1-methyl-3,7-dihydro-1H-purine-2,6-dione",
    "O=c1[nH]cnc2nc[nH]c12":                        "1,7-dihydro-6H-purin-6-one",
    "O=c1[nH]c(=O)c2[nH]c(=O)[nH]c2[nH]1":        "7,9-dihydro-1H-purine-2,6,8(3H)-trione",
    "Nc1nc2[nH]cnc2c(=O)[nH]1":                    "2-amino-3,7-dihydro-1H-purin-6-one",
    # Phase 506: tetrazolo fused bicyclics (IUPAC 2013 P-31.1.3)
    "c1cn2nnnc2cn1":   "tetrazolo[1,5-a]pyrazine",
    "c1cnn2nnnc2c1":   "tetrazolo[1,5-b]pyridazine",
    "c1cnn2nnnc2n1":   "tetrazolo[1,5-b][1,2,4]triazine",
    "c1cc2nnnn2nn1":   "tetrazolo[1,5-f][1,2,3]triazine",
    "c1nncn2nnnc12":   "tetrazolo[1,5-d][1,2,4]triazine",
    "c1ncc2nnnn2n1":   "tetrazolo[1,5-f][1,2,4]triazine",
    # Phase 499: isothiazolo/isoxazolo/oxazolo/thiazolo fused with pyridazine at d-bond (IUPAC 2013 P-31.1.3)
    "c1nncc2sncc12":   "isothiazolo[5,4-d]pyridazine",
    "c1nncc2nscc12":   "isothiazolo[3,4-d]pyridazine",
    "c1nncc2oncc12":   "isoxazolo[5,4-d]pyridazine",
    "c1nncc2nocc12":   "isoxazolo[3,4-d]pyridazine",
    "c1nc2cnncc2o1":   "oxazolo[4,5-d]pyridazine",
    "c1nc2cnncc2s1":   "thiazolo[4,5-d]pyridazine",
    "c1cnc2ccsc2c1":    "thieno[3,2-b]pyridine",   # was mislabelled as [2,3-b]
    "c1cnc2sccc2c1":    "thieno[2,3-b]pyridine",
    "c1cnc2cn[nH]c2c1": "1H-pyrazolo[4,5-b]pyridine",
    "c1cnc2[nH]ccc2c1": "1H-pyrrolo[2,3-b]pyridine",
    "c1cnc2cc[nH]c2c1": "1H-pyrrolo[3,2-b]pyridine",
    "c1cc2ccsc2s1":     "thieno[2,3-b]thiophene",   # was mislabelled as [3,2-b]
    "c1cc2sccc2s1":     "thieno[3,2-b]thiophene",
    "c1cnc2nc[nH]c2c1": "3H-imidazo[4,5-b]pyridine",
    "c1ncc2nc[nH]c2n1": "9H-purine",
    # Phase 141: 単環式ヘテロ芳香族 (Hantzsch-Widman が未対応) IUPAC 2013 P-31.1.3
    # 5-membered, two different heteroatoms
    "c1cnoc1":   "isoxazole",
    "c1cocn1":   "oxazole",
    "c1cscn1":   "thiazole",
    "c1cnsc1":   "isothiazole",
    # 5-membered oxadiazoles
    "c1conn1":   "1,2,3-oxadiazole",
    "c1ncon1":   "1,2,4-oxadiazole",
    "c1cnon1":   "1,2,5-oxadiazole",
    "c1nnco1":   "1,3,4-oxadiazole",
    # 5-membered thiadiazoles
    "c1csnn1":   "1,2,3-thiadiazole",
    "c1ncsn1":   "1,2,4-thiadiazole",
    "c1cnsn1":   "1,2,5-thiadiazole",
    "c1nncs1":   "1,3,4-thiadiazole",
    # 5-membered with adjacent N (pyrazole, triazoles, tetrazole)
    "c1cn[nH]c1":  "1H-pyrazole",
    "Cn1cccn1":    "1-methylpyrazole",
    "c1cn[nH]n1":  "2H-1,2,3-triazole",   # NH between two N's → position 2
    "c1c[nH]nn1":  "1H-1,2,3-triazole",   # NH adjacent to C → position 1
    "c1nc[nH]n1":  "1H-1,2,4-triazole",
    "c1nn[nH]n1":  "1H-tetrazole",
    "c1nnn[nH]1":  "1H-tetrazole",  # alternate canonical form
    # 6-membered triazines and tetrazine
    "c1cnnnc1":    "1,2,3-triazine",
    "c1cnncn1":    "1,2,4-triazine",
    "c1ncncn1":    "1,3,5-triazine",
    "c1nncnn1":    "1,2,4,5-tetrazine",
    # Phase 140: 二ヘテロ原子飽和環 (Hantzsch-Widman, IUPAC 2013 P-31.1.3.4)
    "C1COCCO1":  "1,4-dioxane",
    "C1COCO1":   "1,3-dioxolane",
    "C1COCOC1":  "1,3-dioxane",
    "C1CSCS1":   "1,3-dithiolane",
    "C1CSSC1":   "1,2-dithiolane",
    "C1CSCCS1":  "1,4-dithiane",
    "C1CCSSC1":  "1,2-dithiane",
    "C1CSCCO1":  "1,4-oxathiane",
    "C1CSCO1":   "1,3-oxathiolane",
    # Phase 260: 1,2-ジオキサン・1,2-ジオキソラン (隣接 O-O を持つ環)
    "C1CCOOC1":  "1,2-dioxane",
    "C1COOC1":   "1,2-dioxolane",
    # cyclic carbonates (IUPAC 2013 P-31.1.3)
    "O=C1OCCO1":  "1,3-dioxolan-2-one",
    "O=C1OCCCO1": "1,3-dioxan-2-one",
    # Phase 138: 多環芳香族炭化水素 保留名 (IUPAC 2013 P-31.1.2)
    "c1ccc2cccc-2cc1":                            "azulene",
    "C1=Cc2ccccc2C1":                             "1H-indene",
    "C1=Cc2cccc3cccc1c23":                        "acenaphthylene",
    "c1cc2c3c(cccc3c1)CC2":                       "acenaphthene",
    "c1ccc2c(c1)-c1cccc3cccc-2c13":               "fluoranthene",
    "c1ccc2c(c1)ccc1c3ccccc3ccc21":               "chrysene",
    "c1ccc2c(c1)c1ccccc1c1ccccc21":               "triphenylene",
    "c1cc2ccc3cccc4ccc(c1)c2c34":                 "pyrene",
    "c1cc2cccc3c4cccc5cccc(c(c1)c23)c54":         "perylene",
    "c1ccc2c(c1)cc1ccc3cccc4ccc2c1c34":           "benzo[a]pyrene",
    "c1cc2ccc3ccc4ccc5ccc6ccc1c1c2c3c4c5c61":     "coronene",
    # Phase 441: triazolo-pyrimidine/-pyridazine, pyrazolo-pyrimidine,
    # tetrazolo-pyridine, 1H-imidazo[4,5-b]pyridine (IUPAC 2013 P-31.1.3)
    "c1cnc2ncnn2c1":    "[1,2,4]triazolo[1,5-a]pyrimidine",
    "c1cnn2ncnc2n1":    "[1,2,4]triazolo[1,5-b][1,2,4]triazine",
    "c1cnc2nncn2c1":    "[1,2,4]triazolo[4,3-a]pyrimidine",
    "c1cc2nccnn2c1":    "pyrrolo[1,2-b][1,2,4]triazine",
    "c1ccn2nnnc2c1":    "tetrazolo[1,5-a]pyridine",
    "c1cnc2[nH]cnc2c1": "1H-imidazo[4,5-b]pyridine",
    # Phase 442: imidazo[1,5-a]pyridine, tetrazolo-pyrimidine, pyrrolo-pyridazine,
    # furo/thieno-pyridine isomers, 1H-pyrazolo[3,4-d]pyrimidine (IUPAC 2013 P-31.1.3)
    "c1ccn2cncc2c1":      "imidazo[1,5-a]pyridine",
    "c1cnc2nnnn2c1":      "tetrazolo[1,5-a]pyrimidine",
    "c1cnc2cccn2c1":      "pyrrolo[1,2-a]pyrimidine",
    "c1cc2sccc2cn1":      "thieno[3,2-c]pyridine",
    "c1cnc2ccoc2c1":      "furo[3,2-b]pyridine",
    "c1cnc2cocc2c1":      "furo[3,4-b]pyridine",
    "c1c[nH]c2ncnc-2n1":  "1H-pyrazolo[3,4-d]pyrimidine",
    # Phase 443: pyrazolo/imidazo/pyrrolo-quinoline and pyrazolo-quinoxaline
    "c1ccc2nc3[nH]ncc3cc2c1":  "1H-pyrazolo[3,4-b]quinoline",
    "c1ccc2c(c1)ccc1ccnn12":     "pyrazolo[1,5-a]quinoline",
    "c1ccc2c(c1)ncc1ccnn12":     "pyrazolo[1,5-a]quinoxaline",
    "c1ccc2nc3[nH]cnc3cc2c1":  "1H-imidazo[4,5-b]quinoline",
    "c1ccc2nc3nc[nH]c3cc2c1":  "3H-imidazo[4,5-b]quinoline",
    "c1ccc2nc3cc[nH]c3cc2c1":  "1H-pyrrolo[3,2-b]quinoline",
    # Phase 444: 1H-pyrrolo[2,3-b]quinoline, thieno/furo-quinoline isomers,
    # benzo[g]quinoxaline (IUPAC 2013 P-31.1.3)
    "c1ccc2nc3[nH]ccc3cc2c1":  "1H-pyrrolo[2,3-b]quinoline",
    "c1ccc2nc3sccc3cc2c1":     "thieno[2,3-b]quinoline",
    "c1ccc2nc3cscc3cc2c1":     "thieno[3,4-b]quinoline",
    "c1ccc2nc3ccsc3cc2c1":     "thieno[3,2-b]quinoline",
    "c1ccc2nc3occc3cc2c1":     "furo[2,3-b]quinoline",
    "c1ccc2cc3nccnc3cc2c1":    "benzo[g]quinoxaline",
    # Phase 445: phenanthridine, benzo[h]isoquinoline, phenanthrolines,
    # benzo[f/g]cinnoline (IUPAC 2013 P-31.1.3)
    # benzo[b][1,7]naphthyridine (fix: was mislabeled benzo[c]cinnoline)
    "c1ccc2c(c1)ccc1cccnc12":  "benzo[h]quinoline",
    "c1ccc2c(c1)ccc1ccncc12":  "benzo[h]isoquinoline",
    "c1cnc2c(c1)ccc1cccnc12":  "1,10-phenanthroline",
    "c1cnc2c(c1)ccc1cnccc12":  "1,8-phenanthroline",
    "c1cnc2c(c1)ccc1ncccc12":  "1,7-phenanthroline",
    "c1cnc2cc3ccncc3cc2c1":    "pyrido[3,4-g]quinoline",
    "c1ccc2nc3cnccc3cc2c1":    "benzo[b][1,7]naphthyridine",
    "c1ccc2cc3nnccc3cc2c1":    "benzo[g]cinnoline",
    "c1ccc2c(c1)ccc1ccnnc12":  "benzo[h]cinnoline",
    # Phase 446: benzo[f]quinoxaline, benzo[g/f]phthalazine (IUPAC 2013 P-31.1.3)
    "c1ccc2c(c1)ccc1nccnc12":  "benzo[f]quinoxaline",
    "c1ccc2cc3cnncc3cc2c1":    "benzo[g]phthalazine",
    "c1ccc2c(c1)ccc1cnncc12":  "benzo[f]phthalazine",
    # Phase 447: benzo[h]cinnoline, furo[3,2-b]quinoline, furo[3,4-b]quinoline
    # (IUPAC 2013 P-31.1.3)
    "c1ccc2c(c1)ccc1nnccc12":  "benzo[f]cinnoline",
    "c1ccc2nc3ccoc3cc2c1":     "furo[3,2-b]quinoline",
    "c1ccc2nc3cocc3cc2c1":     "furo[3,4-b]quinoline",
    # Phase 448: benzo[b][1,5/1,6/1,8]naphthyridine, 1H-naphtho[2,3-d]imidazole
    # (IUPAC 2013 P-31.1.2 fusion nomenclature)
    "c1ccc2nc3ncccc3cc2c1":    "benzo[b][1,8]naphthyridine",
    "c1ccc2nc3cccnc3cc2c1":    "benzo[b][1,5]naphthyridine",
    "c1ccc2nc3ccncc3cc2c1":    "benzo[b][1,6]naphthyridine",
    "c1ccc2cc3[nH]cnc3cc2c1":  "1H-naphtho[2,3-d]imidazole",
    # Phase 449: benzo[c]cinnoline (fix), benzo[b][1,7]naphthyridine already in 445,
    # 1H-naphtho[2,3-d]pyrazole, naphtho[2,3-d]oxazole, naphtho[2,3-d]thiazole
    "c1ccc2c(c1)nnc1ccccc12":  "benzo[c]cinnoline",
    "c1ccc2cc3[nH]ncc3cc2c1":  "1H-naphtho[2,3-d]pyrazole",
    "c1ccc2cc3ocnc3cc2c1":     "naphtho[2,3-d]oxazole",
    "c1ccc2cc3scnc3cc2c1":     "naphtho[2,3-d]thiazole",
    # Phase 450: naphtho[2,1-d] series and naphtho[2,3-d]pyrimidine
    # (IUPAC 2013 P-31.1.3 fusion nomenclature)
    "c1ccc2c(c1)ccc1[nH]cnc12":  "1H-naphtho[2,1-d]imidazole",
    "c1ccc2c(c1)ccc1nc[nH]c12":  "1H-naphtho[2,1-d]imidazole",  # 3H tautomeric canonical alias
    "c1ccc2c(c1)ccc1ncoc12":     "naphtho[2,1-d]oxazole",
    "c1ccc2c(c1)ccc1ncsc12":     "naphtho[2,1-d]thiazole",
    "c1ccc2c(c1)ccc1cn[nH]c12":  "1H-naphtho[2,1-d]pyrazole",
    "c1ccc2cc3ncncc3cc2c1":      "naphtho[2,3-d]pyrimidine",
    # Phase 451: naphtho[1,2-d] oxazole/thiazole/pyrimidine and naphtho[1,2-b] furan/thiophene
    # (IUPAC 2013 P-31.1.3 fusion nomenclature)
    "c1ccc2c(c1)ccc1ocnc12":     "naphtho[1,2-d]oxazole",
    "c1ccc2c(c1)ccc1scnc12":     "naphtho[1,2-d]thiazole",
    "c1ccc2c(c1)ccc1ccoc12":     "naphtho[1,2-b]furan",
    "c1ccc2c(c1)ccc1ccsc12":     "naphtho[1,2-b]thiophene",
    "c1ccc2c(c1)ccc1ncncc12":    "naphtho[2,1-d]pyrimidine",
    # Phase 452: naphtho[2,3-b]furan/thiophene and naphtho[2,1-d]pyrimidine
    # (IUPAC 2013 P-31.1.3 fusion nomenclature)
    "c1ccc2cc3occc3cc2c1":        "naphtho[2,3-b]furan",
    "c1ccc2cc3sccc3cc2c1":        "naphtho[2,3-b]thiophene",
    "c1ccc2c(c1)ccc1cncnc12":    "naphtho[1,2-d]pyrimidine",
    # Phase 581: naphtho triazole/oxadiazole/thiadiazole series
    "c1ccc2c(c1)ccc1nn[nH]c12":   "1H-naphtho[1,2-d][1,2,3]triazole",
    "c1ccc2cc3onnc3cc2c1":        "naphtho[2,3-d][1,2,3]oxadiazole",
    "c1ccc2c(c1)ccc1onnc12":      "naphtho[1,2-d][1,2,3]oxadiazole",
    "c1ccc2c(c1)ccc1nnoc12":      "naphtho[2,1-d][1,2,3]oxadiazole",
    "c1ccc2cc3snnc3cc2c1":        "naphtho[2,3-d][1,2,3]thiadiazole",
    "c1ccc2c(c1)ccc1snnc12":      "naphtho[1,2-d][1,2,3]thiadiazole",
    "c1ccc2c(c1)ccc1nnsc12":      "naphtho[2,1-d][1,2,3]thiadiazole",
    # Phase 582: benzo[g]quinoline (= naphtho[2,3-b]pyridine)
    "c1ccc2cc3ncccc3cc2c1":       "benzo[g]quinoline",
    # Phase 580: naphtho isoxazole/isothiazole series
    "c1ccc2cc3oncc3cc2c1":        "naphtho[2,3-d]isoxazole",
    "c1ccc2cc3sncc3cc2c1":        "naphtho[2,3-d]isothiazole",
    "c1ccc2c(c1)ccc1cnoc12":      "naphtho[2,1-d]isoxazole",
    "c1ccc2c(c1)ccc1cnsc12":      "naphtho[2,1-d]isothiazole",
    "c1ccc2c(c1)ccc1oncc12":      "naphtho[1,2-d]isoxazole",
    "c1ccc2c(c1)ccc1sncc12":      "naphtho[1,2-d]isothiazole",
    # Phase 453: pyrido-quinoxaline and pyrido-quinoline tricyclics
    # (IUPAC 2013 P-31.1.3 fusion nomenclature)
    "c1cc2cc3nccnc3cc2cn1":      "pyrido[3,4-g]quinoxaline",
    "c1cnc2cc3nccnc3cc2c1":      "pyrido[2,3-g]quinoxaline",
    "c1ccc2nc3nccnc3cc2c1":      "pyrazino[2,3-b]quinoline",
    "c1cnc2c(c1)cnc1nccnc12":      "pyrazino[2,3-h][1,6]naphthyridine",
    "c1cnc2cc3ncccc3cc2c1":      "pyrido[3,2-g]quinoline",
    # Phase 454: benzo/pyrazino fused naphthyridines and pyrazinoquinoxaline
    "c1cnc2cc3cccnc3cc2c1":      "pyrido[2,3-g]quinoline",
    "c1cnc2cc3nccnc3nc2c1":      "pyrazino[2,3-b][1,5]naphthyridine",
    "c1cnc2nc3nccnc3cc2c1":      "pyrazino[2,3-g][1,8]naphthyridine",
    "c1ccc2nc3nccnc3nc2c1":      "pyrazino[2,3-b]quinoxaline",
    "c1cnc2cc3nccnc3cc2n1":      "pyrazino[2,3-g]quinoxaline",
    # Phase 455: benzo/pyrido fused [1,6]-naphthyridines and pyrido-naphthyridines
    "c1ccc2nc3ncccc3nc2c1":      "pyrido[2,3-b]quinoxaline",
    # Phase 456: pyrido[4,3-g]quinoline, pyrido[3,2-g]quinazoline, pyrimido[4,5-b]quinoline
    "c1cnc2cc3cnccc3cc2c1":      "pyrido[4,3-g]quinoline",
    "c1cnc2cc3ncncc3cc2c1":      "pyrido[3,2-g]quinazoline",
    "c1ccc2nc3ncncc3cc2c1":      "pyrimido[4,5-b]quinoline",
    # Phase 457: pyrido[2,3-g]quinazoline, pyrimido[5,4-g/b] series, pyrido[2,3-g][1,5]naphthyridine
    "c1cnc2cc3cncnc3cc2c1":      "pyrido[2,3-g]quinazoline",
    "c1cnc2cc3ncncc3nc2c1":      "pyrimido[5,4-g][1,5]naphthyridine",
    "c1ccc2nc3cncnc3cc2c1":      "pyrimido[5,4-b]quinoline",
    "c1cnc2nc3ncncc3cc2c1":      "pyrimido[5,4-g][1,8]naphthyridine",
    "c1cnc2nc3cccnc3cc2c1":      "pyrido[2,3-g][1,5]naphthyridine",
    # Phase 458: pyrimido[4,5-b]quinoxaline, pyrido[2,3-b][1,6/7/8]naphthyridines
    "c1ccc2nc3ncncc3nc2c1":      "pyrimido[4,5-b]quinoxaline",
    "c1cnc2nc3ncccc3cc2c1":      "pyrido[2,3-b][1,8]naphthyridine",
    "c1cnc2nc3ccncc3cc2c1":      "pyrido[2,3-b][1,6]naphthyridine",
    "c1cnc2nc3cnccc3cc2c1":      "pyrido[2,3-b][1,7]naphthyridine",
    # Phase 459: benzo[c][1,6]naphthyridine, benzo[c][1,7]naphthyridine
    "c1ccc2c(c1)cnc1ccncc12":    "benzo[c][1,6]naphthyridine",
    "c1ccc2c(c1)cnc1cnccc12":    "benzo[c][1,7]naphthyridine",
    # Phase 460: phenanthroline isomers (IUPAC 2013 P-31.1.3.4)
    "c1cnc2ccc3ccncc3c2c1":      "2,7-phenanthroline",
    "c1ccc2c(c1)ncc1ccncc12":    "2,6-phenanthroline",
    "c1cnc2c(c1)cnc1ccccc12":    "1,6-phenanthroline",
    "c1ccc2c(c1)ncc1cnccc12":    "3,6-phenanthroline",
    "c1ccc2c(c1)cnc1cnccc12":    "3,5-phenanthroline",
    "c1ccc2c(c1)cnc1ncccc12":    "4,5-phenanthroline",
    # Phase 461: benzo-fused 5-membered thiadiazole/oxadiazole/isothiazole (IUPAC 2013 P-31.1.3)
    "c1ccc2snnc2c1":    "1,2,3-benzothiadiazole",
    "c1ccc2nsnc2c1":    "2,1,3-benzothiadiazole",
    "c1ccc2nonc2c1":    "2,1,3-benzoxadiazole",
    "c1ccc2sncc2c1":    "1,2-benzisothiazole",
    "c1ccc2nscc2c1":    "2,1-benzisothiazole",
    # Phase 462: benzo[c]thiophene, benzo[b]selenophene, benzo[c]selenophene,
    #            1,2,3-/2,1,3-benzoselenadiazole (IUPAC 2013 P-31.1.3)
    "c1ccc2cscc2c1":    "benzo[c]thiophene",
    "c1ccc2[se]ccc2c1": "benzo[b]selenophene",
    "c1ccc2c[se]cc2c1": "benzo[c]selenophene",
    "c1ccc2[se]nnc2c1": "1,2,3-benzoselenadiazole",
    "c1ccc2n[se]nc2c1": "2,1,3-benzoselenadiazole",
    # Phase 463: 1,2,4-benzotriazine and pyrido[x,y-z]pyridazine isomers
    #            (IUPAC 2013 P-31.1.3 fusion nomenclature)
    "c1ccc2nncnc2c1":   "1,2,4-benzotriazine",
    "c1cnc2ccnnc2c1":   "pyrido[2,3-e]pyridazine",
    "c1cnc2cnncc2c1":   "pyrido[2,3-d]pyridazine",
    "c1cc2ccnnc2cn1":   "pyrido[3,4-c]pyridazine",
    # Phase 464: remaining pyrido[x,y-e]pyridazine/pyrimidine/pyrazine isomers
    #            (IUPAC 2013 P-31.1.3 fusion nomenclature)
    "c1cc2nnccc2cn1":   "pyrido[3,4-e]pyridazine",
    "c1cnc2nnccc2c1":   "pyrido[2,3-c]pyridazine",
    "c1cnc2cncnc2c1":   "pyrido[2,3-e]pyrimidine",
    "c1cc2ncncc2cn1":   "pyrido[3,4-e]pyrimidine",
    "c1cc2nccnc2cn1":   "pyrido[3,4-e]pyrazine",
    # Phase 465: furo[3,2-c]pyridine, furo[2,3-e]pyridazine,
    #            1H-pyrazolo[4,5-c]pyridine, 1H-pyrrolo[3,2-c]pyridine
    #            (IUPAC 2013 P-31.1.3 fusion nomenclature)
    "c1cc2occc2cn1":    "furo[3,2-c]pyridine",
    "c1cc2occc2nn1":    "furo[2,3-e]pyridazine",
    "c1cc2[nH]ncc2cn1": "1H-pyrazolo[4,5-c]pyridine",
    "c1cc2[nH]ccc2cn1": "1H-pyrrolo[3,2-c]pyridine",
    # Phase 466: thieno/furo/pyrrolo-[x,y-c]pyridine isomers,
    #            thieno/pyrrolo-[x,y-c]pyridazine (IUPAC 2013 P-31.1.3)
    "c1cc2ccsc2cn1":    "thieno[2,3-c]pyridine",
    "c1cc2cscc2cn1":    "thieno[3,4-c]pyridine",
    "c1cc2ccoc2cn1":    "furo[2,3-c]pyridine",
    "c1cc2cocc2cn1":    "furo[3,4-c]pyridine",
    "c1cc2cc[nH]c2cn1": "1H-pyrrolo[2,3-c]pyridine",
    "c1cc2sccc2nn1":    "thieno[3,2-c]pyridazine",
    "c1cc2nccc-2[nH]n1": "1H-pyrrolo[3,2-c]pyridazine",
    # Phase 467: thiazolo/oxazolo/imidazo[x,y-c]pyridine,
    #            1H-[1,2,3]triazolo[4,5-b]pyridine,
    #            pyrazolo[3,4-c] and [4,5-c] pyridine/pyridazine isomers
    #            (IUPAC 2013 P-31.1.3 fusion nomenclature)
    "c1cc2ncsc2cn1":    "thiazolo[5,4-c]pyridine",
    "c1cc2ncoc2cn1":    "oxazolo[5,4-c]pyridine",
    "c1cc2[nH]cnc2cn1": "1H-imidazo[4,5-c]pyridine",
    "c1cnc2nn[nH]c2c1": "1H-[1,2,3]triazolo[4,5-b]pyridine",
    "c1cc2[nH]ncc2nn1": "1H-pyrazolo[4,5-c]pyridazine",
    "c1cc2c[nH]nc2cn1": "1H-pyrazolo[3,4-c]pyridine",
    "c1cc2c[nH]nc2nn1": "1H-pyrazolo[3,4-c]pyridazine",
    # Phase 468: thiazolo[4,5-c], oxazolo[4,5-c], isoxazolo-c pyridine isomers,
    #            2H-[1,2,3]triazolo b/c pyridine, thiazolo-c-pyridazine isomers
    #            (IUPAC 2013 P-31.1.3 fusion nomenclature)
    "c1cc2scnc2cn1":    "thiazolo[4,5-c]pyridine",
    "c1cc2n[nH]nc2cn1": "2H-[1,2,3]triazolo[4,5-c]pyridine",
    "c1cc2nn[nH]c2cn1": "1H-[1,2,3]triazolo[5,4-c]pyridine",
    "c1cnc2n[nH]nc2c1": "2H-[1,2,3]triazolo[4,5-b]pyridine",
    "c1cc2scnc2nn1":    "thiazolo[4,5-c]pyridazine",
    "c1cc2ncsc2nn1":    "thiazolo[5,4-c]pyridazine",
    "c1cc2ocnc2cn1":    "oxazolo[4,5-c]pyridine",
    "c1cc2cnoc2cn1":    "isoxazolo[5,4-c]pyridine",
    "c1cc2oncc2cn1":    "isoxazolo[4,5-c]pyridine",
    "c1cc2nocc2cn1":    "isoxazolo[4,3-c]pyridine",
    # Phase 469: oxazolo/isoxazolo-c-pyridazine, imidazo[4,5-c]pyridazine,
    #            [1,2,3]triazolo-c/b-pyridazine isomers
    #            (IUPAC 2013 P-31.1.3 fusion nomenclature)
    "c1cc2ncoc2nn1":    "oxazolo[5,4-c]pyridazine",
    "c1cc2ocnc2nn1":    "oxazolo[4,5-c]pyridazine",
    "c1cc2cnoc2nn1":    "isoxazolo[5,4-c]pyridazine",
    "c1cc2oncc2nn1":    "isoxazolo[4,5-c]pyridazine",
    "c1cc2nocc2nn1":    "isoxazolo[4,3-c]pyridazine",
    # Phase 479: missing isoxazolo isomers: [3,4-c]pyridine/pyridazine, [4,3-b]pyridine
    "c1cc2conc2cn1":    "isoxazolo[3,4-c]pyridine",
    "c1cc2conc2nn1":    "isoxazolo[3,4-c]pyridazine",
    "c1cc2[nH]cnc2nn1": "1H-imidazo[4,5-c]pyridazine",
    "c1cc2[nH]nnc2nn1": "1H-[1,2,3]triazolo[4,5-c]pyridazine",
    "c1cc2n[nH]nc2nn1": "2H-[1,2,3]triazolo[4,5-c]pyridazine",
    "c1cc2nn[nH]c2nn1": "1H-[1,2,3]triazolo[5,4-c]pyridazine",
    # corrected Phase 474: ring6 is pyrazine (c1cnc2 template), not pyridazine
    "c1cnc2[nH]nnc2n1": "1H-[1,2,3]triazolo[4,5-e]pyrazine",
    "c1cnc2n[nH]nc2n1": "2H-[1,2,3]triazolo[4,5-e]pyrazine",
    # Phase 470: fix Phase 439 b-bond pyridine name swaps (oxazolo/thiazolo [4,5] vs [5,4]),
    #            fix Phase 439 pyrimidine→pyrazine errors, add real pyrimidine d-bond fusions,
    #            isoxazolo-b-pyridine, 1H-pyrazolo[3,4-b]pyridine
    #            (IUPAC 2013 P-31.1.3 fusion nomenclature)
    "c1ncc2ncoc2n1":    "oxazolo[5,4-d]pyrimidine",
    "c1ncc2ocnc2n1":    "oxazolo[4,5-d]pyrimidine",
    "c1ncc2ncsc2n1":    "thiazolo[5,4-d]pyrimidine",
    "c1ncc2scnc2n1":    "thiazolo[4,5-d]pyrimidine",
    "c1cnc2nocc2c1":    "isoxazolo[3,4-b]pyridine",
    "c1cnc2cnoc2c1":    "isoxazolo[4,5-b]pyridine",
    "c1cnc2oncc2c1":    "isoxazolo[5,4-b]pyridine",
    "c1cnc2conc2c1":    "isoxazolo[4,3-b]pyridine",
    "c1cnc2n[nH]cc2c1": "1H-pyrazolo[3,4-b]pyridine",
    # Phase 471: thieno/furo/pyrrolo[x,y-d]pyrimidine (IUPAC 2013 P-31.1.3)
    "c1ncc2ccsc2n1":     "thieno[2,3-d]pyrimidine",
    "c1ncc2cscc2n1":     "thieno[3,4-d]pyrimidine",
    "c1ncc2sccc2n1":     "thieno[3,2-d]pyrimidine",
    "c1ncc2ccoc2n1":     "furo[2,3-d]pyrimidine",
    "c1ncc2cocc2n1":     "furo[3,4-d]pyrimidine",
    "c1ncc2occc2n1":     "furo[3,2-d]pyrimidine",
    "c1ncc2cc[nH]c2n1":  "1H-pyrrolo[2,3-d]pyrimidine",
    "c1ncc2cncc-2[nH]1":  "1H-pyrrolo[3,4-d]pyrimidine",
    "c1cc2[nH]cncc-2n1":  "1H-pyrrolo[3,2-d]pyrimidine",
    # Phase 472: isoxazolo/pyrazolo/[1,2,3]triazolo[x,y-d]pyrimidine (IUPAC 2013 P-31.1.3)
    "c1ncc2nocc2n1":     "isoxazolo[4,3-d]pyrimidine",
    "c1ncc2oncc2n1":     "isoxazolo[4,5-d]pyrimidine",
    "c1ncc2conc2n1":     "isoxazolo[3,4-d]pyrimidine",
    "c1ncc2cnoc2n1":     "isoxazolo[5,4-d]pyrimidine",
    "c1ncc2[nH]ncc2n1":  "1H-pyrazolo[4,5-d]pyrimidine",
    "c1ncc2n[nH]cc2n1":  "1H-pyrazolo[4,3-d]pyrimidine",
    "c1ncc2cn[nH]c2n1":  "1H-pyrazolo[5,4-d]pyrimidine",
    "c1ncc2c[nH]nc2n1":  "1H-pyrazolo[3,4-d]pyrimidine",
    "c1ncc2[nH]nnc2n1":  "1H-[1,2,3]triazolo[4,5-d]pyrimidine",
    "c1ncc2n[nH]nc2n1":  "2H-[1,2,3]triazolo[4,5-d]pyrimidine",
    "c1ncc2nn[nH]c2n1":  "1H-[1,2,3]triazolo[5,4-d]pyrimidine",
}

# atom index in canonical base SMILES → IUPAC locant (None = ring junction, skip)
_FUSED_LOCANT_MAP: dict[str, dict[int, int | None]] = {
    # Phase 612: 5-atom monocyclic heteroaromatics needing locant maps
    "c1cnon1": {0: 3, 1: 3, 2: None, 3: None, 4: None},  # 1,2,5-oxadiazole (symm; C3=C4→loc3)
    "c1cnsn1": {0: 3, 1: 3, 2: None, 3: None, 4: None},  # 1,2,5-thiadiazole (symm)
    "c1conn1": {0: 4, 1: 5, 2: None, 3: None, 4: None},  # 1,2,3-oxadiazole
    "c1csnn1": {0: 4, 1: 5, 2: None, 3: None, 4: None},  # 1,2,3-thiadiazole
    "c1ncon1": {0: 3, 1: None, 2: 5, 3: None, 4: None},  # 1,2,4-oxadiazole
    "c1ncsn1": {0: 3, 1: None, 2: 5, 3: None, 4: None},  # 1,2,4-thiadiazole
    "c1nnco1": {0: 2, 1: None, 2: None, 3: 2, 4: None},  # 1,3,4-oxadiazole (symm)
    "c1nncs1": {0: 2, 1: None, 2: None, 3: 2, 4: None},  # 1,3,4-thiadiazole (symm)
    "c1cc[se]c1": {0: 3, 1: 3, 2: 2, 3: None, 4: 2},  # selenophene (symm; 2,3)
    "c1cc[te]c1": {0: 3, 1: 3, 2: 2, 3: None, 4: 2},  # tellurophene (symm; 2,3)
    # Phase 613: 5- and 6-atom monocyclic heteroaromatics (triazines, tetrazine, isothiazole, isoxazole, oxazole, thiazole)
    "c1cnnnc1": {0: 5, 1: 4, 2: None, 3: None, 4: None, 5: 4},  # 1,2,3-triazine (C4=C6 symm)
    "c1nncnn1": {0: 3, 1: None, 2: None, 3: 3, 4: None, 5: None},  # 1,2,4,5-tetrazine (symm)
    "c1cnncn1": {0: 5, 1: 6, 2: None, 3: None, 4: 3, 5: None},  # 1,2,4-triazine
    "c1ncncn1": {0: 2, 1: None, 2: 2, 3: None, 4: 2, 5: None},  # 1,3,5-triazine (symm)
    "c1cnsc1": {0: 4, 1: 3, 2: None, 3: None, 4: 5},  # isothiazole
    "c1cnoc1": {0: 4, 1: 3, 2: None, 3: None, 4: 5},  # isoxazole
    "c1cocn1": {0: 4, 1: 5, 2: None, 3: 2, 4: None},  # oxazole
    "c1cscn1": {0: 4, 1: 5, 2: None, 3: 2, 4: None},  # thiazole
    # Phase 613: thieno[2,3-b]thiophene naming fix + new thieno[3,2-b]thiophene
    "c1cc2ccsc2s1": {0: 2, 1: 3, 2: None, 3: 3, 4: 2, 5: None, 6: None, 7: None},  # thieno[2,3-b]thiophene (symm 2=5, 3=6)
    "c1cc2sccc2s1": {0: 2, 1: 3, 2: None, 3: None, 4: 2, 5: 3, 6: None, 7: None},  # thieno[3,2-b]thiophene (symm 2=5, 3=6)
    # Phase 613: thieno[3,2-b]pyridine naming fix + new thieno[2,3-b]pyridine
    "c1cnc2ccsc2c1": {0: 6, 1: 5, 2: None, 3: None, 4: 3, 5: 2, 6: None, 7: None, 8: 7},  # thieno[3,2-b]pyridine
    "c1cnc2sccc2c1": {0: 5, 1: 6, 2: None, 3: None, 4: None, 5: 2, 6: 3, 7: None, 8: 4},  # thieno[2,3-b]pyridine
    # Phase 613: azulene (5-7 fused, C2v; 3a/8a junctions), chrysene (4-ring PAH, D2h)
    "c1ccc2cccc-2cc1": {0: 6, 1: 5, 2: 4, 3: None, 4: 1, 5: 2, 6: 1, 7: None, 8: 4, 9: 5},  # azulene (symm: 1=3, 4=8, 5=7)
    "c1ccc2c(c1)ccc1c3ccccc3ccc21": {0: 2, 1: 3, 2: 4, 3: None, 4: None, 5: 1, 6: 6, 7: 5, 8: None, 9: None, 10: 4, 11: 3, 12: 2, 13: 1, 14: None, 15: 6, 16: 5, 17: None},  # chrysene (symm: 1=7, 2=8, 3=9, 4=10, 5=11, 6=12)
    # Phase 614: 13-atom tricyclics (thieno/furo-quinoline, pyrazolo-quinoline/quinoxaline)
    "c1ccc2nc3occc3cc2c1": {0: 6, 1: 7, 2: 8, 3: None, 4: None, 5: None, 6: None, 7: 2, 8: 3, 9: None, 10: 4, 11: None, 12: 5},  # furo[2,3-b]quinoline
    "c1ccc2nc3ccoc3cc2c1": {0: 7, 1: 6, 2: 5, 3: None, 4: None, 5: None, 6: 3, 7: 2, 8: None, 9: None, 10: 9, 11: None, 12: 8},  # furo[3,2-b]quinoline
    "c1ccc2nc3cocc3cc2c1": {0: 7, 1: 6, 2: 5, 3: None, 4: None, 5: None, 6: 3, 7: None, 8: 1, 9: None, 10: 9, 11: None, 12: 8},  # furo[3,4-b]quinoline
    "c1ccc2nc3sccc3cc2c1": {0: 6, 1: 7, 2: 8, 3: None, 4: None, 5: None, 6: None, 7: 2, 8: 3, 9: None, 10: 4, 11: None, 12: 5},  # thieno[2,3-b]quinoline
    "c1ccc2nc3ccsc3cc2c1": {0: 7, 1: 6, 2: 5, 3: None, 4: None, 5: None, 6: 3, 7: 2, 8: None, 9: None, 10: 9, 11: None, 12: 8},  # thieno[3,2-b]quinoline
    "c1ccc2nc3cscc3cc2c1": {0: 7, 1: 6, 2: 5, 3: None, 4: None, 5: None, 6: 3, 7: None, 8: 1, 9: None, 10: 9, 11: None, 12: 8},  # thieno[3,4-b]quinoline
    "c1ccc2c(c1)ccc1ccnn12":  {0: 7, 1: 8, 2: 9, 3: None, 4: None, 5: 6, 6: 5, 7: 4, 8: None, 9: 3, 10: 2, 11: None, 12: None},  # pyrazolo[1,5-a]quinoline
    "c1ccc2c(c1)ncc1ccnn12":  {0: 7, 1: 8, 2: 9, 3: None, 4: None, 5: 6, 6: None, 7: 4, 8: None, 9: 3, 10: 2, 11: None, 12: None},  # pyrazolo[1,5-a]quinoxaline
    # Phase 616: simple NH monocyclics
    "c1cn[nH]c1": {0: 4, 1: 3, 2: None, 3: None, 4: 5},  # 1H-pyrazole
    "c1c[nH]nn1": {0: 4, 1: 5, 2: None, 3: None, 4: None},  # 1H-1,2,3-triazole
    "c1nc[nH]n1": {0: 3, 1: None, 2: 5, 3: None, 4: None},  # 1H-1,2,4-triazole
    "c1nn[nH]n1": {0: 5, 1: None, 2: None, 3: None, 4: None},  # 1H-tetrazole (atom 0 = only C)
    "c1nnn[nH]1": {0: 5, 1: None, 2: None, 3: None, 4: None},  # 1H-tetrazole (alternate canonical)
    "c1cn[nH]n1": {0: 4, 1: 4, 2: None, 3: None, 4: None},  # 2H-1,2,3-triazole (C4=C5 symmetric)
    # Phase 617: NH-containing bicyclic heteroaromatics
    "c1cnc2[nH]nnc2c1": {0: 6, 1: 5, 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: 7},  # 1H-[1,2,3]triazolo[4,5-b]pyridine
    "c1cnc2nn[nH]c2c1": {0: 6, 1: 5, 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: 7},  # 1H-[1,2,3]triazolo[4,5-b]pyridine
    "c1cc2[nH]nnc2nn1": {0: 6, 1: 7, 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: None},  # 1H-[1,2,3]triazolo[4,5-c]pyridazine
    "c1nnnc2nn[nH]c12": {0: 7, 1: None, 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: None},  # 1H-[1,2,3]triazolo[4,5-d][1,2,3]triazine
    "c1ncc2[nH]nnc2n1": {0: 5, 1: None, 2: 7, 3: None, 4: None, 5: None, 6: None, 7: None, 8: None},  # 1H-[1,2,3]triazolo[4,5-d]pyrimidine
    "c1nnc2[nH]nnc2n1": {0: 6, 1: None, 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: None},  # 1H-[1,2,3]triazolo[4,5-e][1,2,4]triazine
    "c1cnc2[nH]nnc2n1": {0: 5, 1: 5, 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: None},  # 1H-[1,2,3]triazolo[4,5-e]pyrazine
    "c1cc2nn[nH]c2nn1": {0: 6, 1: 7, 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: None},  # 1H-[1,2,3]triazolo[5,4-c]pyridazine
    "c1cc2nn[nH]c2cn1": {0: 6, 1: 7, 2: None, 3: None, 4: None, 5: None, 6: None, 7: 4, 8: None},  # 1H-[1,2,3]triazolo[5,4-c]pyridine
    "c1nnnc2[nH]nnc12": {0: 7, 1: None, 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: None},  # 1H-[1,2,3]triazolo[5,4-d][1,2,3]triazine
    "c1ncc2nn[nH]c2n1": {0: 5, 1: None, 2: 7, 3: None, 4: None, 5: None, 6: None, 7: None, 8: None},  # 1H-[1,2,3]triazolo[5,4-d]pyrimidine
    "c1nnc2nn[nH]c2n1": {0: 6, 1: None, 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: None},  # 1H-[1,2,3]triazolo[5,4-e][1,2,4]triazine
    "c1cnc2[nH]cnc2c1": {0: 6, 1: 5, 2: None, 3: None, 4: None, 5: 2, 6: None, 7: None, 8: 7},  # 1H-imidazo[4,5-b]pyridine
    "c1ccc2nc3[nH]cnc3cc2c1": {0: 7, 1: 6, 2: 5, 3: None, 4: None, 5: None, 6: None, 7: 2, 8: None, 9: None, 10: 9, 11: None, 12: 8},  # 1H-imidazo[4,5-b]quinoline
    "c1cc2[nH]cnc2nn1": {0: 3, 1: 4, 2: None, 3: None, 4: 6, 5: None, 6: None, 7: None, 8: None},  # 1H-imidazo[4,5-c]pyridazine
    "c1cc2nc[nH]c2cn1": {0: 6, 1: 7, 2: None, 3: None, 4: 2, 5: None, 6: None, 7: 4, 8: None},  # 1H-imidazo[4,5-c]pyridine
    "c1cc2[nH]cnc2cn1": {0: 6, 1: 7, 2: None, 3: None, 4: 2, 5: None, 6: None, 7: 4, 8: None},  # 1H-imidazo[4,5-c]pyridine
    "c1nc2nnncc2[nH]1": {0: 6, 1: None, 2: None, 3: None, 4: None, 5: None, 6: 4, 7: None, 8: None},  # 1H-imidazo[4,5-d][1,2,3]triazine
    "c1nc2cnncc2[nH]1": {0: 2, 1: None, 2: None, 3: 4, 4: None, 5: None, 6: 4, 7: None, 8: None},  # 1H-imidazo[4,5-d]pyridazine
    "c1nnc2[nH]cnc2n1": {0: 3, 1: None, 2: None, 3: None, 4: None, 5: 6, 6: None, 7: None, 8: None},  # 1H-imidazo[4,5-e][1,2,4]triazine
    "c1cnc2[nH]cnc2n1": {0: 5, 1: 5, 2: None, 3: None, 4: None, 5: 2, 6: None, 7: None, 8: None},  # 1H-imidazo[4,5-e]pyrazine
    "c1nc2cnnnc2[nH]1": {0: 6, 1: None, 2: None, 3: 4, 4: None, 5: None, 6: None, 7: None, 8: None},  # 1H-imidazo[5,4-d][1,2,3]triazine
    "c1nnc2nc[nH]c2n1": {0: 3, 1: None, 2: None, 3: None, 4: None, 5: 6, 6: None, 7: None, 8: None},  # 1H-imidazo[5,4-e][1,2,4]triazine
    "c1cnc2[nH]ncc2c1": {0: 5, 1: 6, 2: None, 3: None, 4: None, 5: None, 6: 3, 7: None, 8: 4},  # 1H-pyrazolo[3,4-b]pyridine
    "c1cnc2n[nH]cc2c1": {0: 5, 1: 6, 2: None, 3: None, 4: None, 5: None, 6: 3, 7: None, 8: 4},  # 1H-pyrazolo[3,4-b]pyridine
    "c1ccc2nc3[nH]ncc3cc2c1": {0: 6, 1: 7, 2: 8, 3: None, 4: None, 5: None, 6: None, 7: None, 8: 3, 9: None, 10: 4, 11: None, 12: 5},  # 1H-pyrazolo[3,4-b]quinoline
    "c1cc2c[nH]nc2nn1": {0: 5, 1: 4, 2: None, 3: 3, 4: None, 5: None, 6: None, 7: None, 8: None},  # 1H-pyrazolo[3,4-c]pyridazine
    "c1cc2c[nH]nc2cn1": {0: 5, 1: 4, 2: None, 3: 3, 4: None, 5: None, 6: None, 7: 7, 8: None},  # 1H-pyrazolo[3,4-c]pyridine
    "c1nnnc2n[nH]cc12": {0: 4, 1: None, 2: None, 3: None, 4: None, 5: None, 6: None, 7: 5, 8: None},  # 1H-pyrazolo[3,4-d][1,2,3]triazine
    "c1nncc2n[nH]cc12": {0: 4, 1: None, 2: None, 3: 7, 4: None, 5: None, 6: None, 7: 3, 8: None},  # 1H-pyrazolo[3,4-d]pyridazine
    "c1ncc2c[nH]nc2n1": {0: 6, 1: None, 2: 4, 3: None, 4: 3, 5: None, 6: None, 7: None, 8: None},  # 1H-pyrazolo[3,4-d]pyrimidine
    "c1n[nH]c2cnnc-2n1": {0: 3, 1: None, 2: None, 3: None, 4: 7, 5: None, 6: None, 7: None, 8: None},  # 1H-pyrazolo[3,4-e][1,2,4]triazine
    "c1cnc2n[nH]cc2n1": {0: 5, 1: 6, 2: None, 3: None, 4: None, 5: None, 6: 3, 7: None, 8: None},  # 1H-pyrazolo[3,4-e]pyrazine
    "c1cnc2c[nH]nc2c1": {0: 6, 1: 5, 2: None, 3: None, 4: 3, 5: None, 6: None, 7: None, 8: 7},  # 1H-pyrazolo[4,3-b]pyridine
    "c1cc2n[nH]cc2nn1": {0: 6, 1: 7, 2: None, 3: None, 4: None, 5: 3, 6: None, 7: None, 8: None},  # 1H-pyrazolo[4,3-c]pyridazine
    "c1nn[nH]c2cnnc1-2": {0: 4, 1: None, 2: None, 3: None, 4: None, 5: 7, 6: None, 7: None, 8: None},  # 1H-pyrazolo[4,3-d][1,2,3]triazine
    "c1ncc2n[nH]cc2n1": {0: 5, 1: None, 2: 7, 3: None, 4: None, 5: None, 6: 3, 7: None, 8: None},  # 1H-pyrazolo[4,3-d]pyrimidine
    "c1nnc2n[nH]cc2n1": {0: 5, 1: None, 2: None, 3: None, 4: None, 5: None, 6: 3, 7: None, 8: None},  # 1H-pyrazolo[4,3-e][1,2,4]triazine
    "c1cnc2cn[nH]c2c1": {0: 6, 1: 5, 2: None, 3: None, 4: 3, 5: None, 6: None, 7: None, 8: 7},  # 1H-pyrazolo[4,5-b]pyridine
    "c1cc2[nH]ncc2nn1": {0: 6, 1: 7, 2: None, 3: None, 4: None, 5: 3, 6: None, 7: None, 8: None},  # 1H-pyrazolo[4,5-c]pyridazine
    "c1cc2[nH]ncc2cn1": {0: 6, 1: 7, 2: None, 3: None, 4: None, 5: 3, 6: None, 7: 4, 8: None},  # 1H-pyrazolo[4,5-c]pyridine
    "c1nncc2[nH]ncc12": {0: 4, 1: None, 2: None, 3: 7, 4: None, 5: None, 6: None, 7: 3, 8: None},  # 1H-pyrazolo[4,5-d]pyridazine
    "c1ncc2[nH]ncc2n1": {0: 5, 1: None, 2: 7, 3: None, 4: None, 5: None, 6: 3, 7: None, 8: None},  # 1H-pyrazolo[4,5-d]pyrimidine
    "c1nnc2[nH]ncc2n1": {0: 5, 1: None, 2: None, 3: None, 4: None, 5: None, 6: 3, 7: None, 8: None},  # 1H-pyrazolo[4,5-e][1,2,4]triazine
    "c1cnc2[nH]ncc2n1": {0: 5, 1: 6, 2: None, 3: None, 4: None, 5: None, 6: 3, 7: None, 8: None},  # 1H-pyrazolo[4,5-e]pyrazine
    "c1cc2cn[nH]c2nn1": {0: 5, 1: 4, 2: None, 3: 3, 4: None, 5: None, 6: None, 7: None, 8: None},  # 1H-pyrazolo[5,4-c]pyridazine
    "c1cc2cn[nH]c2cn1": {0: 5, 1: 4, 2: None, 3: 3, 4: None, 5: None, 6: None, 7: 7, 8: None},  # 1H-pyrazolo[5,4-c]pyridine
    "c1nnnc2[nH]ncc12": {0: 4, 1: None, 2: None, 3: None, 4: None, 5: None, 6: None, 7: 5, 8: None},  # 1H-pyrazolo[5,4-d][1,2,3]triazine
    "c1ncc2cn[nH]c2n1": {0: 6, 1: None, 2: 4, 3: None, 4: 3, 5: None, 6: None, 7: None, 8: None},  # 1H-pyrazolo[5,4-d]pyrimidine
    "c1cnc2[nH]ccc2c1": {0: 5, 1: 6, 2: None, 3: None, 4: None, 5: 2, 6: 3, 7: None, 8: 4},  # 1H-pyrrolo[2,3-b]pyridine
    "c1ccc2nc3[nH]ccc3cc2c1": {0: 6, 1: 7, 2: 8, 3: None, 4: None, 5: None, 6: None, 7: 2, 8: 3, 9: None, 10: 4, 11: None, 12: 5},  # 1H-pyrrolo[2,3-b]quinoline
    "c1cc2cc[nH]c2nn1": {0: 3, 1: 4, 2: None, 3: 5, 4: 6, 5: None, 6: None, 7: None, 8: None},  # 1H-pyrrolo[2,3-c]pyridazine
    "c1cc2cc[nH]c2cn1": {0: 5, 1: 4, 2: None, 3: 3, 4: 2, 5: None, 6: None, 7: 7, 8: None},  # 1H-pyrrolo[2,3-c]pyridine
    "c1cc2cnnnc2[nH]1": {0: 6, 1: 5, 2: None, 3: 4, 4: None, 5: None, 6: None, 7: None, 8: None},  # 1H-pyrrolo[2,3-d][1,2,3]triazine
    "c1cc2cnncc2[nH]1": {0: 2, 1: 3, 2: None, 3: 4, 4: None, 5: None, 6: 7, 7: None, 8: None},  # 1H-pyrrolo[2,3-d]pyridazine
    "c1ncc2cc[nH]c2n1": {0: 2, 1: None, 2: 4, 3: None, 4: 5, 5: 6, 6: None, 7: None, 8: None},  # 1H-pyrrolo[2,3-d]pyrimidine
    "c1cc2[nH]ncnc-2n1": {0: 6, 1: 7, 2: None, 3: None, 4: None, 5: 3, 6: None, 7: None, 8: None},  # 1H-pyrrolo[2,3-e][1,2,4]triazine
    "c1c[nH]c2ccnc-2n1": {0: 3, 1: 2, 2: None, 3: None, 4: 7, 5: 6, 6: None, 7: None, 8: None},  # 1H-pyrrolo[2,3-e]pyrazine
    "c1cnc2cc[nH]c2c1": {0: 6, 1: 5, 2: None, 3: None, 4: 3, 5: 2, 6: None, 7: None, 8: 7},  # 1H-pyrrolo[3,2-b]pyridine
    "c1ccc2nc3cc[nH]c3cc2c1": {0: 7, 1: 6, 2: 5, 3: None, 4: None, 5: None, 6: 3, 7: 2, 8: None, 9: None, 10: 9, 11: None, 12: 8},  # 1H-pyrrolo[3,2-b]quinoline
    "c1cc2nccc-2[nH]n1": {0: 3, 1: 4, 2: None, 3: None, 4: 6, 5: 7, 6: None, 7: None, 8: None},  # 1H-pyrrolo[3,2-c]pyridazine
    "c1cc2[nH]ccc2cn1": {0: 6, 1: 7, 2: None, 3: None, 4: 2, 5: 3, 6: None, 7: 4, 8: None},  # 1H-pyrrolo[3,2-c]pyridine
    "c1cc2[nH]nncc-2n1": {0: 6, 1: 7, 2: None, 3: None, 4: None, 5: None, 6: 4, 7: None, 8: None},  # 1H-pyrrolo[3,2-d][1,2,3]triazine
    "c1cc2[nH]cncc-2n1": {0: 6, 1: 7, 2: None, 3: None, 4: 2, 5: None, 6: 4, 7: None, 8: None},  # 1H-pyrrolo[3,2-d]pyrimidine
    "c1nnc2[nH]ccc2n1": {0: 3, 1: None, 2: None, 3: None, 4: None, 5: 6, 6: 5, 7: None, 8: None},  # 1H-pyrrolo[3,2-e][1,2,4]triazine
    "c1c[nH]c2cncc-2c1": {0: 3, 1: 2, 2: None, 3: None, 4: 7, 5: None, 6: 5, 7: None, 8: 4},  # 1H-pyrrolo[3,4-b]pyridine
    "c1cc2cncc-2[nH]n1": {0: 3, 1: 4, 2: None, 3: 5, 4: None, 5: 7, 6: None, 7: None, 8: None},  # 1H-pyrrolo[3,4-c]pyridazine
    "c1ncc2[nH]nncc1-2": {0: 5, 1: None, 2: 7, 3: None, 4: None, 5: None, 6: None, 7: 4, 8: None},  # 1H-pyrrolo[3,4-d][1,2,3]triazine
    "c1ncc2cncc-2[nH]1": {0: 2, 1: None, 2: 4, 3: None, 4: 5, 5: None, 6: 7, 7: None, 8: None},  # 1H-pyrrolo[3,4-d]pyrimidine
    "c1n[nH]c2cncc-2n1": {0: 3, 1: None, 2: None, 3: None, 4: 7, 5: None, 6: 5, 7: None, 8: None},  # 1H-pyrrolo[3,4-e][1,2,4]triazine
    "c1c[nH]c2cncc-2n1": {0: 3, 1: 2, 2: None, 3: None, 4: 7, 5: None, 6: 5, 7: None, 8: None},  # 1H-pyrrolo[3,4-e]pyrazine
    "c1cnc2n[nH]nc2c1": {0: 6, 1: 5, 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: 7},  # 2H-[1,2,3]triazolo[4,5-b]pyridine
    "c1cc2n[nH]nc2nn1": {0: 6, 1: 7, 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: None},  # 2H-[1,2,3]triazolo[4,5-c]pyridazine
    "c1cc2n[nH]nc2cn1": {0: 6, 1: 7, 2: None, 3: None, 4: None, 5: None, 6: None, 7: 4, 8: None},  # 2H-[1,2,3]triazolo[4,5-c]pyridine
    "c1nnnc2n[nH]nc12": {0: 7, 1: None, 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: None},  # 2H-[1,2,3]triazolo[4,5-d][1,2,3]triazine
    "c1ncc2n[nH]nc2n1": {0: 5, 1: None, 2: 7, 3: None, 4: None, 5: None, 6: None, 7: None, 8: None},  # 2H-[1,2,3]triazolo[4,5-d]pyrimidine
    "c1nnc2n[nH]nc2n1": {0: 6, 1: None, 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: None},  # 2H-[1,2,3]triazolo[4,5-e][1,2,4]triazine
    "c1cnc2n[nH]nc2n1": {0: 5, 1: 5, 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: None},  # 2H-[1,2,3]triazolo[4,5-e]pyrazine
    "c1cnc2nc[nH]c2c1": {0: 6, 1: 5, 2: None, 3: None, 4: None, 5: 2, 6: None, 7: None, 8: 7},  # 3H-imidazo[4,5-b]pyridine
    "c1ccc2nc3nc[nH]c3cc2c1": {0: 7, 1: 6, 2: 5, 3: None, 4: None, 5: None, 6: None, 7: 2, 8: None, 9: None, 10: 9, 11: None, 12: 8},  # 3H-imidazo[4,5-b]quinoline
    "c1ncc2[nH]cnc2n1": {0: 2, 1: None, 2: 6, 3: None, 4: None, 5: 8, 6: None, 7: None, 8: None},  # 7H-purine
    "c1ncc2nc[nH]c2n1": {0: 2, 1: None, 2: 6, 3: None, 4: None, 5: 8, 6: None, 7: None, 8: None},  # 9H-purine
    # Phase 618: non-aromatic and partially aromatic ring locants
    "c1ccc2c(c1)CCCC2": {0: 6, 1: 6, 2: 5, 3: None, 4: None, 5: 5, 6: 1, 7: 2, 8: 2, 9: 1},  # 1,2,3,4-tetrahydronaphthalene
    # Phase 628: 1,2-dihydronaphthalene (C1=C2 sp2; C3,C4 sp3 — note: IUPAC 1,2-dihydro means sp3 at 1,2)
    # Actually atoms 0,1 are the sp2 C=C and atoms 8,9 are the sp3 CH2 positions (C1 and C2)
    "C1=Cc2ccccc2CC1": {0: 3, 1: 4, 2: None, 3: 5, 4: 6, 5: 7, 6: 8, 7: None, 8: 1, 9: 2},
    "C1CCOOC1": {0: 4, 1: 4, 2: 3, 3: None, 4: None, 5: 3},  # 1,2-dioxane
    "C1COOC1": {0: 4, 1: 3, 2: None, 3: None, 4: 3},  # 1,2-dioxolane
    "C1CCSSC1": {0: 4, 1: 4, 2: 3, 3: None, 4: None, 5: 3},  # 1,2-dithiane
    "C1CSSC1": {0: 4, 1: 3, 2: None, 3: None, 4: 3},  # 1,2-dithiolane
    "C1COCOC1": {0: 5, 1: 4, 2: None, 3: 2, 4: None, 5: 4},  # 1,3-dioxane
    "C1COCO1": {0: 4, 1: 4, 2: None, 3: 2, 4: None},  # 1,3-dioxolane
    "C1CSCS1": {0: 4, 1: 4, 2: None, 3: 2, 4: None},  # 1,3-dithiolane
    "C1CSCO1": {0: 5, 1: 4, 2: None, 3: 2, 4: None},  # 1,3-oxathiolane
    "C1COCCO1": {0: 2, 1: 2, 2: None, 3: 2, 4: 2, 5: None},  # 1,4-dioxane
    "C1CSCCS1": {0: 2, 1: 2, 2: None, 3: 2, 4: 2, 5: None},  # 1,4-dithiane
    "C1CSCCO1": {0: 2, 1: 3, 2: None, 3: 3, 4: 2, 5: None},  # 1,4-oxathiane
    "c1ccc2c(c1)CCC2": {0: 5, 1: 5, 2: 4, 3: None, 4: None, 5: 4, 6: 1, 7: 2, 8: 1},  # indane
    "c1ccc2c(c1)CCCS2": {0: 6, 1: 7, 2: 8, 3: None, 4: None, 5: 5, 6: 4, 7: 3, 8: 2, 9: None},  # thiochroman
    # Phase 619: non-aromatic double-bond compounds
    "C1=Cc2ccccc2C1": {0: 2, 1: 3, 2: None, 3: 4, 4: 5, 5: 6, 6: 7, 7: None, 8: 1},  # 1H-indene
    "C1=Cc2ccccc2OC1": {0: 3, 1: 4, 2: None, 3: 5, 4: 6, 5: 7, 6: 8, 7: None, 8: None, 9: 2},  # 2H-chromene
    "C1=COc2ccccc2C1": {0: 3, 1: 2, 2: None, 3: None, 4: 8, 5: 7, 6: 6, 7: 5, 8: None, 9: 4},  # 4H-chromene
    "O=C1CCc2ccccc21": {0: None, 1: None, 2: 2, 3: 3, 4: None, 5: 4, 6: 5, 7: 6, 8: 7, 9: None},  # indan-1-one
    "O=C1CC(=O)c2ccccc21": {0: None, 1: None, 2: 2, 3: None, 4: None, 5: None, 6: 4, 7: 5, 8: 5, 9: 4, 10: None},  # indane-1,3-dione
    "O=C1OC(=O)c2ccccc21": {0: None, 1: None, 2: None, 3: None, 4: None, 5: None, 6: 4, 7: 5, 8: 5, 9: 4, 10: None},  # isobenzofuran-1,3-dione
    "O=C1C=CC(=O)O1": {0: None, 1: None, 2: 3, 3: 3, 4: None, 5: None, 6: None},  # furan-2,5-dione
    "C1=Nc2cccc3cccc(c23)N1": {0: 2, 1: None, 2: None, 3: 4, 4: 5, 5: 6, 6: None, 7: 6, 8: 5, 9: 4, 10: None, 11: None, 12: None},  # perimidine
    # Phase 620: non-aromatic N=C fused heterocycles
    "C1=NCc2ccncc21": {0: 3, 1: None, 2: 1, 3: None, 4: 7, 5: 6, 6: None, 7: 4, 8: None},  # 1H-pyrrolo[3,4-c]pyridine
    "C1=NC=C2CN=NC=C12": {0: 5, 1: None, 2: 7, 3: None, 4: 1, 5: None, 6: None, 7: 4, 8: None},  # 1H-pyrrolo[3,4-d]pyridazine
    "C1=NC2=CCN=NC2=N1": {0: 6, 1: None, 2: None, 3: 4, 4: 3, 5: None, 6: None, 7: None, 8: None},  # 3H-imidazo[4,5-c]pyridazine
    "c1ccc2ncccc2c1":   {0: 6, 1: 7, 2: 8, 3: None, 4: 1, 5: 2, 6: 3, 7: 4, 8: None, 9: 5},
    "c1ccc2cnccc2c1":   {0: 6, 1: 7, 2: 8, 3: None, 4: 1, 5: 2, 6: 3, 7: 4, 8: None, 9: 5},
    # Phase 552: quinoxaline, quinazoline, phthalazine, cinnoline — same 10-atom map
    "c1ccc2nccnc2c1":   {0: 6, 1: 7, 2: 8, 3: None, 4: 1, 5: 2, 6: 3, 7: 4, 8: None, 9: 5},
    "c1ccc2ncncc2c1":   {0: 6, 1: 7, 2: 8, 3: None, 4: 1, 5: 2, 6: 3, 7: 4, 8: None, 9: 5},
    "c1ccc2cnncc2c1":   {0: 6, 1: 7, 2: 8, 3: None, 4: 1, 5: 2, 6: 3, 7: 4, 8: None, 9: 5},
    "c1ccc2nnccc2c1":   {0: 6, 1: 7, 2: 8, 3: None, 4: 1, 5: 2, 6: 3, 7: 4, 8: None, 9: 5},
    "c1ccc2[nH]ccc2c1": {0: 5, 1: 6, 2: 7, 3: None, 4: 1, 5: 2, 6: 3, 7: None, 8: 4},
    "c1ccc2[nH]cnc2c1": {0: 5, 1: 6, 2: 7, 3: None, 4: 1, 5: 2, 6: 3, 7: None, 8: 4},
    "c1ccc2occc2c1":    {0: 5, 1: 6, 2: 7, 3: None, 4: 1, 5: 2, 6: 3, 7: None, 8: 4},
    "c1ccc2sccc2c1":    {0: 5, 1: 6, 2: 7, 3: None, 4: 1, 5: 2, 6: 3, 7: None, 8: 4},
    # Phase 553: 1H-indazole and 1H-benzotriazole — same 9-atom map as indole/benzimidazole
    "c1ccc2[nH]ncc2c1": {0: 5, 1: 6, 2: 7, 3: None, 4: 1, 5: 2, 6: 3, 7: None, 8: 4},
    "c1ccc2[nH]nnc2c1": {0: 5, 1: 6, 2: 7, 3: None, 4: 1, 5: 2, 6: 3, 7: None, 8: 4},
    # Phase 627: 2H-indazole (N1 at pos 1, N-H at pos 2; C3 at 3, C3a/C7a junctions)
    "c1ccc2n[nH]cc2c1": {0: 5, 1: 6, 2: 7, 3: None, 4: None, 5: None, 6: 3, 7: None, 8: 4},
    # Phase 554: 1,3-benzothiazole and 1,3-benzoxazole (S/O at 1, C at 2, N at 3)
    "c1ccc2scnc2c1":    {0: 5, 1: 6, 2: 7, 3: None, 4: None, 5: 2, 6: None, 7: None, 8: 4},
    "c1ccc2ocnc2c1":    {0: 5, 1: 6, 2: 7, 3: None, 4: None, 5: 2, 6: None, 7: None, 8: 4},
    # Phase 554: 1,2-benzisothiazole, 1,2-benzisoxazole, 2,1-benzisothiazole (C at 3, heteroatoms at 1,2)
    "c1ccc2sncc2c1":    {0: 5, 1: 6, 2: 7, 3: None, 4: None, 5: None, 6: 3, 7: None, 8: 4},
    "c1ccc2oncc2c1":    {0: 5, 1: 6, 2: 7, 3: None, 4: None, 5: None, 6: 3, 7: None, 8: 4},
    "c1ccc2nscc2c1":    {0: 5, 1: 6, 2: 7, 3: None, 4: None, 5: None, 6: 3, 7: None, 8: 4},
    # Phase 560: phenazine (14-atom tricyclic, D2h; only 2 unique positions: locant 1 and 2)
    # Ring1: 12(10a,jct)-13(C1)-0(C2)-1(C3)-2(C4)-3(4a,jct); Ring2 central: 3-4(N)-5(4b,jct)-10(9a,jct)-11(N)-12
    # Ring3: 5-6(C6?)-7(C7)-8(C8)-9(C9?)-10; D2h: {C1,C4,C6,C9} equiv, {C2,C3,C7,C8} equiv
    "c1ccc2nc3ccccc3nc2c1": {0: 2, 1: 3, 2: 4, 3: None, 4: None, 5: None, 6: 6, 7: 7, 8: 8, 9: 9, 10: None, 11: None, 12: None, 13: 1},
    # Phase 560: 1,10-phenanthroline (14-atom tricyclic, C2; 4 unique positions: 2,3,4,5)
    # Ring A: 0(C3)-1(C2)-2(N1)-3(10b,jct)-4(4a,jct)-5(C4); Ring B central: 3-4-6(C6)-7(C5)-8(10a,jct)-13(4b,jct)
    # Ring C: 8-9(C7)-10(C8)-11(C9)-12(N10)-13; C2 symmetry: 2↔9,3↔8,4↔7,5↔6
    "c1cnc2c(c1)ccc1cccnc12": {0: 3, 1: 2, 2: None, 3: None, 4: None, 5: 4, 6: 6, 7: 5, 8: None, 9: 7, 10: 8, 11: 9, 12: None, 13: None},
    # Phase 562: remaining phenanthroline isomers (2,7; 2,6; 1,6; 3,6)
    # 2,7-phen: ring containing N2 also holds C1,C9,C10; central ring holds C3,C4; far ring holds N7,C5,C6,C8
    "c1cnc2ccc3ccncc3c2c1": {0: 10, 1: 1, 2: None, 3: None, 4: 3, 5: 4, 6: None, 7: 5, 8: 6, 9: None, 10: 8, 11: None, 12: None, 13: 9},
    # 2,6-phen: N2 in ring with C1,C9,C10; N6 in central ring with C7; outer ring has C3,C4,C5,C8
    "c1ccc2c(c1)ncc1ccncc12": {0: 8, 1: 9, 2: 10, 3: None, 4: None, 5: 7, 6: None, 7: 5, 8: None, 9: 4, 10: 3, 11: None, 12: 1, 13: None},
    # 1,6-phen: ring A has N1,C2,C3,C4; central ring has C5,N6; ring C has C7,C8,C9,C10
    "c1cnc2c(c1)cnc1ccccc12": {0: 3, 1: 2, 2: None, 3: None, 4: None, 5: 4, 6: 5, 7: None, 8: None, 9: 7, 10: 8, 11: 9, 12: 10, 13: None},
    # 3,6-phen: N3,C1,C2,C4 in far ring; N6,C7 in central ring; C5,C8,C9,C10 in outer ring
    "c1ccc2c(c1)ncc1cnccc12": {0: 8, 1: 9, 2: 10, 3: None, 4: None, 5: 7, 6: None, 7: 5, 8: None, 9: 4, 10: None, 11: 2, 12: 1, 13: None},
    # 3,5-phen: peripheral path 1(a5)-2(a6)-N3(a7)-4b-4(a9)-N5(a10)-6(a11)-7(a12)-8a-10a-8(a2)-9(a1)-10(a0)
    "c1ccc2c(c1)cnc1cnccc12": {0: 10, 1: 9, 2: 8, 3: None, 4: None, 5: 1, 6: 2, 7: None, 8: None, 9: 4, 10: None, 11: 6, 12: 7, 13: None},
    # 4,5-phen: peripheral path 1(a0)-2(a5)-3(a6)-N4(a7)-4b-N5(a9)-6(a10)-7(a11)-8(a12)-8a-10a-9(a2)-10(a1)
    "c1ccc2c(c1)cnc1ncccc12": {0: 1, 1: 10, 2: 9, 3: None, 4: None, 5: 2, 6: 3, 7: None, 8: None, 9: None, 10: 6, 11: 7, 12: 8, 13: None},
    # Phase 561: 4,7-phenanthroline (C2 symmetric; unique positions: 1,2,3,5,6)
    "c1cnc2ccc3ncccc3c2c1": {0: 2, 1: 3, 2: None, 3: None, 4: 6, 5: 5, 6: None, 7: None, 8: 8, 9: 9, 10: 10, 11: None, 12: None, 13: 1},
    # Phase 561: 1,7-phenanthroline (no symmetry; 8 unique C positions: 2,3,4,5,6,8,9,10)
    "c1cnc2c(c1)ccc1ncccc12": {0: 3, 1: 2, 2: None, 3: None, 4: None, 5: 4, 6: 5, 7: 6, 8: None, 9: None, 10: 8, 11: 9, 12: 10, 13: None},
    # Phase 561: 1,8-phenanthroline (no symmetry; 8 unique C positions: 2,3,4,5,6,7,9,10)
    "c1cnc2c(c1)ccc1cnccc12": {0: 3, 1: 2, 2: None, 3: None, 4: None, 5: 4, 6: 5, 7: 6, 8: None, 9: 7, 10: None, 11: 9, 12: 10, 13: None},
    # Phase 559: acridine (14-atom tricyclic; C2v symmetry)
    # Ring A: 10a(3)-C1(2)-C2(1)-C3(0)-C4(13)-4a(12); Ring B central: 10a(3)-N10(4)-4b(5)-8a(10)-C9(11)-4a(12)
    # Ring C: 4b(5)-C8(6)-C7(7)-C6(8)-C5(9)-8a(10); GetSubstructMatch picks lower-locant atom for symmetric pairs
    "c1ccc2nc3ccccc3cc2c1": {0: 3, 1: 2, 2: 1, 3: None, 4: None, 5: None, 6: 8, 7: 7, 8: 6, 9: 5, 10: None, 11: 9, 12: None, 13: 4},
    # Phase 558: 2,x-naphthyridines (N2=atom9, C8a=atom7, C4a=atom2, C1=atom8, C3=atom0, C4=atom1)
    # Both 2,6 and 2,7 have C2 symmetry; only lower-locant half included
    "c1cc2cnccc2cn1": {0: 3, 1: 4, 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: 1, 9: None},  # 2,6
    "c1cc2ccncc2cn1": {0: 3, 1: 4, 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: 1, 9: None},  # 2,7
    # Phase 557: 1,x-naphthyridines (10-atom bicyclic diazines; 8a junction adj to N1, 4a junction opposite)
    # For all: N1=atom2(None), 8a=atom3(None), C2=atom1, C3=atom0, C4=atom9, 4a=atom8(None)
    # 1,5 and 1,8 have C2 symmetry (C2≡C7 etc.); only include lower-locant half to avoid duplicate SMILES
    "c1cnc2cccnc2c1": {0: 3, 1: 2, 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: None, 9: 4},  # 1,5
    "c1cnc2ccncc2c1": {0: 3, 1: 2, 2: None, 3: None, 4: 8, 5: 7, 6: None, 7: 5, 8: None, 9: 4},  # 1,6
    "c1cnc2cnccc2c1": {0: 3, 1: 2, 2: None, 3: None, 4: 8, 5: None, 6: 6, 7: 5, 8: None, 9: 4},  # 1,7
    "c1cnc2ncccc2c1": {0: 3, 1: 2, 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: None, 9: 4},  # 1,8
    # Phase 556: 1,3-benzodioxole (C5 substitutable; symmetry at C5/C7 but canonical gives lowest locant)
    "c1ccc2c(c1)OCO2":  {0: 5, 1: 6, 2: 7, 3: None, 4: None, 5: 4, 6: None, 7: None, 8: None},
    # Phase 555: benzo-fused thia/oxa-diazoles and selenadiazoles (only benzo C4-C7 substitutable)
    "c1ccc2snnc2c1":    {0: 5, 1: 6, 2: 7, 3: None, 4: None, 5: None, 6: None, 7: None, 8: 4},
    "c1ccc2nsnc2c1":    {0: 5, 1: 6, 2: 7, 3: None, 4: None, 5: None, 6: None, 7: None, 8: 4},
    "c1ccc2nonc2c1":    {0: 5, 1: 6, 2: 7, 3: None, 4: None, 5: None, 6: None, 7: None, 8: 4},
    "c1ccc2n[se]nc2c1": {0: 5, 1: 6, 2: 7, 3: None, 4: None, 5: None, 6: None, 7: None, 8: 4},
    "c1ccc2[se]nnc2c1": {0: 5, 1: 6, 2: 7, 3: None, 4: None, 5: None, 6: None, 7: None, 8: 4},
    # Phase 158: 9H-carbazole (N at atom 6 = IUPAC position 9)
    "c1ccc2c(c1)[nH]c1ccccc12": {
        0: 2, 1: 3, 2: 4, 3: None, 4: None, 5: 1,
        6: 9, 7: None, 8: 8, 9: 7, 10: 6, 11: 5, 12: None,
    },
    # Phase 563: dibenzofuran/dibenzothiophene (6+5+6; O/S→None; positions 1-4,6-9; pairs sum to 10)
    # peripheral path: a2(1)-a1(2)-a0(3)-jct3-jct12-a11(9)-a10(8)-a9(7)-a8(6)-jct7-O/S(None)-jct4-a5(4)
    "c1ccc2c(c1)oc1ccccc12": {0: 3, 1: 2, 2: 1, 3: None, 4: None, 5: 4, 6: None, 7: None, 8: 6, 9: 7, 10: 8, 11: 9, 12: None},
    "c1ccc2c(c1)sc1ccccc12": {0: 3, 1: 2, 2: 1, 3: None, 4: None, 5: 4, 6: None, 7: None, 8: 6, 9: 7, 10: 8, 11: 9, 12: None},
    # Phase 563: fluorene (6+5+6; CH2→9; positions 1-8; pairs sum to 9)
    # peripheral path: a5(1)-a0(2)-a1(3)-a2(4)-jct3-jct12-a11(5)-a10(6)-a9(7)-a8(8)-jct7-CH2(a6,9)-jct4
    "c1ccc2c(c1)Cc1ccccc1-2": {0: 2, 1: 3, 2: 4, 3: None, 4: None, 5: 1, 6: 9, 7: None, 8: 8, 9: 7, 10: 6, 11: 5, 12: None},
    # Phase 563: xanthene/thioxanthene (6+6+6; CH2→9, O/S→None; positions 1-9; pairs sum to 9)
    # peripheral path: a5(1)-a0(2)-a1(3)-a2(4)-jct3-O/S(a13,None)-jct12-a11(5)-a10(6)-a9(7)-a8(8)-jct7-CH2(a6,9)-jct4
    "c1ccc2c(c1)Cc1ccccc1O2": {0: 2, 1: 3, 2: 4, 3: None, 4: None, 5: 1, 6: 9, 7: None, 8: 8, 9: 7, 10: 6, 11: 5, 12: None, 13: None},
    "c1ccc2c(c1)Cc1ccccc1S2": {0: 2, 1: 3, 2: 4, 3: None, 4: None, 5: 1, 6: 9, 7: None, 8: 8, 9: 7, 10: 6, 11: 5, 12: None, 13: None},
    # Phase 563: phenoxazine/phenothiazine (6+6+6; N→10, O/S→None; positions 1-4,6-9; pairs sum to 10)
    # peripheral path: a5(1)-a0(2)-a1(3)-a2(4)-jct3-O/S(a13,None)-jct12-a11(6)-a10(7)-a9(8)-a8(9)-jct7-N(a6,10)-jct4
    "c1ccc2c(c1)Nc1ccccc1O2": {0: 2, 1: 3, 2: 4, 3: None, 4: None, 5: 1, 6: 10, 7: None, 8: 9, 9: 8, 10: 7, 11: 6, 12: None, 13: None},
    "c1ccc2c(c1)Nc1ccccc1S2": {0: 2, 1: 3, 2: 4, 3: None, 4: None, 5: 1, 6: 10, 7: None, 8: 9, 9: 8, 10: 7, 11: 6, 12: None, 13: None},
    # Phase 563: phenoxathiin (6+6+6; C2-symmetric via O/S axis; O+S→None; positions 1-4,6-9; pairs sum to 10)
    # peripheral path: a2(1)-a1(2)-a0(3)-jct3-S(a13,None)-jct12-a11(9)-a10(8)-a9(7)-a8(6)-jct7-O(a6,None)-jct4-a5(4)
    "c1ccc2c(c1)Oc1ccccc1S2": {0: 3, 1: 2, 2: 1, 3: None, 4: None, 5: 4, 6: None, 7: None, 8: 6, 9: 7, 10: 8, 11: 9, 12: None, 13: None},
    # Phase 563: thianthrene (6+6+6; D2h symmetry: 4-fold degeneracy → only 2 unique C environments)
    # Group A (adj to any S junction): {a2,a5,a8,a11} → all map to locant 1 (lowest)
    # Group B (non-junction-adj): {a0,a1,a9,a10} → all map to locant 2 (lowest)
    "c1ccc2c(c1)Sc1ccccc1S2": {0: 2, 1: 2, 2: 1, 3: None, 4: None, 5: 1, 6: None, 7: None, 8: 1, 9: 2, 10: 2, 11: 1, 12: None, 13: None},
    # Phase 626: 2,3-dihydrobenzofuran: O(1)-C(2)-C(3)-C(3a,junc)-C(4)-C(5)-C(6)-C(7)-C(7a,junc)
    "c1ccc2c(c1)CCO2": {0: 5, 1: 6, 2: 7, 3: None, 4: None, 5: 4, 6: 3, 7: 2, 8: None},
    # 1,3-dihydro-2-benzofuran (symmetric): C(1)-O(2)-C(3) symmetric; C(4)-C(5)-C(6)-C(7) benzene
    "c1ccc2c(c1)COC2": {0: 5, 1: 5, 2: 4, 3: None, 4: None, 5: 4, 6: 1, 7: None, 8: 1},
    # 2,3-dihydrobenzothiophene: S(1)-C(2)-C(3)-C(3a,junc)-C(4)-C(5)-C(6)-C(7)-C(7a,junc)
    "c1ccc2c(c1)CCS2": {0: 5, 1: 6, 2: 7, 3: None, 4: None, 5: 4, 6: 3, 7: 2, 8: None},
    # 1,3-dihydro-2-benzothiophene (symmetric): C(1)-S(2)-C(3) symmetric; C(4)-C(7) benzene
    "c1ccc2c(c1)CSC2": {0: 5, 1: 5, 2: 4, 3: None, 4: None, 5: 4, 6: 1, 7: None, 8: 1},
    # Phase 404: 部分飽和縮合環 (N/O heteroatom, fused with benzene)
    # indoline (2,3-dihydro-1H-indole): N(1)-C(2)-C(3)-C(3a,junc)-C(4)-C(5)-C(6)-C(7)-C(7a,junc)
    "c1ccc2c(c1)CCN2": {0: 5, 1: 6, 2: 7, 3: None, 4: None, 5: 4, 6: 3, 7: 2, 8: 1},
    # Phase 635: 2,3-dihydrobenzo[d] heterocycles (all share the same atom layout)
    # benzo[d]oxazole: O(1,None)-C(2,sp3)-N(3,NH)-C(3a,junc)-C(4)-C(5)-C(6)-C(7)-C(7a,junc)
    "c1ccc2c(c1)NCO2":  {0: 5, 1: 6, 2: 7, 3: None, 4: None, 5: 4, 6: 3, 7: 2, 8: None},
    "c1ccc2c(c1)NCS2":  {0: 5, 1: 6, 2: 7, 3: None, 4: None, 5: 4, 6: 3, 7: 2, 8: None},
    # benzo[d]isoxazole: O(1,None)-N(2,NH)-C(3,sp3)-C(3a,junc)-C(4)-C(5)-C(6)-C(7)-C(7a,junc)
    "c1ccc2c(c1)CNO2":  {0: 5, 1: 6, 2: 7, 3: None, 4: None, 5: 4, 6: 3, 7: 2, 8: None},
    "c1ccc2c(c1)CNS2":  {0: 5, 1: 6, 2: 7, 3: None, 4: None, 5: 4, 6: 3, 7: 2, 8: None},
    # 2,3-dihydro-1H-benzimidazole: C2-symmetric (N1≡N3→1, C4≡C7→4, C5≡C6→5)
    "c1ccc2c(c1)NCN2":  {0: 5, 1: 5, 2: 4, 3: None, 4: None, 5: 4, 6: 1, 7: 2, 8: 1},
    # 2,3-dihydro-1H-indazole: N(1)-N(2)-C(3,sp3)-C(3a,junc)-C(4)-C(5)-C(6)-C(7)-C(7a,junc)
    "c1ccc2c(c1)CNN2":  {0: 5, 1: 6, 2: 7, 3: None, 4: None, 5: 4, 6: 3, 7: 2, 8: 1},
    # 1,2,3,4-tetrahydroquinoline: N(1)-C(2)-C(3)-C(4)-C(4a,junc)-C(5)-C(6)-C(7)-C(8)-C(8a,junc)
    "c1ccc2c(c1)CCCN2": {0: 6, 1: 7, 2: 8, 3: None, 4: None, 5: 5, 6: 4, 7: 3, 8: 2, 9: 1},
    # Phase 627: 1,2,3,4-tetrahydroisoquinoline: C(1)-N(2)-C(3)-C(4)-C(4a,junc)-C(5)-C(6)-C(7)-C(8)-C(8a,junc)
    "c1ccc2c(c1)CCNC2": {0: 6, 1: 7, 2: 8, 3: None, 4: None, 5: 5, 6: 4, 7: 3, 8: None, 9: 1},
    # Phase 627: 1,2,3,4-tetrahydroquinoxaline (C2-symmetric: C5↔C8→5, C6↔C7→6, C2↔C3→2)
    "c1ccc2c(c1)NCCN2": {0: 6, 1: 6, 2: 5, 3: None, 4: None, 5: 5, 6: None, 7: 2, 8: 2, 9: None},
    # Phase 628: 2,3-dihydro-1,4-benzodioxine (C2-symmetric: C5↔C8→5, C6↔C7→6, C2↔C3→2)
    "c1ccc2c(c1)OCCO2": {0: 6, 1: 6, 2: 5, 3: None, 4: None, 5: 5, 6: None, 7: 2, 8: 2, 9: None},
    # Phase 628: 3,4-dihydro-2H-1,4-benzoxazine (O at 1, N at 4; C2-C3 aliphatic; C5-C8 benzo)
    "c1ccc2c(c1)NCCO2": {0: 6, 1: 7, 2: 8, 3: None, 4: None, 5: 5, 6: None, 7: 3, 8: 2, 9: None},
    # Phase 628: 3,4-dihydro-2H-1,4-benzothiazine (S at 1, N at 4; C2-C3 aliphatic; C5-C8 benzo)
    "c1ccc2c(c1)NCCS2": {0: 6, 1: 7, 2: 8, 3: None, 4: None, 5: 5, 6: None, 7: 3, 8: 2, 9: None},
    # Phase 628: 3,4-dihydroisoquinoline (C1=N2; C3-C4 aliphatic; C5-C8 benzo)
    "C1=NCCc2ccccc21": {0: 1, 1: None, 2: 3, 3: 4, 4: None, 5: 5, 6: 6, 7: 7, 8: 8, 9: None},
    # Phase 628: 3,4-dihydroquinoline (N1=C2; C3-C4 aliphatic; C5-C8 benzo)
    "C1=Nc2ccccc2CC1": {0: 2, 1: None, 2: None, 3: 8, 4: 7, 5: 6, 6: 5, 7: None, 8: 4, 9: 3},
    # Phase 628: 3,4-dihydroquinoxaline (N1=C2; N4; C3 aliphatic; C5-C8 benzo)
    "C1=Nc2ccccc2NC1": {0: 2, 1: None, 2: None, 3: 8, 4: 7, 5: 6, 6: 5, 7: None, 8: None, 9: 3},
    # Phase 629: 1,2-dihydroquinoline (N1; C2 sp3; C3=C4 sp2; C4a-C8a benzo junctions)
    "C1=Cc2ccccc2NC1": {0: 3, 1: 4, 2: None, 3: 5, 4: 6, 5: 7, 6: 8, 7: None, 8: None, 9: 2},
    # Phase 629: 3,4-dihydroquinoxalin-2(1H)-one (C2=O; C3 sp3; N4; C5-C8 benzo)
    "O=C1CNc2ccccc2N1": {0: None, 1: None, 2: 3, 3: None, 4: None, 5: 5, 6: 6, 7: 7, 8: 8, 9: None, 10: None},
    # Phase 630: benzo-fused 7-membered saturated/partially unsaturated rings
    # 6,7,8,9-tetrahydro-5H-benzo[7]annulene (C2-symmetric: 1≡4→1, 2≡3→2, 5≡9→5, 6≡8→6)
    "c1ccc2c(c1)CCCCC2": {0: 2, 1: 2, 2: 1, 3: None, 4: None, 5: 1, 6: 5, 7: 6, 8: 7, 9: 6, 10: 5},
    # 2,3,4,5-tetrahydro-1-benzazepine: N1-C2-C3-C4-C5-C5a(junc)-C6-C7-C8-C9-C9a(junc)
    "c1ccc2c(c1)CCCCN2": {0: 7, 1: 8, 2: 9, 3: None, 4: None, 5: 6, 6: 5, 7: 4, 8: 3, 9: 2, 10: 1},
    # 2,3,4,5-tetrahydro-1H-2-benzazepine: C1-N2-C3-C4-C5-C5a(junc)-C6-C7-C8-C9-C9a(junc)
    "c1ccc2c(c1)CCCNC2": {0: 7, 1: 8, 2: 9, 3: None, 4: None, 5: 6, 6: 5, 7: 4, 8: 3, 9: 2, 10: 1},
    # 2,3,4,5-tetrahydro-1H-3-benzazepine (C2-symmetric: 1≡5→1, 2≡4→2, 6≡9→6, 7≡8→7)
    "c1ccc2c(c1)CCNCC2": {0: 7, 1: 7, 2: 6, 3: None, 4: None, 5: 6, 6: 1, 7: 2, 8: 3, 9: 2, 10: 1},
    # 2,3,4,5-tetrahydro-1-benzoxepine: O1-C2-C3-C4-C5-C5a(junc)-C6-C7-C8-C9-C9a(junc)
    "c1ccc2c(c1)CCCCO2": {0: 7, 1: 8, 2: 9, 3: None, 4: None, 5: 6, 6: 5, 7: 4, 8: 3, 9: 2, 10: None},
    # 2,3,4,5-tetrahydro-1H-2-benzoxepine: C1-O2-C3-C4-C5-C5a(junc)-C6-C7-C8-C9-C9a(junc)
    "c1ccc2c(c1)CCCOC2": {0: 7, 1: 8, 2: 9, 3: None, 4: None, 5: 6, 6: 5, 7: 4, 8: 3, 9: None, 10: 1},
    # 2,3,4,5-tetrahydro-1H-3-benzoxepine (C2-symmetric: 1≡5→1, 2≡4→2, 6≡9→6, 7≡8→7)
    "c1ccc2c(c1)CCOCC2": {0: 7, 1: 7, 2: 6, 3: None, 4: None, 5: 6, 6: 1, 7: 2, 8: None, 9: 2, 10: 1},
    # 2,3,4,5-tetrahydro-1-benzothiepine: S1-C2-C3-C4-C5-C5a(junc)-C6-C7-C8-C9-C9a(junc)
    "c1ccc2c(c1)CCCCS2": {0: 7, 1: 8, 2: 9, 3: None, 4: None, 5: 6, 6: 5, 7: 4, 8: 3, 9: 2, 10: None},
    # 2,3,4,5-tetrahydro-1H-1,5-benzodiazepine (C2-symmetric: N1≡N5→1, C2≡C4→2, C6≡C9→6, C7≡C8→7)
    "c1ccc2c(c1)NCCCN2": {0: 7, 1: 7, 2: 6, 3: None, 4: None, 5: 6, 6: 1, 7: 2, 8: 3, 9: 2, 10: 1},
    # 2,3-dihydro-1H-1,5-benzodiazepine: N1-C2-C3-C4=N5-C5a(junc)-C6-C7-C8-C9-C9a(junc)
    "C1=Nc2ccccc2NCC1": {0: 4, 1: None, 2: None, 3: 6, 4: 7, 5: 8, 6: 9, 7: None, 8: 1, 9: 2, 10: 3},
    # 3,4-dihydro-1H-1-benzazepin-2(5H)-one: N1-C2(=O)-C3-C4-C5-C5a(junc)-C6-C7-C8-C9-C9a(junc)
    "O=C1CCCc2ccccc2N1": {0: None, 1: None, 2: 3, 3: 4, 4: 5, 5: None, 6: 6, 7: 7, 8: 8, 9: 9, 10: None, 11: 1},
    # 2,3-dihydro-1H-1,4-benzodiazepin-5(4H)-one: N1-C2-C3-N4-C5(=O)-C5a(junc)-C6-C7-C8-C9-C9a(junc)
    "O=C1NCCNc2ccccc21": {0: None, 1: None, 2: 4, 3: 3, 4: 2, 5: 1, 6: None, 7: 9, 8: 8, 9: 7, 10: 6, 11: None},
    # Phase 631: tricyclic partially saturated ring systems
    # 1,2,3,4-tetrahydrophenanthrene: C1-C2-C3-C4-C4a(junc)-C4b(junc)-C5-C6-C7-C8-C8a(junc)-C9-C10-C10a(junc)
    "c1ccc2c3c(ccc2c1)CCCC3": {0: 7, 1: 6, 2: 5, 3: None, 4: None, 5: None, 6: 10, 7: 9, 8: None, 9: 8, 10: 1, 11: 2, 12: 3, 13: 4},
    # 9,10-dihydrophenanthrene (C2-symmetric: 1≡6, 2≡5, 3≡4, 9≡10)
    "c1ccc2c(c1)CCc1ccccc1-2": {0: 2, 1: 3, 2: 4, 3: None, 4: None, 5: 1, 6: 9, 7: 9, 8: None, 9: 1, 10: 2, 11: 3, 12: 4, 13: None},
    # 1,2-dihydrophenanthrene: C1-C2 sp3; C3=C4 sp2; benzo at C4a-C10a
    "C1=Cc2c(ccc3ccccc23)CC1": {0: 3, 1: 4, 2: None, 3: None, 4: 10, 5: 9, 6: None, 7: 8, 8: 7, 9: 6, 10: 5, 11: None, 12: 1, 13: 2},
    # 1,2,3,4-tetrahydroanthracene (C2-symmetric: 1≡4, 2≡3, 5≡9, 6≡... wait use OPSIN values)
    "c1ccc2cc3c(cc2c1)CCCC3": {0: 6, 1: 6, 2: 5, 3: None, 4: 9, 5: None, 6: None, 7: 9, 8: None, 9: 5, 10: 1, 11: 2, 12: 2, 13: 1},
    # 9,10-dihydroanthracene (C2-symmetric: 1≡8, 2≡7, 3≡6, 4≡5, 9≡10)
    "c1ccc2c(c1)Cc1ccccc1C2": {0: 2, 1: 2, 2: 1, 3: None, 4: None, 5: 1, 6: 9, 7: None, 8: 1, 9: 2, 10: 2, 11: 1, 12: None, 13: 9},
    # 1,2-dihydroanthracene: C1-C2 sp3; C3=C4 sp2; benzo rings at C4a-C9a
    "C1=Cc2cc3ccccc3cc2CC1": {0: 3, 1: 4, 2: None, 3: 10, 4: None, 5: 5, 6: 6, 7: 7, 8: 8, 9: None, 10: 9, 11: None, 12: 1, 13: 2},
    # 2,3-dihydro-1H-phenalene (C2-symmetric via C3v: 1≡3, 4≡9, 5≡8, 6≡7→6)
    "c1cc2c3c(cccc3c1)CCC2": {0: 5, 1: 4, 2: None, 3: None, 4: None, 5: 4, 6: 5, 7: 6, 8: None, 9: 6, 10: 1, 11: 2, 12: 1},
    # Phase 632: 5,6,7,8-tetrahydro N-heterocycles (aromatic ring intact, benzo ring saturated)
    # 5,6,7,8-tetrahydroquinoline: N(1)-C(2)-C(3)-C(4)-C(4a,junc)-C(5)-C(6)-C(7)-C(8)-C(8a,junc)
    "c1cnc2c(c1)CCCC2": {0: 3, 1: 2, 2: None, 3: None, 4: None, 5: 4, 6: 5, 7: 6, 8: 7, 9: 8},
    # 5,6,7,8-tetrahydroisoquinoline: C(1)-N(2)-C(3)-C(4)-C(4a,junc)-C(5)-C(6)-C(7)-C(8)-C(8a,junc)
    "c1cc2c(cn1)CCCC2": {0: 3, 1: 4, 2: None, 3: None, 4: 1, 5: None, 6: 8, 7: 7, 8: 6, 9: 5},
    # 5,6,7,8-tetrahydroquinoxaline (C2-symmetric: 2≡3, 5≡8, 6≡7)
    "c1cnc2c(n1)CCCC2": {0: 2, 1: 2, 2: None, 3: None, 4: None, 5: None, 6: 5, 7: 6, 8: 6, 9: 5},
    # 5,6,7,8-tetrahydroquinazoline: N(1)-C(2)-N(3)-C(4)-C(4a,junc)-C(5)-C(6)-C(7)-C(8)-C(8a,junc)
    "c1ncc2c(n1)CCCC2": {0: 2, 1: None, 2: 4, 3: None, 4: None, 5: None, 6: 8, 7: 7, 8: 6, 9: 5},
    # 5,6,7,8-tetrahydrocinnoline: N(1)-N(2)-C(3)-C(4)-C(4a,junc)-C(5)-C(6)-C(7)-C(8)-C(8a,junc)
    "c1cc2c(nn1)CCCC2": {0: 3, 1: 4, 2: None, 3: None, 4: None, 5: None, 6: 8, 7: 7, 8: 6, 9: 5},
    # 5,6,7,8-tetrahydrophthalazine (C2-symmetric: 1≡4, 5≡8, 6≡7)
    "c1nncc2c1CCCC2": {0: 1, 1: None, 2: None, 3: 1, 4: None, 5: None, 6: 5, 7: 6, 8: 6, 9: 5},
    # 1,2,3,4-tetrahydroquinazoline: N(1)-C(2)-N(3)-C(4)-C(4a,junc)-C(5)-C(6)-C(7)-C(8)-C(8a,junc)
    "c1ccc2c(c1)CNCN2": {0: 6, 1: 7, 2: 8, 3: None, 4: None, 5: 5, 6: 4, 7: 3, 8: 2, 9: 1},
    # quinolizine: C(1)-C(2)-C(3)-C(4)-N(4a,junc)-C(5)-C(6)-C(7)-C(8)-C(9)-C(9a,junc)
    "C1=CCN2C=CC=CC2=C1": {0: 2, 1: 3, 2: 4, 3: None, 4: 6, 5: 7, 6: 8, 7: 9, 8: None, 9: 1},
    # Phase 636: tetrahydronaphthyridines
    # 1,2,3,4-tetrahydro-1,5-naphthyridine: N(5,None)-C(6)-C(7)-C(8)-C(8a,junc)-C(4a,junc)-C(4)-C(3)-C(2)-N(1,H)
    "c1cnc2c(c1)NCCC2": {0: 7, 1: 6, 2: None, 3: None, 4: None, 5: 8, 6: 1, 7: 2, 8: 3, 9: 4},
    # 5,6,7,8-tetrahydro-1,6-naphthyridine: N(1,None)-C(2)-C(3)-C(4)-C(4a,junc)-C(8a,junc)-C(5)-N(6,H)-C(7)-C(8)
    "c1cnc2c(c1)CNCC2": {0: 3, 1: 2, 2: None, 3: None, 4: None, 5: 4, 6: 5, 7: 6, 8: 7, 9: 8},
    # 5,6,7,8-tetrahydro-1,7-naphthyridine: N(1,None)-C(2)-C(3)-C(4)-C(4a,junc)-C(8a,junc)-C(5)-C(6)-N(7,H)-C(8)
    "c1cnc2c(c1)CCNC2": {0: 3, 1: 2, 2: None, 3: None, 4: None, 5: 4, 6: 5, 7: 6, 8: 7, 9: 8},
    # 1,2,3,4-tetrahydro-1,8-naphthyridine: N(8,None)-C(7)-C(6)-C(5)-C(4a,junc)-C(8a,junc)-C(4)-C(3)-C(2)-N(1,H)
    "c1cnc2c(c1)CCCN2": {0: 6, 1: 7, 2: None, 3: None, 4: None, 5: 5, 6: 4, 7: 3, 8: 2, 9: 1},
    # 1,2,3,4-tetrahydro-1,6-naphthyridine: N(6,None)-C(7)-C(8)-C(8a,junc)-C(5)-C(4a,junc)-C(4)-C(3)-C(2)-N(1,H)
    "c1cc2c(cn1)CCCN2": {0: 7, 1: 8, 2: None, 3: None, 4: 5, 5: None, 6: 4, 7: 3, 8: 2, 9: 1},
    # 1,2,3,4-tetrahydro-1,7-naphthyridine: atom5=N(7,aromatic,None), atom6=N(1,sp3,H)
    "c1cc2c(cn1)NCCC2": {0: 6, 1: 5, 2: None, 3: None, 4: 8, 5: None, 6: 1, 7: 2, 8: 3, 9: 4},
    # Phase 637: 1,2,3,4-tetrahydro-2,6-naphthyridine: N(6,None)-C(7)-C(8)-C(8a,junc)-C(5)-C(4a,junc)-C(4)-C(3)-N(2,H)-C(1)
    "c1cc2c(cn1)CCNC2": {0: 7, 1: 8, 2: None, 3: None, 4: 5, 5: None, 6: 4, 7: 3, 8: 2, 9: 1},
    # 1,2,3,4-tetrahydro-2,7-naphthyridine: atom5=N(7,aromatic,None), atom6=C(1), atom7=N(2,H)
    "c1cc2c(cn1)CNCC2": {0: 6, 1: 5, 2: None, 3: None, 4: 8, 5: None, 6: 1, 7: 2, 8: 3, 9: 4},
    # Phase 638: tricyclic acridine/phenanthridine derivatives (14 atoms)
    # 1,2,3,4-tetrahydroacridine: sp3 ring C1-C4, aromatic ring with N10, C9
    "c1ccc2nc3c(cc2c1)CCCC3": {0: 7, 1: 6, 2: 5, 3: None, 4: None, 5: None, 6: None, 7: 9, 8: None, 9: 8, 10: 1, 11: 2, 12: 3, 13: 4},
    # 1,2,3,4-tetrahydrophenanthridine: sp3 ring C1-C4, aromatic ring with N5, C6
    "c1ccc2c3c(ncc2c1)CCCC3": {0: 8, 1: 9, 2: 10, 3: None, 4: None, 5: None, 6: None, 7: 6, 8: None, 9: 7, 10: 4, 11: 3, 12: 2, 13: 1},
    # 9,10-dihydroacridine: C2-symmetric (1≡5→1, 2≡6→2, 3≡7→3, 4≡8→4); C9(sp3,CH2), N10(sp3,NH)
    "c1ccc2c(c1)Cc1ccccc1N2": {0: 2, 1: 3, 2: 4, 3: None, 4: None, 5: 1, 6: 9, 7: None, 8: 1, 9: 2, 10: 3, 11: 4, 12: None, 13: 10},
    # 1,2,3,4,5,6,7,8-octahydroacridine: C2-symmetric (1≡5→1, 2≡6→2, 3≡7→3, 4≡8→4); C9(ar,CH)
    "c1c2c(nc3c1CCCC3)CCCC2": {0: 9, 1: None, 2: None, 3: None, 4: None, 5: None, 6: 1, 7: 2, 8: 3, 9: 4, 10: 4, 11: 3, 12: 2, 13: 1},
    # Phase 639: 1,2,3,4-tetrahydrobenzo[g/f/h]quinoline (14 atoms; N in terminal ring)
    # benzo[g]: N(1,sp3,NH)-C2-C3-C4-C(4a,junc)-C5-C(5a,junc)-C6-C7-C8-C9-C(9a,junc)-C10-C(10a,junc)
    "c1ccc2cc3c(cc2c1)CCCN3": {0: 7, 1: 8, 2: 9, 3: None, 4: 10, 5: None, 6: None, 7: 5, 8: None, 9: 6, 10: 4, 11: 3, 12: 2, 13: 1},
    # benzo[f]: C1-C2-C3-N(4,sp3,NH)-C(4a,junc)-C5-C6-C7-C8-C(8a,junc)-C9-C10-C(10a,junc)-C(10b,junc)
    "c1ccc2c3c(ccc2c1)NCCC3": {0: 8, 1: 9, 2: 10, 3: None, 4: None, 5: None, 6: 5, 7: 6, 8: None, 9: 7, 10: 4, 11: 3, 12: 2, 13: 1},
    # benzo[h]: N(1,sp3,NH)-C2-C3-C4-C(4a,junc)-C5-C6-C7-C8-C(8a,junc)-C9-C10-C(10a,junc)-C(10b,junc)
    "c1ccc2c3c(ccc2c1)CCCN3": {0: 8, 1: 9, 2: 10, 3: None, 4: None, 5: None, 6: 5, 7: 6, 8: None, 9: 7, 10: 4, 11: 3, 12: 2, 13: 1},
    # Phase 640: 9,10-dihydrophenanthridine and 5,6-dihydrophenanthridine (14 atoms)
    # 9,10-dihydro: N5(arom,None)-C(4a,junc,None)-C(4b,junc,None)-C(8a,junc,None)-C(10a,junc,None); C9,C10 sp3
    "C1=Cc2cnc3ccccc3c2CC1": {0: 8, 1: 7, 2: None, 3: 6, 4: None, 5: None, 6: 4, 7: 3, 8: 2, 9: 1, 10: None, 11: None, 12: 10, 13: 9},
    # 5,6-dihydro: N5(sp3,NH)->5, C6(sp3,CH2)->6; C(4a,junc)-C(4b,junc)-C(8a,junc)-C(10a,junc) None
    "c1ccc2c(c1)CNc1ccccc1-2": {0: 8, 1: 9, 2: 10, 3: None, 4: None, 5: 7, 6: 6, 7: 5, 8: None, 9: 4, 10: 3, 11: 2, 12: 1, 13: None},
    # Phase 641: 1,2,3,4-tetrahydrobenzo[f/g/h]isoquinoline (14 atoms; N in terminal ring)
    # benzo[f]: C1-C2-N(3,sp3,NH)-C4-C(4a,junc)-C5-C6-C7-C8-C(8a,junc)-C9-C10-C(10a,junc)-C(10b,junc)
    "c1ccc2c3c(ccc2c1)CNCC3": {0: 8, 1: 9, 2: 10, 3: None, 4: None, 5: None, 6: 5, 7: 6, 8: None, 9: 7, 10: 4, 11: 3, 12: 2, 13: 1},
    # benzo[g]: C1-C2-C3-C4-C(4a,junc)-C5-C(5a,junc)-C6-C7-C8-C9-C(9a,junc)-C10-C(10a,junc); N at 2
    "c1ccc2cc3c(cc2c1)CCNC3": {0: 7, 1: 8, 2: 9, 3: None, 4: 10, 5: None, 6: None, 7: 5, 8: None, 9: 6, 10: 4, 11: 3, 12: 2, 13: 1},
    # benzo[h]: C1-C2-C3-C4-C(4a,junc)-C5-C6-C7-C8-C(8a,junc)-C9-C10-C(10a,junc)-C(10b,junc); N at 3
    "c1ccc2c3c(ccc2c1)CCNC3": {0: 8, 1: 9, 2: 10, 3: None, 4: None, 5: None, 6: 5, 7: 6, 8: None, 9: 7, 10: 4, 11: 3, 12: 2, 13: 1},
    # Phase 642: 1,2-dihydroacridine (14 atoms; N10 aromatic, C1/C2 sp3, C3=C4 double bond)
    "C1=Cc2nc3ccccc3cc2CC1": {0: 3, 1: 4, 2: None, 3: None, 4: None, 5: 5, 6: 6, 7: 7, 8: 8, 9: None, 10: 9, 11: None, 12: 1, 13: 2},
    # Phase 643: 3,4-dihydroacridine (14 atoms; N10 aromatic, C3/C4 sp3, C1=C2 double bond)
    "C1=Cc2cc3ccccc3nc2CC1": {0: 2, 1: 1, 2: None, 3: 9, 4: None, 5: 8, 6: 7, 7: 6, 8: 5, 9: None, 10: None, 11: None, 12: 4, 13: 3},
    # Phase 643: 1,2-dihydrophenanthridine (14 atoms; N5 aromatic, C1/C2 sp3, C3=C4 double bond)
    "C1=Cc2ncc3ccccc3c2CC1": {0: 3, 1: 4, 2: None, 3: None, 4: 6, 5: None, 6: 7, 7: 8, 8: 9, 9: 10, 10: None, 11: None, 12: 1, 13: 2},
    # Phase 644: 3,4-dihydrophenanthridine (14 atoms; N5 aromatic, C3/C4 sp3, C1=C2 double bond)
    "C1=Cc2c(ncc3ccccc23)CC1": {0: 2, 1: 1, 2: None, 3: None, 4: None, 5: 6, 6: None, 7: 7, 8: 8, 9: 9, 10: 10, 11: None, 12: 4, 13: 3},
    # Phase 633: 4,5,6,7-tetrahydrobenzo-fused 5-membered aromatic heterocycles
    # 4,5,6,7-tetrahydrobenzofuran: O(1)-C(2)-C(3)-C(3a,junc)-C(4)-C(5)-C(6)-C(7)-C(7a,junc)
    "c1cc2c(o1)CCCC2":    {0: 2, 1: 3, 2: None, 3: None, 4: None, 5: 7, 6: 6, 7: 5, 8: 4},
    "c1cc2c(s1)CCCC2":    {0: 2, 1: 3, 2: None, 3: None, 4: None, 5: 7, 6: 6, 7: 5, 8: 4},
    "c1cc2c([nH]1)CCCC2": {0: 2, 1: 3, 2: None, 3: None, 4: 1, 5: 7, 6: 6, 7: 5, 8: 4},
    "c1nc2c([nH]1)CCCC2": {0: 2, 1: 3, 2: None, 3: None, 4: 1, 5: 4, 6: 5, 7: 5, 8: 4},
    "c1n[nH]c2c1CCCC2":   {0: 3, 1: None, 2: 1, 3: None, 4: None, 5: 4, 6: 5, 7: 6, 8: 7},
    # Phase 634: 4,5,6,7-tetrahydrobenzo[d] 5-membered heterocycles
    # benzo[d]oxazole: O(1)-C(2)-N(3,tertiary=None)-C(3a,junc)-C(4)-C(5)-C(6)-C(7)-C(7a,junc)
    "c1nc2c(o1)CCCC2":  {0: 2, 1: None, 2: None, 3: None, 4: None, 5: 7, 6: 6, 7: 5, 8: 4},
    "c1nc2c(s1)CCCC2":  {0: 2, 1: None, 2: None, 3: None, 4: None, 5: 7, 6: 6, 7: 5, 8: 4},
    # benzo[d]isoxazole: O(1)-N(2,None)-C(3)-C(3a,junc)-C(4)-C(5)-C(6)-C(7)-C(7a,junc)
    "c1noc2c1CCCC2":    {0: 3, 1: None, 2: None, 3: None, 4: None, 5: 4, 6: 5, 7: 6, 8: 7},
    "c1nsc2c1CCCC2":    {0: 3, 1: None, 2: None, 3: None, 4: None, 5: 4, 6: 5, 7: 6, 8: 7},
    # 1H-benzo[d][1,2,3]triazole: atom3=C7a(None), atom4=N1(1,H), atom7=C3a(None)
    # C7=atom2(7), C6=atom1(6), C5=atom0(5), C4=atom8(4)
    "C1CCc2[nH]nnc2C1": {0: 5, 1: 6, 2: 7, 3: None, 4: 1, 5: None, 6: None, 7: None, 8: 4},
    # 2H-benzo[d][1,2,3]triazole: C2-symmetric (4≡7 → 4, 5≡6 → 5), N-2(H) at atom5
    "C1CCc2n[nH]nc2C1": {0: 5, 1: 5, 2: 4, 3: None, 4: None, 5: 2, 6: None, 7: None, 8: 4},
    # chromane: O(1)-C(2)-C(3)-C(4)-C(4a,junc)-C(5)-C(6)-C(7)-C(8)-C(8a,junc)
    "c1ccc2c(c1)CCCO2": {0: 6, 1: 7, 2: 8, 3: None, 4: None, 5: 5, 6: 4, 7: 3, 8: 2, 9: 1},
    # isochromane: C(1)-O(2)-C(3)-C(4)-C(4a,junc)-C(5)-C(6)-C(7)-C(8)-C(8a,junc)
    "c1ccc2c(c1)CCOC2": {0: 6, 1: 7, 2: 8, 3: None, 4: None, 5: 5, 6: 4, 7: 3, 8: 2, 9: 1},
    # Phase 407: indolin-2-one: N(1)-C(2,=O)-C(3)-C(3a,junc)-C(4)-C(5)-C(6)-C(7)-C(7a,junc)
    # atom 0=O(exo,None), 1=C2(2), 2=C3(3), 3=C3a(None), 4=C4(4)..7=C7(7), 8=C7a(None), 9=N1(1)
    "O=C1Cc2ccccc2N1":  {0: None, 1: 2, 2: 3, 3: None, 4: 4, 5: 5, 6: 6, 7: 7, 8: None, 9: 1},
    # Phase 407: isoindolin-1-one: C(1,=O)-N(2)-C(3)-C(3a,junc)-C(4)-C(5)-C(6)-C(7)-C(7a,junc)
    # atom 0=O(exo,None), 1=C1(1), 2=C3(3), 3=N2(2), 4=C3a(None), 5=C4(4)..8=C7(7), 9=C7a(None)
    "O=C1NCc2ccccc21":  {0: None, 1: 1, 2: 2, 3: 3, 4: None, 5: 4, 6: 5, 7: 6, 8: 7, 9: None},
    # Phase 409: 3,4-dihydroquinolin-2(1H)-one: N(1)-C(2,=O)-C(3)-C(4)-C(4a,junc)-...-C(8a,junc)
    # atom 0=O(exo,None), 1=C2(2), 2=C3(3), 3=C4(4), 4=C4a(None), 5=C5..8=C8, 9=C8a(None), 10=N1(1)
    "O=C1CCc2ccccc2N1": {0: None, 1: 2, 2: 3, 3: 4, 4: None, 5: 5, 6: 6, 7: 7, 8: 8, 9: None, 10: 1},
    # Phase 409: 3,4-dihydroisoquinolin-1(2H)-one: C(1,=O)-N(2)-C(3)-C(4)-C(4a,junc)-...-C(8a,junc)
    # atom 0=O(exo,None), 1=C1(1), 2=N2(2), 3=C3(3), 4=C4(4), 5=C4a(None), 6=C5..9=C8, 10=C8a(None)
    "O=C1NCCc2ccccc21": {0: None, 1: 1, 2: 2, 3: 3, 4: 4, 5: None, 6: 5, 7: 6, 8: 7, 9: 8, 10: None},
    # Phase 410: 3,4-dihydronaphthalen-1(2H)-one (β-tetralone)
    # atom 0=O(exo,None), 1=C1(1), 2=C2(2), 3=C3(3), 4=C4(4), 5=C4a(None), 6=C5..9=C8, 10=C8a(None)
    "O=C1CCCc2ccccc21": {0: None, 1: 1, 2: 2, 3: 3, 4: 4, 5: None, 6: 5, 7: 6, 8: 7, 9: 8, 10: None},
    # Phase 628: 3,4-dihydronaphthalen-2(1H)-one (beta-tetralone)
    "O=C1CCc2ccccc2C1": {0: None, 1: None, 2: 3, 3: 4, 4: None, 5: 5, 6: 6, 7: 7, 8: 8, 9: None, 10: 1},
    # Phase 410: isochroman-1-one: C(1,=O)-O(2)-C(3)-C(4)-C(4a,junc)-...-C(8a,junc)
    # atom 0=O(exo,None), 1=C1(1), 2=O2(2), 3=C3(3), 4=C4(4), 5=C4a(None), 6=C5..9=C8, 10=C8a(None)
    "O=C1OCCc2ccccc21": {0: None, 1: 1, 2: 2, 3: 3, 4: 4, 5: None, 6: 5, 7: 6, 8: 7, 9: 8, 10: None},
    # Phase 410: chroman-2-one: O(1)-C(2,=O)-C(3)-C(4)-C(4a,junc)-...-C(8a,junc)
    # atom 0=O(exo,None), 1=C2(2), 2=C3(3), 3=C4(4), 4=C4a(None), 5=C5..8=C8, 9=C8a(None), 10=O1(1)
    "O=C1CCc2ccccc2O1": {0: None, 1: 2, 2: 3, 3: 4, 4: None, 5: 5, 6: 6, 7: 7, 8: 8, 9: None, 10: 1},
    # Phase 410: chroman-4-one: O(1)-C(2)-C(3)-C(4,=O)-C(4a,junc)-...-C(8a,junc)
    # atom 0=O(exo,None), 1=C4(4), 2=C3(3), 3=C2(2), 4=O1(1), 5=C8a(None), 6=C8..9=C5, 10=C4a(None)
    "O=C1CCOc2ccccc21": {0: None, 1: 4, 2: 3, 3: 2, 4: 1, 5: None, 6: 8, 7: 7, 8: 6, 9: 5, 10: None},
    # Phase 411: chromone (4H-chromen-4-one): same locant structure as chroman-4-one
    # atom 0=O(exo,None), 1=C4(4), 2=C3(3), 3=C2(2), 4=O1(1), 5=C8a(None), 6=C8..9=C5, 10=C4a(None)
    "O=c1ccoc2ccccc12":  {0: None, 1: 4, 2: 3, 3: 2, 4: 1, 5: None, 6: 8, 7: 7, 8: 6, 9: 5, 10: None},
    # Phase 549: coumarin (2H-chromen-2-one)
    # atom 0=O(exo,None), 1=C2(2), 2=C3(3), 3=C4(4), 4=C4a(None), 5..8=C5..8, 9=C8a(None), 10=O1(None)
    "O=c1ccc2ccccc2o1":  {0: None, 1: 2, 2: 3, 3: 4, 4: None, 5: 5, 6: 6, 7: 7, 8: 8, 9: None, 10: None},
    # Phase 550: isocoumarin (1H-2-benzopyran-1-one)
    # atom 0=O(exo,None), 1=C1(None,carbonyl no-H), 2=O2(None), 3=C3(3), 4=C4(4), 5=C4a(None), 6..9=C5..8, 10=C8a(None)
    "O=c1occc2ccccc12":  {0: None, 1: None, 2: None, 3: 3, 4: 4, 5: None, 6: 5, 7: 6, 8: 7, 9: 8, 10: None},
    # Phase 411: indan-2-one: C1(1)-C2(2,=O)-C3(3)-C3a(None)-C4..C7-C7a(None)
    # atom 0=O(exo,None), 1=C2(2), 2=C1(1), 3=C3a(None), 4=C4..7=C7, 8=C7a(None), 9=C3(3)
    "O=C1Cc2ccccc2C1":   {0: None, 1: 2, 2: 1, 3: None, 4: 4, 5: 5, 6: 6, 7: 7, 8: None, 9: 3},
    # Phase 412: benzazolone C2=O retained names
    # atom 0=O(exo,None), 1=C2(2), 2=N3(3) or heteroatom3, 3=C3a(None), 4..7=C4-C7, 8=C7a(None), 9=N1/S1/O1(1)
    "O=c1[nH]c2ccccc2[nH]1": {0: None, 1: 2, 2: 3, 3: None, 4: 4, 5: 5, 6: 6, 7: 7, 8: None, 9: 1},
    "O=c1[nH]c2ccccc2s1":    {0: None, 1: 2, 2: 3, 3: None, 4: 4, 5: 5, 6: 6, 7: 7, 8: None, 9: 1},
    "O=c1[nH]c2ccccc2o1":    {0: None, 1: 2, 2: 3, 3: None, 4: 4, 5: 5, 6: 6, 7: 7, 8: None, 9: 1},
    # Phase 413: acridin-9(10H)-one: C9(1)-C4b(None)-C1..C4(aromatic)-C4a(None)-N10(10)-C8a(None)-C5..C8-C9a(None)
    # atom 0=O(exo,None), 1=C9(9), 2=C4b(None), 3=C1..6=C4, 7=C4a(None), 8=N10(10), 9=C8a(None), 10=C5..13=C8, 14=C9a(None)
    "O=c1c2ccccc2[nH]c2ccccc12": {0: None, 1: 9, 2: None, 3: 1, 4: 2, 5: 3, 6: 4, 7: None, 8: 10, 9: None, 10: 5, 11: 6, 12: 7, 13: 8, 14: None},
    # Phase 413: quinoxalin-2(1H)-one: N1(10,H)-C2(2,=O)-C3(3)-N4(4)-C4a(None)-C5..C8-C8a(None)
    # atom 0=O(exo,None), 1=C2(2), 2=C3(3), 3=N4(4), 4=C4a(None), 5=C5..8=C8, 9=C8a(None), 10=N1(1,H)
    "O=c1cnc2ccccc2[nH]1":       {0: None, 1: 2, 2: 3, 3: 4, 4: None, 5: 5, 6: 6, 7: 7, 8: 8, 9: None, 10: 1},
    # Phase 414: phthalazin-1(2H)-one: C1(1,=O)-N2(2,H)-N3(3)-C4(4)-C4a(None)-C5..C8-C8a(None)
    "O=c1[nH]ncc2ccccc12": {0: None, 1: 1, 2: 2, 3: 3, 4: 4, 5: None, 6: 5, 7: 6, 8: 7, 9: 8, 10: None},
    # Phase 414: cinnolin-4(1H)-one: C4(1,=O)-C3(3)-N2(2)-N1(1,H)-C8a(None)-C8..C5-C4a(None)
    # atom 0=O(exo,None), 1=C4(4), 2=C3(3), 3=N2(2), 4=N1(1,H), 5=C8a(None), 6=C8..9=C5, 10=C4a(None)
    "O=c1cn[nH]c2ccccc12": {0: None, 1: 4, 2: 3, 3: 2, 4: 1, 5: None, 6: 8, 7: 7, 8: 6, 9: 5, 10: None},
    # Phase 414: quinazolin-4(3H)-one: C4(4,=O)-N3(3,H)-C2(2)-N1(1)-C8a(None)-C8..C5-C4a(None)
    # atom 0=O(exo,None), 1=C4(4), 2=N3(3,H), 3=C2(2), 4=N1(1), 5=C8a(None), 6=C8..9=C5, 10=C4a(None)
    "O=c1[nH]cnc2ccccc12": {0: None, 1: 4, 2: 3, 3: 2, 4: 1, 5: None, 6: 8, 7: 7, 8: 6, 9: 5, 10: None},
    # Phase 415: anthracen-9(10H)-one: 15-atom map; atoms 2,7,9,14 are junctions (None)
    # 0=O(exo), 1=C9, 2=C9a(junc), 3..6=C1..4, 7=C4a(junc), 8=C10, 9=C4b(junc), 10..13=C5..8, 14=C8a(junc)
    "O=C1c2ccccc2Cc2ccccc21": {0: None, 1: 9, 2: None, 3: 1, 4: 2, 5: 3, 6: 4, 7: None, 8: 10, 9: None, 10: 5, 11: 6, 12: 7, 13: 8, 14: None},
    # Phase 551: xanthen-9-one: same 15-atom map; atom8=O (bridge, None), atom1=C9 (carbonyl, None)
    # 0=O(exo), 1=C9(None), 2=C9a(junc), 3..6=C1..4, 7=C4a(junc), 8=O(ring,None), 9=C4b(junc), 10..13=C5..8, 14=C8a(junc)
    "O=c1c2ccccc2oc2ccccc12": {0: None, 1: None, 2: None, 3: 1, 4: 2, 5: 3, 6: 4, 7: None, 8: None, 9: None, 10: 5, 11: 6, 12: 7, 13: 8, 14: None},
    # Phase 551: thioxanthen-9-one: same map, S at atom8
    "O=c1c2ccccc2sc2ccccc12": {0: None, 1: None, 2: None, 3: 1, 4: 2, 5: 3, 6: 4, 7: None, 8: None, 9: None, 10: 5, 11: 6, 12: 7, 13: 8, 14: None},
    # Phase 415: quinolin-2(1H)-one: same map structure as quinoxalin-2(1H)-one
    # 0=O(exo), 1=C2(2), 2=C3(3), 3=C4(4), 4=C4a(junc), 5..8=C5..8, 9=C8a(junc), 10=N1(1,H)
    "O=c1ccc2ccccc2[nH]1":         {0: None, 1: 2, 2: 3, 3: 4, 4: None, 5: 5, 6: 6, 7: 7, 8: 8, 9: None, 10: 1},
    # Phase 415: quinazoline-2,4(1H,3H)-dione: 12-atom map (two exo oxygens)
    # 0=O(exo,C2=O), 1=C2(2), 2=N3(3,H), 3=C4(4), 4=O(exo,C4=O), 5=C4a(junc), 6..9=C5..8, 10=C8a(junc), 11=N1(1,H)
    "O=c1[nH]c(=O)c2ccccc2[nH]1": {0: None, 1: 2, 2: 3, 3: 4, 4: None, 5: None, 6: 5, 7: 6, 8: 7, 9: 8, 10: None, 11: 1},
    # Phase 416: quinolin-4(1H)-one
    # 0=O(exo), 1=C4(4), 2=C3(3), 3=C2(2), 4=N1(1,H), 5=C8a(junc), 6..9=C8..5, 10=C4a(junc)
    "O=c1cc[nH]c2ccccc12": {0: None, 1: 4, 2: 3, 3: 2, 4: 1, 5: None, 6: 8, 7: 7, 8: 6, 9: 5, 10: None},
    # Phase 416: isoquinolin-1(2H)-one
    # 0=O(exo), 1=C1(1), 2=N2(2,H), 3=C3(3), 4=C4(4), 5=C4a(junc), 6..9=C5..8, 10=C8a(junc)
    "O=c1[nH]ccc2ccccc12": {0: None, 1: 1, 2: 2, 3: 3, 4: 4, 5: None, 6: 5, 7: 6, 8: 7, 9: 8, 10: None},
    # Phase 416: isoquinolin-3(2H)-one
    # 0=O(exo), 1=C3(3), 2=C4(4), 3=C4a(junc), 4..7=C5..8, 8=C8a(junc), 9=C1(1), 10=N2(2,H)
    "O=c1cc2ccccc2c[nH]1": {0: None, 1: 3, 2: 4, 3: None, 4: 5, 5: 6, 6: 7, 7: 8, 8: None, 9: 1, 10: 2},
    # Phase 417: 1H-indazol-3(2H)-one (10 atoms, 5-membered ring)
    # 0=O(exo), 1=C3(3), 2=N2(2,H), 3=N1(1,H), 4=C3a(junc), 5..8=C4..7, 9=C7a(junc)
    "O=c1[nH][nH]c2ccccc12": {0: None, 1: 3, 2: 2, 3: 1, 4: None, 5: 4, 6: 5, 7: 6, 8: 7, 9: None},
    # Phase 417: phenanthridin-6(5H)-one (15 atoms, tricyclic)
    # 0=O(exo), 1=C6(6), 2=N5(5,H), 3=C4b(junc), 4..7=C4..1, 8=C4a(junc), 9=C10a(junc), 10..13=C10..7, 14=C6a(junc)
    "O=c1[nH]c2ccccc2c2ccccc12": {0: None, 1: 6, 2: 5, 3: None, 4: 4, 5: 3, 6: 2, 7: 1, 8: None, 9: None, 10: 10, 11: 9, 12: 8, 13: 7, 14: None},
    # Phase 417: phenanthridine (14 atoms, tricyclic)
    # c1ccc2c(c1)cnc1ccccc12: 0=C3, 1=C2, 2=C1, 3=C4a(junc), 4=C4b(junc), 5=C4, 6=C6, 7=N5, 8=C6a(junc), 9..12=C7..10, 13=C10a(junc)
    "c1ccc2c(c1)cnc1ccccc12": {0: 3, 1: 2, 2: 1, 3: None, 4: None, 5: 4, 6: 6, 7: 5, 8: None, 9: 7, 10: 8, 11: 9, 12: 10, 13: None},
    # Phase 564: benzo[f]quinoline (N→5; path 1(a2)-2(a12)-3(a11)-4(a10)-N5(a9)-6(a7)-7(a6)-8(a5)-9(a0)-10(a1))
    "c1ccc2c(c1)ccc1ncccc12": {0: 9, 1: 10, 2: 1, 3: None, 4: None, 5: 8, 6: 7, 7: 6, 8: None, 9: 5, 10: 4, 11: 3, 12: 2, 13: None},
    # Phase 564: benzo[h]quinoline (N→10; path 1(a2)-2(a1)-3(a0)-4(a5)-5(a6)-6(a7)-7(a9)-8(a10)-9(a11)-N10(a12))
    "c1ccc2c(c1)ccc1cccnc12": {0: 3, 1: 2, 2: 1, 3: None, 4: None, 5: 4, 6: 5, 7: 6, 8: None, 9: 7, 10: 8, 11: 9, 12: 10, 13: None},
    # Phase 564: benzo[f]isoquinoline (N→6; path 1(a0)-2(a1)-3(a2)-4(a12)-5(a11)-N6(a10)-7(a9)-8(a7)-9(a6)-10(a5))
    "c1ccc2c(c1)ccc1cnccc12": {0: 1, 1: 2, 2: 3, 3: None, 4: None, 5: 10, 6: 9, 7: 8, 8: None, 9: 7, 10: 6, 11: 5, 12: 4, 13: None},
    # Phase 564: benzo[h]isoquinoline (N→7; path 1(a6)-2(a5)-3(a0)-4(a1)-5(a2)-6(a12)-N7(a11)-8(a10)-9(a9)-10(a7))
    "c1ccc2c(c1)ccc1ccncc12": {0: 3, 1: 4, 2: 5, 3: None, 4: None, 5: 2, 6: 1, 7: 10, 8: None, 9: 9, 10: 8, 11: 7, 12: 6, 13: None},
    # Phase 565: benz[a]acridine (18-atom 4×6 tetracyclic; N→10; path start a9→a10→…→a4(N)→a6→a7)
    "c1ccc2nc3ccc4ccccc4c3cc2c1": {0: 7, 1: 8, 2: 9, 3: None, 4: 10, 5: None, 6: 11, 7: 12, 8: None, 9: 1, 10: 2, 11: 3, 12: 4, 13: None, 14: None, 15: 5, 16: None, 17: 6},
    # Phase 565: benz[c]acridine (18-atom 4×6 tetracyclic; N→7; path start a15→a17→…→a6(N)→a8→…→a13)
    "c1ccc2cc3nc4ccccc4cc3cc2c1": {0: 3, 1: 4, 2: 5, 3: None, 4: 6, 5: None, 6: 7, 7: None, 8: 8, 9: 9, 10: 10, 11: 11, 12: None, 13: 12, 14: None, 15: 1, 16: None, 17: 2},
    # Phase 566: benzo[c]cinnoline (N@6→5, N@7→6; Dir1 from a2; C at {1,2,3,4,7,8,9,10})
    "c1ccc2c(c1)nnc1ccccc12": {0: 3, 1: 2, 2: 1, 3: None, 4: None, 5: 4, 6: 5, 7: 6, 8: None, 9: 7, 10: 8, 11: 9, 12: 10, 13: None},
    # Phase 566: benzo[f]cinnoline (N@9→5, N@10→4; Dir2 from a2; C at {1,2,3,6,7,8,9,10})
    "c1ccc2c(c1)ccc1nnccc12": {0: 9, 1: 10, 2: 1, 3: None, 4: None, 5: 8, 6: 7, 7: 6, 8: None, 9: 5, 10: 4, 11: 3, 12: 2, 13: None},
    # Phase 566: benzo[h]cinnoline (N@11→9, N@12→10; Dir1 from a2; C at {1,2,3,4,5,6,7,8})
    "c1ccc2c(c1)ccc1ccnnc12": {0: 3, 1: 2, 2: 1, 3: None, 4: None, 5: 4, 6: 5, 7: 6, 8: None, 9: 7, 10: 8, 11: 9, 12: 10, 13: None},
    # Phase 566: benzo[f]phthalazine (N@10→4, N@11→3; Dir2 from a2; C at {1,2,5,6,7,8,9,10})
    "c1ccc2c(c1)ccc1cnncc12": {0: 9, 1: 10, 2: 1, 3: None, 4: None, 5: 8, 6: 7, 7: 6, 8: None, 9: 5, 10: 4, 11: 3, 12: 2, 13: None},
    # Phase 566: benzo[f]quinoxaline (N@9→5, N@12→2; Dir2 from a2; C at {1,3,4,6,7,8,9,10})
    "c1ccc2c(c1)ccc1nccnc12": {0: 9, 1: 10, 2: 1, 3: None, 4: None, 5: 8, 6: 7, 7: 6, 8: None, 9: 5, 10: 4, 11: 3, 12: 2, 13: None},
    # Phase 566: benzo[g]cinnoline (linear; N@7→7, N@6→8; acridine pattern; C at {1,2,3,4,5,6,9,10})
    "c1ccc2cc3nnccc3cc2c1": {0: 3, 1: 2, 2: 1, 3: None, 4: 10, 5: None, 6: 8, 7: 7, 8: 6, 9: 5, 10: None, 11: 9, 12: None, 13: 4},
    # Phase 566: benzo[g]phthalazine (linear; N@7→6, N@8→7; acridine pattern; C at {1,2,3,4,5,8,9,10})
    "c1ccc2cc3cnncc3cc2c1": {0: 2, 1: 3, 2: 4, 3: None, 4: 10, 5: None, 6: 5, 7: 6, 8: 7, 9: 8, 10: None, 11: 9, 12: None, 13: 1},
    # Phase 566: benzo[g]quinoxaline (linear; N@6→5, N@9→8; acridine pattern; C at {1,2,3,4,6,7,9,10})
    "c1ccc2cc3nccnc3cc2c1": {0: 2, 1: 3, 2: 4, 3: None, 4: 10, 5: None, 6: 5, 7: 6, 8: 7, 9: 8, 10: None, 11: 9, 12: None, 13: 1},
    # Phase 567: pyrene (D2h; 16 atoms; equivalence classes {1,3,6,8}→pref1, {2,7}→pref2, {4,5,9,10}→pref4)
    # canonical: c1cc2ccc3cccc4ccc(c1)c2c34
    # base_idx atoms 0,7→class{2,7}; 1,6,8,13→class{1,3,6,8}; 3,4,10,11→class{4,5,9,10}; junctions→None
    # All atoms in each class get the minimum (preferred) locant so min(matches) always picks the right one.
    "c1cc2ccc3cccc4ccc(c1)c2c34": {0: 2, 1: 1, 2: None, 3: 4, 4: 4, 5: None, 6: 1, 7: 2, 8: 1, 9: None, 10: 4, 11: 4, 12: None, 13: 1, 14: None, 15: None},
    # Phase 567: triphenylene (D3h; 18 atoms; unique positions {1,5,9}→b5, {2,6,10}→b0)
    # canonical: c1ccc2c(c1)c1ccccc1c1ccccc21
    "c1ccc2c(c1)c1ccccc1c1ccccc21": {0: 2, 1: 3, 2: 4, 3: None, 4: None, 5: 1, 6: None, 7: 12, 8: 11, 9: 10, 10: 9, 11: None, 12: None, 13: 8, 14: 7, 15: 6, 16: 5, 17: None},
    # Phase 567: tetracene (D2h; 18 atoms; unique positions {1,4,7,10}→b17, {2,3,8,9}→b0, {5,6,11,12}→b15)
    # canonical: c1ccc2cc3cc4ccccc4cc3cc2c1
    "c1ccc2cc3cc4ccccc4cc3cc2c1": {0: 2, 1: 3, 2: 4, 3: None, 4: 6, 5: None, 6: 11, 7: None, 8: 7, 9: 8, 10: 9, 11: 10, 12: None, 13: 12, 14: None, 15: 5, 16: None, 17: 1},
    # Phase 579: 1H-naphtho[2,3-d]imidazole (N1H@b6, N3@b8; 7 sub C: 2,4-9)
    # canonical: c1ccc2cc3[nH]cnc3cc2c1; junctions: b3,b5,b9,b11
    "c1ccc2cc3[nH]cnc3cc2c1": {0: 6, 1: 7, 2: 8, 3: None, 4: 9, 5: None, 6: None, 7: 2, 8: None, 9: None, 10: 4, 11: None, 12: 5},
    # Phase 579: 1H-naphtho[2,3-d]pyrazole (N1H@b6, N2@b7; 7 sub C: 3-9)
    # canonical: c1ccc2cc3[nH]ncc3cc2c1
    "c1ccc2cc3[nH]ncc3cc2c1": {0: 6, 1: 7, 2: 8, 3: None, 4: 9, 5: None, 6: None, 7: None, 8: 3, 9: None, 10: 4, 11: None, 12: 5},
    # Phase 579: 1H-naphtho[2,1-d]imidazole (N1H@b10, N3@b11; 7 sub C: 2,4-9)
    # canonical: c1ccc2c(c1)ccc1[nH]cnc12; junctions: b3,b4,b8,b12
    "c1ccc2c(c1)ccc1[nH]cnc12": {0: 7, 1: 8, 2: 9, 3: None, 4: None, 5: 6, 6: 5, 7: 4, 8: None, 9: None, 10: 2, 11: None, 12: None},
    # 3H tautomeric canonical alias (MolFragmentToSmiles gives this form for many substituted cases)
    "c1ccc2c(c1)ccc1nc[nH]c12": {0: 7, 1: 8, 2: 9, 3: None, 4: None, 5: 6, 6: 5, 7: 4, 8: None, 9: None, 10: 2, 11: None, 12: None},
    # Phase 579: 1H-naphtho[2,1-d]pyrazole (N1H@b11, N2@b10; 7 sub C: 3-9)
    # canonical: c1ccc2c(c1)ccc1cn[nH]c12
    "c1ccc2c(c1)ccc1cn[nH]c12": {0: 7, 1: 8, 2: 9, 3: None, 4: None, 5: 6, 6: 5, 7: 4, 8: None, 9: 3, 10: None, 11: None, 12: None},
    # Phase 579: naphtho[1,2-d]oxazole (O@1→b9, N@3→b11; 7 sub C: 2,4-9)
    # canonical: c1ccc2c(c1)ccc1ocnc12; junctions: b3,b4,b8,b12
    "c1ccc2c(c1)ccc1ocnc12": {0: 7, 1: 8, 2: 9, 3: None, 4: None, 5: 6, 6: 5, 7: 4, 8: None, 9: None, 10: 2, 11: None, 12: None},
    # Phase 579: naphtho[1,2-d]thiazole (S@1→b9, N@3→b11; same skeleton)
    # canonical: c1ccc2c(c1)ccc1scnc12
    "c1ccc2c(c1)ccc1scnc12": {0: 7, 1: 8, 2: 9, 3: None, 4: None, 5: 6, 6: 5, 7: 4, 8: None, 9: None, 10: 2, 11: None, 12: None},
    # Phase 581: 1H-naphtho[1,2-d][1,2,3]triazole (N1H@b11, N2@b10, N3@b9; 6 sub C: 4-9)
    # canonical: c1ccc2c(c1)ccc1nn[nH]c12; junctions: b3,b4,b8,b12
    "c1ccc2c(c1)ccc1nn[nH]c12": {0: 7, 1: 8, 2: 9, 3: None, 4: None, 5: 6, 6: 5, 7: 4, 8: None, 9: None, 10: None, 11: None, 12: None},
    # Phase 581: naphtho[2,3-d][1,2,3]oxadiazole (O@b6, N@b7, N@b8; 6 sub C: 4-9)
    "c1ccc2cc3onnc3cc2c1": {0: 6, 1: 7, 2: 8, 3: None, 4: 9, 5: None, 6: None, 7: None, 8: None, 9: None, 10: 4, 11: None, 12: 5},
    # Phase 581: naphtho[1,2-d][1,2,3]oxadiazole (O@b9, N@b10, N@b11; 6 sub C: 4-9)
    "c1ccc2c(c1)ccc1onnc12": {0: 7, 1: 8, 2: 9, 3: None, 4: None, 5: 6, 6: 5, 7: 4, 8: None, 9: None, 10: None, 11: None, 12: None},
    # Phase 581: naphtho[2,1-d][1,2,3]oxadiazole (N@b9, N@b10, O@b11; 6 sub C: 4-9)
    "c1ccc2c(c1)ccc1nnoc12": {0: 7, 1: 8, 2: 9, 3: None, 4: None, 5: 6, 6: 5, 7: 4, 8: None, 9: None, 10: None, 11: None, 12: None},
    # Phase 581: naphtho[2,3-d][1,2,3]thiadiazole (S@b6, N@b7, N@b8; 6 sub C: 4-9)
    "c1ccc2cc3snnc3cc2c1": {0: 6, 1: 7, 2: 8, 3: None, 4: 9, 5: None, 6: None, 7: None, 8: None, 9: None, 10: 4, 11: None, 12: 5},
    # Phase 581: naphtho[1,2-d][1,2,3]thiadiazole (S@b9, N@b10, N@b11; 6 sub C: 4-9)
    "c1ccc2c(c1)ccc1snnc12": {0: 7, 1: 8, 2: 9, 3: None, 4: None, 5: 6, 6: 5, 7: 4, 8: None, 9: None, 10: None, 11: None, 12: None},
    # Phase 581: naphtho[2,1-d][1,2,3]thiadiazole (N@b9, N@b10, S@b11; 6 sub C: 4-9)
    "c1ccc2c(c1)ccc1nnsc12": {0: 7, 1: 8, 2: 9, 3: None, 4: None, 5: 6, 6: 5, 7: 4, 8: None, 9: None, 10: None, 11: None, 12: None},
    # Phase 580: naphtho[2,3-d]isoxazole (O@1→b6, N@2→b7; 7 sub C: 3-9)
    # canonical: c1ccc2cc3oncc3cc2c1; junctions: b3,b5,b9,b11
    "c1ccc2cc3oncc3cc2c1": {0: 6, 1: 7, 2: 8, 3: None, 4: 9, 5: None, 6: None, 7: None, 8: 3, 9: None, 10: 4, 11: None, 12: 5},
    # Phase 580: naphtho[2,3-d]isothiazole (S@1→b6, N@2→b7; same skeleton)
    "c1ccc2cc3sncc3cc2c1": {0: 6, 1: 7, 2: 8, 3: None, 4: 9, 5: None, 6: None, 7: None, 8: 3, 9: None, 10: 4, 11: None, 12: 5},
    # Phase 580: naphtho[2,1-d]isoxazole (C3@b9, N@2→b10, O@1→b11; 7 sub C: 3-9)
    # canonical: c1ccc2c(c1)ccc1cnoc12; junctions: b3,b4,b8,b12
    "c1ccc2c(c1)ccc1cnoc12": {0: 7, 1: 8, 2: 9, 3: None, 4: None, 5: 6, 6: 5, 7: 4, 8: None, 9: 3, 10: None, 11: None, 12: None},
    # Phase 580: naphtho[2,1-d]isothiazole (C3@b9, N@2→b10, S@1→b11; same skeleton)
    "c1ccc2c(c1)ccc1cnsc12": {0: 7, 1: 8, 2: 9, 3: None, 4: None, 5: 6, 6: 5, 7: 4, 8: None, 9: 3, 10: None, 11: None, 12: None},
    # Phase 580: naphtho[1,2-d]isoxazole (O@1→b9, N@2→b10, C3@b11; 7 sub C: 3-9)
    # canonical: c1ccc2c(c1)ccc1oncc12; junctions: b3,b4,b8,b12
    "c1ccc2c(c1)ccc1oncc12": {0: 7, 1: 8, 2: 9, 3: None, 4: None, 5: 6, 6: 5, 7: 4, 8: None, 9: None, 10: None, 11: 3, 12: None},
    # Phase 580: naphtho[1,2-d]isothiazole (S@1→b9, N@2→b10, C3@b11; same skeleton)
    "c1ccc2c(c1)ccc1sncc12": {0: 7, 1: 8, 2: 9, 3: None, 4: None, 5: 6, 6: 5, 7: 4, 8: None, 9: None, 10: None, 11: 3, 12: None},
    # Phase 582: benzo[g]quinoline (N@b6; 9 sub C: 2-10)
    "c1ccc2cc3ncccc3cc2c1": {0: 7, 1: 8, 2: 9, 3: None, 4: 10, 5: None, 6: None, 7: 2, 8: 3, 9: 4, 10: None, 11: 5, 12: None, 13: 6},
    # Phase 583: naphtho[2,3-d]pyrimidine (N@a6/a8; 8 sub C: 2,4-10)
    "c1ccc2cc3ncncc3cc2c1": {0: 7, 1: 8, 2: 9, 3: None, 4: 10, 5: None, 6: None, 7: 2, 8: None, 9: 4, 10: None, 11: 5, 12: None, 13: 6},
    # Phase 583: naphtho[2,1-d]pyrimidine (N@a9/a11; 8 sub C: 1,3,5-10)
    "c1ccc2c(c1)ccc1ncncc12": {0: 8, 1: 9, 2: 10, 3: None, 4: None, 5: 7, 6: 6, 7: 5, 8: None, 9: None, 10: 3, 11: None, 12: 1, 13: None},
    # Phase 583: naphtho[1,2-d]pyrimidine (N@a10/a12; 8 sub C: 2,4-10)
    "c1ccc2c(c1)ccc1cncnc12": {0: 8, 1: 9, 2: 10, 3: None, 4: None, 5: 7, 6: 6, 7: 5, 8: None, 9: 4, 10: None, 11: 2, 12: None, 13: None},
    # Phase 584: 1H-benzo[f]indole (N@a6; 8 sub C: 2-9)
    "c1ccc2cc3[nH]ccc3cc2c1": {0: 6, 1: 7, 2: 8, 3: None, 4: 9, 5: None, 6: None, 7: 2, 8: 3, 9: None, 10: 4, 11: None, 12: 5},
    # Phase 584: 1H-benzo[g]indole (N@a11; 8 sub C: 2-9)
    "c1ccc2c(c1)ccc1cc[nH]c12": {0: 7, 1: 8, 2: 9, 3: None, 4: None, 5: 6, 6: 5, 7: 4, 8: None, 9: 3, 10: 2, 11: None, 12: None},
    # Phase 585: benzo[b][1,7]naphthyridine (N@a4/a7; 8 sub C: 1,3-9)
    "c1ccc2nc3cnccc3cc2c1": {0: 7, 1: 8, 2: 9, 3: None, 4: None, 5: None, 6: 1, 7: None, 8: 3, 9: 4, 10: None, 11: 5, 12: None, 13: 6},
    # Phase 585: benzo[b][1,8]naphthyridine (N@a4/a6; 8 sub C: 2-9)
    "c1ccc2nc3ncccc3cc2c1": {0: 7, 1: 8, 2: 9, 3: None, 4: None, 5: None, 6: None, 7: 2, 8: 3, 9: 4, 10: None, 11: 5, 12: None, 13: 6},
    # Phase 585: benzo[b][1,5]naphthyridine (N@a4/a9; 8 sub C: 2-4,6-10)
    "c1ccc2nc3cccnc3cc2c1": {0: 8, 1: 7, 2: 6, 3: None, 4: None, 5: None, 6: 4, 7: 3, 8: 2, 9: None, 10: None, 11: 10, 12: None, 13: 9},
    # Phase 585: benzo[b][1,6]naphthyridine (N@a4/a8; 8 sub C: 1,3,4,6-10)
    "c1ccc2nc3ccncc3cc2c1": {0: 8, 1: 7, 2: 6, 3: None, 4: None, 5: None, 6: 4, 7: 3, 8: None, 9: 1, 10: None, 11: 10, 12: None, 13: 9},
    # Phase 585: benzo[c][1,6]naphthyridine (N@a7/a11; 8 sub C: 1,3,4,6-10)
    "c1ccc2c(c1)cnc1ccncc12": {0: 8, 1: 9, 2: 10, 3: None, 4: None, 5: 7, 6: 6, 7: None, 8: None, 9: 4, 10: 3, 11: None, 12: 1, 13: None},
    # Phase 586: pyrido[3,4-g]quinoline (N@a2/a8; 8 sub C: 2-6,8-10)
    "c1cnc2cc3ccncc3cc2c1": {0: 3, 1: 2, 2: None, 3: None, 4: 10, 5: None, 6: 9, 7: 8, 8: None, 9: 6, 10: None, 11: 5, 12: None, 13: 4},
    # Phase 586: pyrido[3,2-g]quinoline (N@a2/a6 symm; 5 unique sub C: 2-5,10)
    "c1cnc2cc3ncccc3cc2c1": {0: 3, 1: 2, 2: None, 3: None, 4: 10, 5: None, 6: None, 7: 2, 8: 3, 9: 4, 10: None, 11: 5, 12: None, 13: 4},
    # Phase 586: pyrido[2,3-g]quinoline (N@a2/a9 symm; 4 unique sub C: 2-5)
    "c1cnc2cc3cccnc3cc2c1": {0: 3, 1: 2, 2: None, 3: None, 4: 5, 5: None, 6: 4, 7: 3, 8: 2, 9: None, 10: None, 11: 5, 12: None, 13: 4},
    # Phase 586: pyrido[4,3-g]quinoline (N@a2/a7; 8 sub C: 2-7,9-10)
    "c1cnc2cc3cnccc3cc2c1": {0: 3, 1: 2, 2: None, 3: None, 4: 10, 5: None, 6: 9, 7: None, 8: 7, 9: 6, 10: None, 11: 5, 12: None, 13: 4},
    # Phase 587: pyrido[3,4-g]quinoxaline (N@a2/a8/a9; 7 sub C: 2,3,5,6,8-10)
    "c1cc2cc3nccnc3cc2cn1": {0: 8, 1: 9, 2: None, 3: 10, 4: None, 5: None, 6: 2, 7: 3, 8: None, 9: None, 10: 5, 11: None, 12: 6, 13: None},
    # Phase 587: pyrido[2,3-g]quinoxaline (N@a2/a5/a6; 7 sub C: 2,3,5,7-10)
    "c1cnc2cc3nccnc3cc2c1": {0: 8, 1: 7, 2: None, 3: None, 4: 5, 5: None, 6: None, 7: 3, 8: 2, 9: None, 10: None, 11: 10, 12: None, 13: 9},
    # Phase 587: pyrazino[2,3-b]quinoline (N@a3/a4/a6; 7 sub C: 2,3,6-10)
    "c1ccc2nc3nccnc3cc2c1": {0: 8, 1: 7, 2: 6, 3: None, 4: None, 5: None, 6: None, 7: 3, 8: 2, 9: None, 10: None, 11: 10, 12: None, 13: 9},
    # Phase 588: pyrido[3,2-g]quinazoline (N@a2/a6/a8; 7 sub C: 2,4-8,10)
    "c1cnc2cc3ncncc3cc2c1": {0: 7, 1: 8, 2: None, 3: None, 4: 10, 5: None, 6: None, 7: 2, 8: None, 9: 4, 10: None, 11: 5, 12: None, 13: 6},
    # Phase 588: pyrido[2,3-g]quinazoline (N@a2/a7/a9; 7 sub C: 2,4-5,7-10)
    "c1cnc2cc3cncnc3cc2c1": {0: 8, 1: 7, 2: None, 3: None, 4: 5, 5: None, 6: 4, 7: None, 8: 2, 9: None, 10: None, 11: 10, 12: None, 13: 9},
    # Phase 588: pyrimido[4,5-b]quinoline (N@a3/a6/a8; 7 sub C: 2,4-9)
    "c1ccc2nc3ncncc3cc2c1": {0: 7, 1: 8, 2: 9, 3: None, 4: None, 5: None, 6: None, 7: 2, 8: None, 9: 4, 10: None, 11: 5, 12: None, 13: 6},
    # Phase 588: pyrimido[5,4-b]quinoline (N@a3/a4/a7; 7 sub C: 2,4,6-10)
    "c1ccc2nc3cncnc3cc2c1": {0: 8, 1: 7, 2: 6, 3: None, 4: None, 5: None, 6: 4, 7: None, 8: 2, 9: None, 10: None, 11: 10, 12: None, 13: 9},
    # Phase 589: pyrimido[5,4-g][1,5]naphthyridine (N@a2/a6/a8/a11; 6 sub C: 2,4,6-8,10)
    "c1cnc2cc3ncncc3nc2c1": {0: 7, 1: 8, 2: None, 3: None, 4: 10, 5: None, 6: None, 7: 2, 8: None, 9: 4, 10: None, 11: None, 12: None, 13: 6},
    # Phase 589: pyrimido[5,4-g][1,8]naphthyridine (N@a2/a4/a8/a11; 6 sub C: 2,4-8)
    "c1cnc2nc3ncncc3cc2c1": {0: 7, 1: 8, 2: None, 3: None, 4: None, 5: None, 6: None, 7: 2, 8: None, 9: 4, 10: None, 11: 5, 12: None, 13: 6},
    # Phase 589: pyrido[2,3-g][1,5]naphthyridine (N@a2/a4/a9; 7 sub C: 2-4,7-10)
    "c1cnc2nc3cccnc3cc2c1": {0: 8, 1: 7, 2: None, 3: None, 4: None, 5: None, 6: 4, 7: 3, 8: 2, 9: None, 10: None, 11: 10, 12: None, 13: 9},
    # Phase 590: pyrazino[2,3-g][1,8]naphthyridine (N@a2/a4/a5/a6; 6 sub C: 2,3,7-10)
    "c1cnc2nc3nccnc3cc2c1": {0: 8, 1: 7, 2: None, 3: None, 4: None, 5: None, 6: None, 7: 3, 8: 2, 9: None, 10: None, 11: 10, 12: None, 13: 9},
    # Phase 590: pyrazino[2,3-b][1,5]naphthyridine (N@a2/a4/a5/a9; 6 sub C: 2,3,6-8,10)
    "c1cnc2cc3nccnc3nc2c1": {0: 7, 1: 8, 2: None, 3: None, 4: 10, 5: None, 6: None, 7: 2, 8: 3, 9: None, 10: None, 11: None, 12: None, 13: 6},
    # Phase 590: pyrazino[2,3-b]quinoxaline (N@a3/a4/a5/a11; 3 unique sub C: 2,6,7)
    "c1ccc2nc3nccnc3nc2c1": {0: 7, 1: 7, 2: 6, 3: None, 4: None, 5: None, 6: None, 7: 2, 8: 2, 9: None, 10: None, 11: None, 12: None, 13: 6},
    # Phase 590: pyrazino[2,3-g]quinoxaline (N@a2/a6/a9/a13; 2 unique sub C: 2,5)
    "c1cnc2cc3nccnc3cc2n1": {0: 2, 1: 2, 2: None, 3: None, 4: 5, 5: None, 6: None, 7: 2, 8: 2, 9: None, 10: None, 11: 5, 12: None, 13: None},
    # Phase 591: pyrido[2,3-b]quinoxaline (N@a4/a6/a11; 7 sub C: 2-4,6-9)
    "c1ccc2nc3ncccc3nc2c1": {0: 7, 1: 8, 2: 9, 3: None, 4: None, 5: None, 6: None, 7: 2, 8: 3, 9: 4, 10: None, 11: None, 12: None, 13: 6},
    # Phase 591: pyrido[2,3-b][1,7]naphthyridine (N@a2/a4/a7; 7 sub C: 2-7,9)
    "c1cnc2nc3cnccc3cc2c1": {0: 3, 1: 2, 2: None, 3: None, 4: None, 5: None, 6: 9, 7: None, 8: 7, 9: 6, 10: None, 11: 5, 12: None, 13: 4},
    # Phase 591: pyrido[2,3-b][1,6]naphthyridine (N@a2/a4/a8; 7 sub C: 2-6,8-9)
    "c1cnc2nc3ccncc3cc2c1": {0: 3, 1: 2, 2: None, 3: None, 4: None, 5: None, 6: 9, 7: 8, 8: None, 9: 6, 10: None, 11: 5, 12: None, 13: 4},
    # Phase 591: pyrido[2,3-b][1,8]naphthyridine (N@a2/a4/a6; symm; 4 unique sub C: 2-5)
    "c1cnc2nc3ncccc3cc2c1": {0: 3, 1: 2, 2: None, 3: None, 4: None, 5: None, 6: None, 7: 2, 8: 3, 9: 4, 10: None, 11: 5, 12: None, 13: 4},
    # Phase 591: pyrazino[2,3-h][1,6]naphthyridine (N@a2/a7/a9/a12; 6 sub C: 2,3,6-9)
    "c1cnc2c(c1)cnc1nccnc12": {0: 8, 1: 9, 2: None, 3: None, 4: None, 5: 7, 6: 6, 7: None, 8: None, 9: None, 10: 3, 11: 2, 12: None, 13: None},
    # Phase 592: pyrimido[4,5-b]quinoxaline (N@a4/a6/a8/a11; 6 sub C: 2,4,6-9)
    "c1ccc2nc3ncncc3nc2c1": {0: 7, 1: 8, 2: 9, 3: None, 4: None, 5: None, 6: None, 7: 2, 8: None, 9: 4, 10: None, 11: None, 12: None, 13: 6},
    # Phase 593: pteridine (N@a2/a4/a6/a9; 4 sub C: 2,4,6,7)
    "c1cnc2ncncc2n1": {0: 6, 1: 7, 2: None, 3: None, 4: None, 5: 2, 6: None, 7: 4, 8: None, 9: None},
    # Phase 593: pyrido[2,3-d]pyrimidine (N@a2/a4/a8; 5 sub C: 2,4-7)
    "c1cnc2ncncc2c1": {0: 6, 1: 7, 2: None, 3: None, 4: None, 5: 2, 6: None, 7: 4, 8: None, 9: 5},
    # Phase 593: pyrido[3,4-d]pyrimidine (N@a2/a4/a7; 5 sub C: 2,4-6,8)
    "c1cc2cncnc2cn1": {0: 6, 1: 5, 2: None, 3: 4, 4: None, 5: 2, 6: None, 7: None, 8: 8, 9: None},
    # Phase 593: pyrido[2,3-e]pyrimidine (N@a2/a5/a7; 5 sub C: 2,4,6-8)
    "c1cnc2cncnc2c1": {0: 7, 1: 6, 2: None, 3: None, 4: 4, 5: None, 6: 2, 7: None, 8: None, 9: 8},
    # Phase 593: pyrido[3,4-e]pyrimidine (N@a2/a5/a9; 5 sub C: 2,4,5,7,8)
    "c1cc2ncncc2cn1": {0: 7, 1: 8, 2: None, 3: None, 4: 2, 5: None, 6: 4, 7: None, 8: 5, 9: None},
    # Phase 594: pyrido[2,3-e]pyridazine (N@a2/a6/a7; 5 sub C: 3-4,6-8)
    "c1cnc2ccnnc2c1": {0: 7, 1: 6, 2: None, 3: None, 4: 4, 5: 3, 6: None, 7: None, 8: None, 9: 8},
    # Phase 594: pyrido[2,3-d]pyridazine (N@a2/a5/a6; 5 sub C: 2-5,8)
    "c1cnc2cnncc2c1": {0: 3, 1: 2, 2: None, 3: None, 4: 8, 5: None, 6: None, 7: 5, 8: None, 9: 4},
    # Phase 594: pyrido[3,4-c]pyridazine (N@a5/a6/a9; 5 sub C: 3-6,8)
    "c1cc2ccnnc2cn1": {0: 6, 1: 5, 2: None, 3: 4, 4: 3, 5: None, 6: None, 7: None, 8: 8, 9: None},
    # Phase 594: pyrido[3,4-e]pyridazine (N@a2/a3/a6; 5 sub C: 3-5,7-8)
    "c1cc2nnccc2cn1": {0: 7, 1: 8, 2: None, 3: None, 4: None, 5: 3, 6: 4, 7: None, 8: 5, 9: None},
    # Phase 594: pyrido[2,3-c]pyridazine (N@a2/a3/a8; 5 sub C: 3-7)
    "c1cnc2nnccc2c1": {0: 6, 1: 7, 2: None, 3: None, 4: None, 5: None, 6: 3, 7: 4, 8: None, 9: 5},
    # Phase 594: 1,2,4-benzotriazine (N@a4/a5/a7; 5 sub C: 3,5-8)
    "c1ccc2nncnc2c1": {0: 6, 1: 7, 2: 8, 3: None, 4: None, 5: None, 6: 3, 7: None, 8: None, 9: 5},
    # Phase 594: pyrido[3,4-e]pyrazine (N@a2/a3/a7; 5 sub C: 2,3,5,7,8)
    "c1cc2nccnc2cn1": {0: 7, 1: 8, 2: None, 3: None, 4: 2, 5: 3, 6: None, 7: None, 8: 5, 9: None},
    # Phase 595: indolizine (N@a3/a7; 7 sub C: 1-3,5-8)
    "c1ccn2cccc2c1": {0: 7, 1: 6, 2: 5, 3: None, 4: 3, 5: 2, 6: 1, 7: None, 8: 8},
    # Phase 595: imidazo[1,2-a]pyridine (N@a3/a7/a8; 6 sub C: 2,3,5-8)
    "c1ccn2ccnc2c1": {0: 7, 1: 6, 2: 5, 3: None, 4: 3, 5: 2, 6: None, 7: None, 8: 8},
    # Phase 595: imidazo[1,5-a]pyridine (N@a3/a5/a7; 6 sub C: 1,3,5-8)
    "c1ccn2cncc2c1": {0: 7, 1: 6, 2: 5, 3: None, 4: 3, 5: None, 6: 1, 7: None, 8: 8},
    # Phase 595: pyrazolo[1,5-a]pyridine (N@a3/a4/a7; 5 sub C: 2,3,5-7)
    "c1ccn2nccc2c1": {0: 5, 1: 6, 2: 7, 3: None, 4: None, 5: 2, 6: 3, 7: None, 8: None},
    # Phase 595: [1,2,4]triazolo[4,3-a]pyridine (N@a3/a5/a6/a7; 5 sub C: 3,5-8)
    "c1ccn2cnnc2c1": {0: 7, 1: 6, 2: 5, 3: None, 4: 3, 5: None, 6: None, 7: None, 8: 8},
    # Phase 595: [1,2,4]triazolo[1,5-a]pyridine (N@a3/a4/a5/a7; 5 sub C: 2,5-8)
    "c1ccn2ncnc2c1": {0: 7, 1: 6, 2: 5, 3: None, 4: None, 5: 2, 6: None, 7: None, 8: 8},
    # Phase 595: [1,2,3]triazolo[1,5-a]pyridine (N@a3/a4/a5/a7; 4 sub C: 3,5-7)
    "c1ccn2nncc2c1": {0: 5, 1: 6, 2: 7, 3: None, 4: None, 5: None, 6: 3, 7: None, 8: None},
    # Phase 595: tetrazolo[1,5-a]pyridine (N@a3/a4/a5/a6/a7; 4 sub C: 5-8)
    "c1ccn2nnnc2c1": {0: 7, 1: 6, 2: 5, 3: None, 4: None, 5: None, 6: None, 7: None, 8: 8},
    # Phase 596: imidazo[1,2-b]pyridazine (N@a2/a3/a7/a8; 5 sub C: 2,3,6-8)
    "c1cnn2ccnc2c1": {0: 7, 1: 6, 2: None, 3: None, 4: 3, 5: 2, 6: None, 7: None, 8: 8},
    # Phase 596: imidazo[1,5-b]pyridazine (N@a2/a3/a5/a7; 4 sub C: 2,3,5,7)
    "c1cnn2cncc2c1": {0: 3, 1: 2, 2: None, 3: None, 4: 7, 5: None, 6: 5, 7: None, 8: None},
    # Phase 596: pyrazolo[1,5-b]pyridazine (N@a2/a3/a4/a7; 4 sub C: 2,3,5,6)
    "c1cnn2nccc2c1": {0: 5, 1: 6, 2: None, 3: None, 4: None, 5: 2, 6: 3, 7: None, 8: None},
    # Phase 596: [1,2,4]triazolo[1,5-b]pyridazine (N@a2/a3/a5/a6/a7; 4 sub C: 2,6-8)
    "c1cnn2ncnc2c1": {0: 7, 1: 6, 2: None, 3: None, 4: None, 5: 2, 6: None, 7: None, 8: 8},
    # Phase 596: [1,2,3]triazolo[1,5-b]pyridazine (N@a2/a3/a4/a5/a7; 3 sub C: 3,5,6)
    "c1cnn2nncc2c1": {0: 5, 1: 6, 2: None, 3: None, 4: None, 5: None, 6: 3, 7: None, 8: None},
    # Phase 596: tetrazolo[1,5-b]pyridazine (N@a2/a3/a4/a5/a6/a7; 3 sub C: 6-8)
    "c1cnn2nnnc2c1": {0: 7, 1: 6, 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: 8},
    # Phase 596: imidazo[1,2-a]pyrimidine (N@a2/a3/a8/a9; 5 sub C: 2,3,5-7)
    "c1cnc2nccn2c1": {0: 6, 1: 7, 2: None, 3: None, 4: None, 5: 2, 6: 3, 7: None, 8: 5},
    # Phase 596: imidazo[1,5-a]pyrimidine (N@a2/a5/a8/a9; 4 sub C: 2,3,6,8)
    "c1cnc2cncn2c1": {0: 3, 1: 2, 2: None, 3: None, 4: 8, 5: None, 6: 6, 7: None, 8: None},
    # Phase 597: pyrazolo[1,5-a]pyrimidine (N@a3/a4/a8; 5 sub C: 2,3,5,6,7)
    "c1cnc2ccnn2c1": {0: 6, 1: 5, 2: None, 3: None, 4: 3, 5: 2, 6: None, 7: None, 8: 7},
    # Phase 597: [1,2,4]triazolo[4,3-a]pyrimidine (N@a3/a6/a7/a8; 4 sub C: 3,5,6,7)
    "c1cnc2nncn2c1": {0: 6, 1: 7, 2: None, 3: None, 4: None, 5: None, 6: 3, 7: None, 8: 5},
    # Phase 597: [1,2,4]triazolo[1,5-a]pyrimidine (N@a3/a5/a6/a8; 4 sub C: 2,5,6,7)
    "c1cnc2ncnn2c1": {0: 6, 1: 5, 2: None, 3: None, 4: None, 5: 2, 6: None, 7: None, 8: 7},
    # Phase 597: pyrrolo[1,2-a]pyrazine (N@a4/a8; 6 sub C: 1,3,4,6,7,8)
    "c1cc2cnccn2c1": {0: 7, 1: 8, 2: None, 3: 1, 4: None, 5: 3, 6: 4, 7: None, 8: 6},
    # Phase 597: imidazo[1,2-a]pyrazine (N@a4/a5/a8; 5 sub C: 2,3,5,6,8)
    "c1cn2ccnc2cn1": {0: 6, 1: 5, 2: None, 3: 3, 4: 2, 5: None, 6: None, 7: 8, 8: None},
    # Phase 597: pyrrolo[1,2-b]pyridazine (N@a2/a3/a8; 6 sub C: 2-7)
    "c1cnn2cccc2c1": {0: 3, 1: 2, 2: None, 3: None, 4: 7, 5: 6, 6: 5, 7: None, 8: 4},
    # Phase 597: pyrrolo[1,2-a]pyrimidine (N@a3/a4/a8; 6 sub C: 2-4,6-8)
    "c1cnc2cccn2c1": {0: 3, 1: 2, 2: None, 3: None, 4: 8, 5: 7, 6: 6, 7: None, 8: 4},
    # Phase 598: thieno/furo-pyrazine, [1,2,3]triazolo-pyrimidine, isoxazolo/isothiazolo-pyrimidine/pyrazine
    "c1cnc2sccc2n1": {0: 2, 1: 3, 2: None, 3: None, 4: None, 5: 6, 6: 7, 7: None, 8: None},  # thieno[2,3-e]pyrazine
    "c1cnc2occc2n1": {0: 2, 1: 3, 2: None, 3: None, 4: None, 5: 6, 6: 7, 7: None, 8: None},  # furo[2,3-e]pyrazine
    "c1cnc2cocc2n1": {0: 2, 1: 2, 2: None, 3: None, 4: 5, 5: None, 6: 5, 7: None, 8: None},  # furo[3,4-e]pyrazine
    "c1cnc2occc2c1": {0: 5, 1: 6, 2: None, 3: None, 4: None, 5: 2, 6: 3, 7: None, 8: 4},  # furo[2,3-b]pyridine
    "c1cnc2cscc2n1": {0: 2, 1: 2, 2: None, 3: None, 4: 5, 5: None, 6: 5, 7: None, 8: None},  # thieno[3,4-e]pyrazine
    "c1cnc2cnnn2c1": {0: 6, 1: 5, 2: None, 3: None, 4: 3, 5: None, 6: None, 7: None, 8: 7},  # [1,2,3]triazolo[1,5-a]pyrimidine
    "c1cnc2nocc2n1": {0: 5, 1: 6, 2: None, 3: None, 4: None, 5: None, 6: 3, 7: None, 8: None},  # isoxazolo[3,4-e]pyrazine
    "c1cnc2oncc2n1": {0: 5, 1: 6, 2: None, 3: None, 4: None, 5: None, 6: 3, 7: None, 8: None},  # isoxazolo[4,5-e]pyrazine
    "c1ncc2sncc2n1": {0: 5, 1: None, 2: 7, 3: None, 4: None, 5: None, 6: 3, 7: None, 8: None},  # isothiazolo[4,5-d]pyrimidine
    "c1ncc2nscc2n1": {0: 5, 1: None, 2: 7, 3: None, 4: None, 5: None, 6: 3, 7: None, 8: None},  # isothiazolo[4,3-d]pyrimidine
    "c1ncc2cnsc2n1": {0: 6, 1: None, 2: 4, 3: None, 4: 3, 5: None, 6: None, 7: None, 8: None},  # isothiazolo[5,4-d]pyrimidine
    "c1ncc2csnc2n1": {0: 6, 1: None, 2: 4, 3: None, 4: 3, 5: None, 6: None, 7: None, 8: None},  # isothiazolo[3,4-d]pyrimidine
    "c1cnc2sncc2n1": {0: 5, 1: 6, 2: None, 3: None, 4: None, 5: None, 6: 3, 7: None, 8: None},  # isothiazolo[4,5-e]pyrazine
    "c1cnc2nscc2n1": {0: 5, 1: 6, 2: None, 3: None, 4: None, 5: None, 6: 3, 7: None, 8: None},  # isothiazolo[3,4-e]pyrazine
    # Phase 599: isothiazolo fused with pyridine and pyridazine
    "c1cnc2nscc2c1": {0: 5, 1: 6, 2: None, 3: None, 4: None, 5: None, 6: 3, 7: None, 8: 4},  # isothiazolo[3,4-b]pyridine
    "c1cnc2csnc2c1": {0: 6, 1: 5, 2: None, 3: None, 4: 3, 5: None, 6: None, 7: None, 8: 7},  # isothiazolo[4,3-b]pyridine
    "c1cnc2cnsc2c1": {0: 6, 1: 5, 2: None, 3: None, 4: 3, 5: None, 6: None, 7: None, 8: 7},  # isothiazolo[4,5-b]pyridine
    "c1cnc2sncc2c1": {0: 5, 1: 6, 2: None, 3: None, 4: None, 5: None, 6: 3, 7: None, 8: 4},  # isothiazolo[5,4-b]pyridine
    "c1cc2csnc2cn1": {0: 5, 1: 4, 2: None, 3: 3, 4: None, 5: None, 6: None, 7: 7, 8: None},  # isothiazolo[3,4-c]pyridine
    "c1cc2nscc2cn1": {0: 6, 1: 7, 2: None, 3: None, 4: None, 5: 3, 6: None, 7: 4, 8: None},  # isothiazolo[4,3-c]pyridine
    "c1cc2sncc2cn1": {0: 6, 1: 7, 2: None, 3: None, 4: None, 5: 3, 6: None, 7: 4, 8: None},  # isothiazolo[4,5-c]pyridine
    "c1cc2cnsc2cn1": {0: 5, 1: 4, 2: None, 3: 3, 4: None, 5: None, 6: None, 7: 7, 8: None},  # isothiazolo[5,4-c]pyridine
    "c1cc2csnc2nn1": {0: 5, 1: 4, 2: None, 3: 3, 4: None, 5: None, 6: None, 7: None, 8: None},  # isothiazolo[3,4-c]pyridazine
    "c1cc2nscc2nn1": {0: 6, 1: 7, 2: None, 3: None, 4: None, 5: 3, 6: None, 7: None, 8: None},  # isothiazolo[4,3-c]pyridazine
    "c1cc2sncc2nn1": {0: 6, 1: 7, 2: None, 3: None, 4: None, 5: 3, 6: None, 7: None, 8: None},  # isothiazolo[4,5-c]pyridazine
    "c1cc2cnsc2nn1": {0: 5, 1: 4, 2: None, 3: 3, 4: None, 5: None, 6: None, 7: None, 8: None},  # isothiazolo[5,4-c]pyridazine
    # Phase 600: oxazolo/thiazolo/thieno fused with pyridine/pyrazine, imidazo[1,2-c]pyrimidine
    "c1cnc2cscc2c1": {0: 3, 1: 2, 2: None, 3: None, 4: 7, 5: None, 6: 5, 7: None, 8: 4},  # thieno[3,4-b]pyridine
    "c1cnc2ocnc2c1": {0: 6, 1: 5, 2: None, 3: None, 4: None, 5: 2, 6: None, 7: None, 8: 7},  # oxazolo[5,4-b]pyridine
    "c1cnc2ncoc2c1": {0: 6, 1: 5, 2: None, 3: None, 4: None, 5: 2, 6: None, 7: None, 8: 7},  # oxazolo[4,5-b]pyridine
    "c1cnc2scnc2c1": {0: 6, 1: 5, 2: None, 3: None, 4: None, 5: 2, 6: None, 7: None, 8: 7},  # thiazolo[5,4-b]pyridine
    "c1cnc2ncsc2c1": {0: 6, 1: 5, 2: None, 3: None, 4: None, 5: 2, 6: None, 7: None, 8: 7},  # thiazolo[4,5-b]pyridine
    "c1cnc2ocnc2n1": {0: 5, 1: 6, 2: None, 3: None, 4: None, 5: 2, 6: None, 7: None, 8: None},  # oxazolo[4,5-e]pyrazine
    "c1cnc2scnc2n1": {0: 5, 1: 6, 2: None, 3: None, 4: None, 5: 2, 6: None, 7: None, 8: None},  # thiazolo[4,5-e]pyrazine
    "c1cc2cscc2cn1": {0: 6, 1: 7, 2: None, 3: 1, 4: None, 5: 3, 6: None, 7: 4, 8: None},  # thieno[3,4-c]pyridine
    "c1cc2nccn2cn1": {0: 7, 1: 8, 2: None, 3: None, 4: 2, 5: 3, 6: None, 7: 5, 8: None},  # imidazo[1,2-c]pyrimidine
    # Phase 601: pyrazolo/triazolo/imidazo fused with pyrazine/[1,2,4]triazine (N-bridged 5+6)
    "c1cn2nccc2cn1": {0: 6, 1: 7, 2: None, 3: None, 4: 2, 5: 3, 6: None, 7: 4, 8: None},  # pyrazolo[1,5-a]pyrazine
    "c1cnn2cnnc2n1": {0: 7, 1: 6, 2: None, 3: None, 4: 3, 5: None, 6: None, 7: None, 8: None},  # [1,2,4]triazolo[4,3-b][1,2,4]triazine
    "c1cn2ncnc2cn1": {0: 6, 1: 5, 2: None, 3: None, 4: 2, 5: None, 6: None, 7: 8, 8: None},  # [1,2,4]triazolo[1,5-a]pyrazine
    "c1cnn2cncc2n1": {0: 2, 1: 3, 2: None, 3: None, 4: 6, 5: None, 6: 8, 7: None, 8: None},  # imidazo[1,5-b][1,2,4]triazine
    "c1cnn2ccnc2n1": {0: 3, 1: 2, 2: None, 3: None, 4: 7, 5: 6, 6: None, 7: None, 8: None},  # imidazo[3,2-b][1,2,4]triazine
    "c1cn2cncc2cn1": {0: 6, 1: 5, 2: None, 3: 3, 4: None, 5: 1, 6: None, 7: 8, 8: None},  # imidazo[1,5-a]pyrazine
    "c1cn2cnnc2cn1": {0: 6, 1: 5, 2: None, 3: 3, 4: None, 5: None, 6: None, 7: 8, 8: None},  # [1,2,4]triazolo[4,3-a]pyrazine
    "c1cnn2nccc2n1": {0: 2, 1: 3, 2: None, 3: None, 4: None, 5: 7, 6: 8, 7: None, 8: None},  # pyrazolo[1,5-b][1,2,4]triazine
    # Phase 602: thieno/furo fused with pyridazine
    "c1cc2ccsc2nn1": {0: 3, 1: 4, 2: None, 3: 5, 4: 6, 5: None, 6: None, 7: None, 8: None},  # thieno[2,3-c]pyridazine
    "c1cc2ccoc2nn1": {0: 3, 1: 4, 2: None, 3: 5, 4: 6, 5: None, 6: None, 7: None, 8: None},  # furo[2,3-c]pyridazine
    "c1cc2cscc2nn1": {0: 3, 1: 4, 2: None, 3: 5, 4: None, 5: 7, 6: None, 7: None, 8: None},  # thieno[3,4-c]pyridazine
    "c1cc2cocc2nn1": {0: 3, 1: 4, 2: None, 3: 5, 4: None, 5: 7, 6: None, 7: None, 8: None},  # furo[3,4-c]pyridazine
    "c1cc2cnncc2s1": {0: 2, 1: 3, 2: None, 3: 4, 4: None, 5: None, 6: 7, 7: None, 8: None},  # thieno[3,2-d]pyridazine
    "c1nncc2cscc12": {0: 1, 1: None, 2: None, 3: 1, 4: None, 5: 5, 6: None, 7: 5, 8: None},  # thieno[3,4-d]pyridazine
    "c1cc2cnncc2o1": {0: 2, 1: 3, 2: None, 3: 4, 4: None, 5: None, 6: 7, 7: None, 8: None},  # furo[3,2-d]pyridazine
    "c1nncc2cocc12": {0: 1, 1: None, 2: None, 3: 1, 4: None, 5: 5, 6: None, 7: 5, 8: None},  # furo[3,4-d]pyridazine
    # Phase 603: 5-membered heteroaromatics fused with [1,2,4]triazine
    "c1nnc2ccsc2n1": {0: 3, 1: None, 2: None, 3: None, 4: 7, 5: 6, 6: None, 7: None, 8: None},  # thieno[2,3-e][1,2,4]triazine
    "c1nnc2sccc2n1": {0: 3, 1: None, 2: None, 3: None, 4: None, 5: 6, 6: 5, 7: None, 8: None},  # thieno[3,2-e][1,2,4]triazine
    "c1nnc2cscc2n1": {0: 3, 1: None, 2: None, 3: None, 4: 7, 5: None, 6: 5, 7: None, 8: None},  # thieno[3,4-e][1,2,4]triazine
    "c1nnc2ccoc2n1": {0: 3, 1: None, 2: None, 3: None, 4: 7, 5: 6, 6: None, 7: None, 8: None},  # furo[2,3-e][1,2,4]triazine
    "c1nnc2occc2n1": {0: 3, 1: None, 2: None, 3: None, 4: None, 5: 6, 6: 5, 7: None, 8: None},  # furo[3,2-e][1,2,4]triazine
    "c1nnc2cocc2n1": {0: 3, 1: None, 2: None, 3: None, 4: 7, 5: None, 6: 5, 7: None, 8: None},  # furo[3,4-e][1,2,4]triazine
    "c1nnc2csnc2n1": {0: 3, 1: None, 2: None, 3: None, 4: 7, 5: None, 6: None, 7: None, 8: None},  # isothiazolo[3,4-e][1,2,4]triazine
    "c1nnc2nscc2n1": {0: 5, 1: None, 2: None, 3: None, 4: None, 5: None, 6: 3, 7: None, 8: None},  # isothiazolo[4,3-e][1,2,4]triazine
    "c1nnc2cnsc2n1": {0: 3, 1: None, 2: None, 3: None, 4: 7, 5: None, 6: None, 7: None, 8: None},  # isothiazolo[5,4-e][1,2,4]triazine
    "c1nnc2conc2n1": {0: 3, 1: None, 2: None, 3: None, 4: 7, 5: None, 6: None, 7: None, 8: None},  # isoxazolo[3,4-e][1,2,4]triazine
    "c1nnc2nocc2n1": {0: 5, 1: None, 2: None, 3: None, 4: None, 5: None, 6: 3, 7: None, 8: None},  # isoxazolo[4,3-e][1,2,4]triazine
    "c1nnc2cnoc2n1": {0: 3, 1: None, 2: None, 3: None, 4: 7, 5: None, 6: None, 7: None, 8: None},  # isoxazolo[5,4-e][1,2,4]triazine
    "c1nnc2oncc2n1": {0: 5, 1: None, 2: None, 3: None, 4: None, 5: None, 6: 3, 7: None, 8: None},  # isoxazolo[4,5-e][1,2,4]triazine
    "c1nnc2sncc2n1": {0: 5, 1: None, 2: None, 3: None, 4: None, 5: None, 6: 3, 7: None, 8: None},  # isothiazolo[4,5-e][1,2,4]triazine
    "c1nnc2ocnc2n1": {0: 3, 1: None, 2: None, 3: None, 4: None, 5: 6, 6: None, 7: None, 8: None},  # oxazolo[4,5-e][1,2,4]triazine
    "c1nnc2ncoc2n1": {0: 3, 1: None, 2: None, 3: None, 4: None, 5: 6, 6: None, 7: None, 8: None},  # oxazolo[5,4-e][1,2,4]triazine
    "c1nnc2scnc2n1": {0: 3, 1: None, 2: None, 3: None, 4: None, 5: 6, 6: None, 7: None, 8: None},  # thiazolo[4,5-e][1,2,4]triazine
    "c1nnc2ncsc2n1": {0: 3, 1: None, 2: None, 3: None, 4: None, 5: 6, 6: None, 7: None, 8: None},  # thiazolo[5,4-e][1,2,4]triazine
    # Phase 604: 5-membered heteroaromatics fused with [1,2,3]triazine
    "c1cc2cnnnc2s1": {0: 6, 1: 5, 2: None, 3: 4, 4: None, 5: None, 6: None, 7: None, 8: None},  # thieno[2,3-d][1,2,3]triazine
    "c1cc2nnncc2s1": {0: 6, 1: 7, 2: None, 3: None, 4: None, 5: None, 6: 4, 7: None, 8: None},  # thieno[3,2-d][1,2,3]triazine
    "c1nnnc2cscc12": {0: 4, 1: None, 2: None, 3: None, 4: None, 5: 7, 6: None, 7: 5, 8: None},  # thieno[3,4-d][1,2,3]triazine
    "c1cc2cnnnc2o1": {0: 6, 1: 5, 2: None, 3: 4, 4: None, 5: None, 6: None, 7: None, 8: None},  # furo[2,3-d][1,2,3]triazine
    "c1cc2nnncc2o1": {0: 6, 1: 7, 2: None, 3: None, 4: None, 5: None, 6: 4, 7: None, 8: None},  # furo[3,2-d][1,2,3]triazine
    "c1nnnc2cocc12": {0: 4, 1: None, 2: None, 3: None, 4: None, 5: 7, 6: None, 7: 5, 8: None},  # furo[3,4-d][1,2,3]triazine
    "c1snc2cnnnc12": {0: 7, 1: None, 2: None, 3: None, 4: 4, 5: None, 6: None, 7: None, 8: None},  # isothiazolo[4,3-d][1,2,3]triazine
    "c1nnnc2nscc12": {0: 4, 1: None, 2: None, 3: None, 4: None, 5: None, 6: None, 7: 5, 8: None},  # isothiazolo[3,4-d][1,2,3]triazine
    "c1onc2cnnnc12": {0: 7, 1: None, 2: None, 3: None, 4: 4, 5: None, 6: None, 7: None, 8: None},  # isoxazolo[4,3-d][1,2,3]triazine
    "c1nnnc2nocc12": {0: 4, 1: None, 2: None, 3: None, 4: None, 5: None, 6: None, 7: 5, 8: None},  # isoxazolo[3,4-d][1,2,3]triazine
    "c1nc2cnnnc2s1": {0: 6, 1: None, 2: None, 3: 4, 4: None, 5: None, 6: None, 7: None, 8: None},  # thiazolo[5,4-d][1,2,3]triazine
    "c1nc2nnncc2s1": {0: 6, 1: None, 2: None, 3: None, 4: None, 5: None, 6: 4, 7: None, 8: None},  # thiazolo[4,5-d][1,2,3]triazine
    "c1nc2cnnnc2o1": {0: 6, 1: None, 2: None, 3: 4, 4: None, 5: None, 6: None, 7: None, 8: None},  # oxazolo[5,4-d][1,2,3]triazine
    "c1nc2nnncc2o1": {0: 6, 1: None, 2: None, 3: None, 4: None, 5: None, 6: 4, 7: None, 8: None},  # oxazolo[4,5-d][1,2,3]triazine
    "c1nnnc2sncc12": {0: 4, 1: None, 2: None, 3: None, 4: None, 5: None, 6: None, 7: 5, 8: None},  # isothiazolo[5,4-d][1,2,3]triazine
    "c1noc2cnnnc12": {0: 7, 1: None, 2: None, 3: None, 4: 4, 5: None, 6: None, 7: None, 8: None},  # isoxazolo[4,5-d][1,2,3]triazine
    "c1nsc2cnnnc12": {0: 7, 1: None, 2: None, 3: None, 4: 4, 5: None, 6: None, 7: None, 8: None},  # isothiazolo[4,5-d][1,2,3]triazine
    "c1nnnc2oncc12": {0: 4, 1: None, 2: None, 3: None, 4: None, 5: None, 6: None, 7: 5, 8: None},  # isoxazolo[5,4-d][1,2,3]triazine
    # Phase 605: [1,2,3]oxadiazolo and [1,2,3]thiadiazolo fused bicyclics
    "c1cnc2nnoc2c1": {0: 6, 1: 5, 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: 7},  # [1,2,3]oxadiazolo[4,5-b]pyridine
    "c1cc2onnc2cn1": {0: 6, 1: 7, 2: None, 3: None, 4: None, 5: None, 6: None, 7: 4, 8: None},  # [1,2,3]oxadiazolo[4,5-c]pyridine
    "c1cnc2onnc2n1": {0: 5, 1: 6, 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: None},  # [1,2,3]oxadiazolo[4,5-e]pyrazine
    "c1nncc2onnc12": {0: 4, 1: None, 2: None, 3: 7, 4: None, 5: None, 6: None, 7: None, 8: None},  # [1,2,3]oxadiazolo[5,4-d]pyridazine
    "c1nnc2onnc2n1": {0: 6, 1: None, 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: None},  # [1,2,3]oxadiazolo[4,5-e][1,2,4]triazine
    "c1nnnc2onnc12": {0: 7, 1: None, 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: None},  # [1,2,3]oxadiazolo[5,4-d][1,2,3]triazine
    "c1nnnc2nnoc12": {0: 7, 1: None, 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: None},  # [1,2,3]oxadiazolo[4,5-d][1,2,3]triazine
    "c1cnc2nnsc2c1": {0: 6, 1: 5, 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: 7},  # [1,2,3]thiadiazolo[4,5-b]pyridine
    "c1cc2snnc2cn1": {0: 6, 1: 7, 2: None, 3: None, 4: None, 5: None, 6: None, 7: 4, 8: None},  # [1,2,3]thiadiazolo[4,5-c]pyridine
    "c1cnc2snnc2n1": {0: 5, 1: 6, 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: None},  # [1,2,3]thiadiazolo[4,5-e]pyrazine
    "c1nncc2snnc12": {0: 4, 1: None, 2: None, 3: 7, 4: None, 5: None, 6: None, 7: None, 8: None},  # [1,2,3]thiadiazolo[5,4-d]pyridazine
    "c1nnc2snnc2n1": {0: 6, 1: None, 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: None},  # [1,2,3]thiadiazolo[4,5-e][1,2,4]triazine
    "c1nnnc2snnc12": {0: 7, 1: None, 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: None},  # [1,2,3]thiadiazolo[5,4-d][1,2,3]triazine
    "c1nnnc2nnsc12": {0: 7, 1: None, 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: None},  # [1,2,3]thiadiazolo[4,5-d][1,2,3]triazine
    # Phase 606: [1,2,5]oxadiazolo and [1,2,5]thiadiazolo fused bicyclics
    "c1cnc2nonc2c1": {0: 6, 1: 5, 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: 7},  # [1,2,5]oxadiazolo[3,4-b]pyridine
    "c1cnc2nonc2n1": {0: 5, 1: 5, 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: None},  # [1,2,5]oxadiazolo[3,4-e]pyrazine
    "c1nncc2nonc12": {0: 4, 1: None, 2: None, 3: 4, 4: None, 5: None, 6: None, 7: None, 8: None},  # [1,2,5]oxadiazolo[3,4-d]pyridazine
    "c1nnc2nonc2n1": {0: 6, 1: None, 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: None},  # [1,2,5]oxadiazolo[3,4-e][1,2,4]triazine
    "c1cnc2nsnc2c1": {0: 6, 1: 5, 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: 7},  # [1,2,5]thiadiazolo[3,4-b]pyridine
    "c1cnc2nsnc2n1": {0: 5, 1: 5, 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: None},  # [1,2,5]thiadiazolo[3,4-e]pyrazine
    "c1nncc2nsnc12": {0: 4, 1: None, 2: None, 3: 4, 4: None, 5: None, 6: None, 7: None, 8: None},  # [1,2,5]thiadiazolo[3,4-d]pyridazine
    "c1nnc2nsnc2n1": {0: 6, 1: None, 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: None},  # [1,2,5]thiadiazolo[3,4-e][1,2,4]triazine
    # Phase 607: tetrazolo fusions; isothiazolo/isoxazolo/oxazolo/thiazolo-pyridazine; pyrrolo/triazolo-triazine
    "c1cn2nnnc2cn1": {0: 6, 1: 5, 2: None, 3: None, 4: None, 5: None, 6: None, 7: 8, 8: None},  # tetrazolo[1,5-a]pyrazine
    "c1cnn2nnnc2n1": {0: 7, 1: 6, 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: None},  # tetrazolo[1,5-b][1,2,4]triazine
    "c1nncn2nnnc12": {0: 8, 1: None, 2: None, 3: 5, 4: None, 5: None, 6: None, 7: None, 8: None},  # tetrazolo[1,5-d][1,2,4]triazine
    "c1nncc2sncc12": {0: 4, 1: None, 2: None, 3: 7, 4: None, 5: None, 6: None, 7: 3, 8: None},  # isothiazolo[5,4-d]pyridazine
    "c1nncc2nscc12": {0: 4, 1: None, 2: None, 3: 7, 4: None, 5: None, 6: None, 7: 3, 8: None},  # isothiazolo[3,4-d]pyridazine
    "c1nncc2oncc12": {0: 4, 1: None, 2: None, 3: 7, 4: None, 5: None, 6: None, 7: 3, 8: None},  # isoxazolo[5,4-d]pyridazine
    "c1nncc2nocc12": {0: 4, 1: None, 2: None, 3: 7, 4: None, 5: None, 6: None, 7: 3, 8: None},  # isoxazolo[3,4-d]pyridazine
    "c1nc2cnncc2o1": {0: 2, 1: None, 2: None, 3: 4, 4: None, 5: None, 6: 7, 7: None, 8: None},  # oxazolo[4,5-d]pyridazine
    "c1nc2cnncc2s1": {0: 2, 1: None, 2: None, 3: 4, 4: None, 5: None, 6: 7, 7: None, 8: None},  # thiazolo[4,5-d]pyridazine
    "c1cnn2ncnc2n1": {0: 7, 1: 6, 2: None, 3: None, 4: None, 5: 2, 6: None, 7: None, 8: None},  # [1,2,4]triazolo[1,5-b][1,2,4]triazine
    "c1cc2nccnn2c1": {0: 7, 1: 8, 2: None, 3: None, 4: 2, 5: 3, 6: None, 7: None, 8: 6},  # pyrrolo[1,2-b][1,2,4]triazine
    "c1cnc2nnnn2c1": {0: 6, 1: 5, 2: None, 3: None, 4: None, 5: None, 6: None, 7: None, 8: 7},  # tetrazolo[1,5-a]pyrimidine
    # Phase 608: thieno/furo fused with pyridine (b/c-fusion), benzo[c]thiophene, benzo-selenophene
    "c1cc2sccc2cn1": {0: 6, 1: 7, 2: None, 3: None, 4: 2, 5: 3, 6: None, 7: 4, 8: None},  # thieno[3,2-c]pyridine
    "c1cnc2ccoc2c1": {0: 6, 1: 5, 2: None, 3: None, 4: 3, 5: 2, 6: None, 7: None, 8: 7},  # furo[3,2-b]pyridine
    "c1cnc2cocc2c1": {0: 3, 1: 2, 2: None, 3: None, 4: 7, 5: None, 6: 5, 7: None, 8: 4},  # furo[3,4-b]pyridine
    "c1ccc2cscc2c1": {0: 5, 1: 5, 2: 4, 3: None, 4: 1, 5: None, 6: 1, 7: None, 8: 4},  # benzo[c]thiophene
    "c1ccc2[se]ccc2c1": {0: 5, 1: 6, 2: 7, 3: None, 4: None, 5: 2, 6: 3, 7: None, 8: 4},  # benzo[b]selenophene
    "c1ccc2c[se]cc2c1": {0: 5, 1: 5, 2: 4, 3: None, 4: 1, 5: None, 6: 1, 7: None, 8: 4},  # benzo[c]selenophene
    "c1cc2occc2cn1": {0: 6, 1: 7, 2: None, 3: None, 4: 2, 5: 3, 6: None, 7: 4, 8: None},  # furo[3,2-c]pyridine
    "c1cc2occc2nn1": {0: 3, 1: 4, 2: None, 3: None, 4: 6, 5: 7, 6: None, 7: None, 8: None},  # furo[2,3-e]pyridazine
    "c1cc2ccsc2cn1": {0: 5, 1: 4, 2: None, 3: 3, 4: 2, 5: None, 6: None, 7: 7, 8: None},  # thieno[2,3-c]pyridine
    "c1cc2ccoc2cn1": {0: 5, 1: 4, 2: None, 3: 3, 4: 2, 5: None, 6: None, 7: 7, 8: None},  # furo[2,3-c]pyridine
    "c1cc2cocc2cn1": {0: 6, 1: 7, 2: None, 3: 1, 4: None, 5: 3, 6: None, 7: 4, 8: None},  # furo[3,4-c]pyridine
    "c1cc2sccc2nn1": {0: 3, 1: 4, 2: None, 3: None, 4: 6, 5: 7, 6: None, 7: None, 8: None},  # thieno[3,2-c]pyridazine
    # Phase 609: thiazolo/oxazolo/isoxazolo fused with pyridine and pyridazine (c-fused)
    "c1cc2ncsc2cn1": {0: 6, 1: 7, 2: None, 3: None, 4: 2, 5: None, 6: None, 7: 4, 8: None},  # thiazolo[5,4-c]pyridine
    "c1cc2ncoc2cn1": {0: 6, 1: 7, 2: None, 3: None, 4: 2, 5: None, 6: None, 7: 4, 8: None},  # oxazolo[5,4-c]pyridine
    "c1cc2scnc2cn1": {0: 6, 1: 7, 2: None, 3: None, 4: 2, 5: None, 6: None, 7: 4, 8: None},  # thiazolo[4,5-c]pyridine
    "c1cc2scnc2nn1": {0: 3, 1: 4, 2: None, 3: None, 4: 6, 5: None, 6: None, 7: None, 8: None},  # thiazolo[4,5-c]pyridazine
    "c1cc2ncsc2nn1": {0: 3, 1: 4, 2: None, 3: None, 4: 6, 5: None, 6: None, 7: None, 8: None},  # thiazolo[5,4-c]pyridazine
    "c1cc2ocnc2cn1": {0: 6, 1: 7, 2: None, 3: None, 4: 2, 5: None, 6: None, 7: 4, 8: None},  # oxazolo[4,5-c]pyridine
    "c1cc2cnoc2cn1": {0: 5, 1: 4, 2: None, 3: 3, 4: None, 5: None, 6: None, 7: 7, 8: None},  # isoxazolo[5,4-c]pyridine
    "c1cc2oncc2cn1": {0: 6, 1: 7, 2: None, 3: None, 4: None, 5: 3, 6: None, 7: 4, 8: None},  # isoxazolo[4,5-c]pyridine
    "c1cc2nocc2cn1": {0: 6, 1: 7, 2: None, 3: None, 4: None, 5: 3, 6: None, 7: 4, 8: None},  # isoxazolo[4,3-c]pyridine
    "c1cc2ncoc2nn1": {0: 3, 1: 4, 2: None, 3: None, 4: 6, 5: None, 6: None, 7: None, 8: None},  # oxazolo[5,4-c]pyridazine
    "c1cc2ocnc2nn1": {0: 3, 1: 4, 2: None, 3: None, 4: 6, 5: None, 6: None, 7: None, 8: None},  # oxazolo[4,5-c]pyridazine
    "c1cc2cnoc2nn1": {0: 5, 1: 4, 2: None, 3: 3, 4: None, 5: None, 6: None, 7: None, 8: None},  # isoxazolo[5,4-c]pyridazine
    "c1cc2oncc2nn1": {0: 6, 1: 7, 2: None, 3: None, 4: None, 5: 3, 6: None, 7: None, 8: None},  # isoxazolo[4,5-c]pyridazine
    "c1cc2nocc2nn1": {0: 6, 1: 7, 2: None, 3: None, 4: None, 5: 3, 6: None, 7: None, 8: None},  # isoxazolo[4,3-c]pyridazine
    "c1cc2conc2cn1": {0: 5, 1: 4, 2: None, 3: 3, 4: None, 5: None, 6: None, 7: 7, 8: None},  # isoxazolo[3,4-c]pyridine
    "c1cc2conc2nn1": {0: 5, 1: 4, 2: None, 3: 3, 4: None, 5: None, 6: None, 7: None, 8: None},  # isoxazolo[3,4-c]pyridazine
    # Phase 610: oxazolo/thiazolo fused with pyrimidine (d); isoxazolo fused with pyridine (b)
    "c1ncc2ncoc2n1": {0: 5, 1: None, 2: 7, 3: None, 4: None, 5: 2, 6: None, 7: None, 8: None},  # oxazolo[5,4-d]pyrimidine
    "c1ncc2ocnc2n1": {0: 5, 1: None, 2: 7, 3: None, 4: None, 5: 2, 6: None, 7: None, 8: None},  # oxazolo[4,5-d]pyrimidine
    "c1ncc2ncsc2n1": {0: 5, 1: None, 2: 7, 3: None, 4: None, 5: 2, 6: None, 7: None, 8: None},  # thiazolo[5,4-d]pyrimidine
    "c1ncc2scnc2n1": {0: 5, 1: None, 2: 7, 3: None, 4: None, 5: 2, 6: None, 7: None, 8: None},  # thiazolo[4,5-d]pyrimidine
    "c1cnc2nocc2c1": {0: 5, 1: 6, 2: None, 3: None, 4: None, 5: None, 6: 3, 7: None, 8: 4},  # isoxazolo[3,4-b]pyridine
    "c1cnc2cnoc2c1": {0: 6, 1: 5, 2: None, 3: None, 4: 3, 5: None, 6: None, 7: None, 8: 7},  # isoxazolo[4,5-b]pyridine
    "c1cnc2oncc2c1": {0: 5, 1: 6, 2: None, 3: None, 4: None, 5: None, 6: 3, 7: None, 8: 4},  # isoxazolo[5,4-b]pyridine
    "c1cnc2conc2c1": {0: 6, 1: 5, 2: None, 3: None, 4: 3, 5: None, 6: None, 7: None, 8: 7},  # isoxazolo[4,3-b]pyridine
    # Phase 611: thieno/furo/isoxazolo fused with pyrimidine (d)
    "c1ncc2ccsc2n1": {0: 2, 1: None, 2: 4, 3: None, 4: 5, 5: 6, 6: None, 7: None, 8: None},  # thieno[2,3-d]pyrimidine
    "c1ncc2cscc2n1": {0: 2, 1: None, 2: 4, 3: None, 4: 5, 5: None, 6: 7, 7: None, 8: None},  # thieno[3,4-d]pyrimidine
    "c1ncc2sccc2n1": {0: 2, 1: None, 2: 4, 3: None, 4: None, 5: 6, 6: 7, 7: None, 8: None},  # thieno[3,2-d]pyrimidine
    "c1ncc2ccoc2n1": {0: 2, 1: None, 2: 4, 3: None, 4: 5, 5: 6, 6: None, 7: None, 8: None},  # furo[2,3-d]pyrimidine
    "c1ncc2cocc2n1": {0: 2, 1: None, 2: 4, 3: None, 4: 5, 5: None, 6: 7, 7: None, 8: None},  # furo[3,4-d]pyrimidine
    "c1ncc2occc2n1": {0: 2, 1: None, 2: 4, 3: None, 4: None, 5: 6, 6: 7, 7: None, 8: None},  # furo[3,2-d]pyrimidine
    "c1ncc2nocc2n1": {0: 5, 1: None, 2: 7, 3: None, 4: None, 5: None, 6: 3, 7: None, 8: None},  # isoxazolo[4,3-d]pyrimidine
    "c1ncc2oncc2n1": {0: 5, 1: None, 2: 7, 3: None, 4: None, 5: None, 6: 3, 7: None, 8: None},  # isoxazolo[4,5-d]pyrimidine
    "c1ncc2conc2n1": {0: 6, 1: None, 2: 4, 3: None, 4: 3, 5: None, 6: None, 7: None, 8: None},  # isoxazolo[3,4-d]pyrimidine
    "c1ncc2cnoc2n1": {0: 6, 1: None, 2: 4, 3: None, 4: 3, 5: None, 6: None, 7: None, 8: None},  # isoxazolo[5,4-d]pyrimidine
    # Phase 578: naphtho[2,3-d]oxazole (O@1→b6, N@3→b8; 7 sub C: 2,4-9)
    # canonical: c1ccc2cc3ocnc3cc2c1; junctions: b3,b5,b9,b11
    "c1ccc2cc3ocnc3cc2c1": {0: 6, 1: 7, 2: 8, 3: None, 4: 9, 5: None, 6: None, 7: 2, 8: None, 9: None, 10: 4, 11: None, 12: 5},
    # Phase 578: naphtho[2,3-d]thiazole (S@1→b6, N@3→b8; same skeleton)
    # canonical: c1ccc2cc3scnc3cc2c1
    "c1ccc2cc3scnc3cc2c1": {0: 6, 1: 7, 2: 8, 3: None, 4: 9, 5: None, 6: None, 7: 2, 8: None, 9: None, 10: 4, 11: None, 12: 5},
    # Phase 578: naphtho[2,1-d]oxazole (N@3→b9, O@1→b11; 7 sub C: 2,4-9)
    # canonical: c1ccc2c(c1)ccc1ncoc12; junctions: b3,b4,b8,b12
    "c1ccc2c(c1)ccc1ncoc12": {0: 7, 1: 8, 2: 9, 3: None, 4: None, 5: 6, 6: 5, 7: 4, 8: None, 9: None, 10: 2, 11: None, 12: None},
    # Phase 578: naphtho[2,1-d]thiazole (N@3→b9, S@1→b11; same skeleton)
    # canonical: c1ccc2c(c1)ccc1ncsc12
    "c1ccc2c(c1)ccc1ncsc12": {0: 7, 1: 8, 2: 9, 3: None, 4: None, 5: 6, 6: 5, 7: 4, 8: None, 9: None, 10: 2, 11: None, 12: None},
    # Phase 577: naphtho[1,2-b]furan (O@1→b11; 8 substitutable C at locants 2-9)
    # canonical: c1ccc2c(c1)ccc1ccoc12; junctions: b3,b4,b8,b12; O: b11
    "c1ccc2c(c1)ccc1ccoc12": {0: 7, 1: 8, 2: 9, 3: None, 4: None, 5: 6, 6: 5, 7: 4, 8: None, 9: 3, 10: 2, 11: None, 12: None},
    # Phase 577: naphtho[1,2-b]thiophene (S@1→b11; same skeleton as furan)
    # canonical: c1ccc2c(c1)ccc1ccsc12
    "c1ccc2c(c1)ccc1ccsc12": {0: 7, 1: 8, 2: 9, 3: None, 4: None, 5: 6, 6: 5, 7: 4, 8: None, 9: 3, 10: 2, 11: None, 12: None},
    # Phase 577: naphtho[2,3-b]furan (O@1→b6; 8 substitutable C at locants 2-9)
    # canonical: c1ccc2cc3occc3cc2c1; junctions: b3,b5,b9,b11; O: b6
    "c1ccc2cc3occc3cc2c1": {0: 6, 1: 7, 2: 8, 3: None, 4: 9, 5: None, 6: None, 7: 2, 8: 3, 9: None, 10: 4, 11: None, 12: 5},
    # Phase 577: naphtho[2,3-b]thiophene (S@1→b6; same skeleton as furan)
    # canonical: c1ccc2cc3sccc3cc2c1
    "c1ccc2cc3sccc3cc2c1": {0: 6, 1: 7, 2: 8, 3: None, 4: 9, 5: None, 6: None, 7: 2, 8: 3, 9: None, 10: 4, 11: None, 12: 5},
    # Phase 576: naphtho[2,1-b]furan (O@3→b9; 8 substitutable C at locants 1,2,4-9)
    # canonical: c1ccc2c(c1)ccc1occc12; junctions: b3,b4,b8,b12; O: b9
    "c1ccc2c(c1)ccc1occc12": {0: 7, 1: 8, 2: 9, 3: None, 4: None, 5: 6, 6: 5, 7: 4, 8: None, 9: None, 10: 2, 11: 1, 12: None},
    # Phase 576: naphtho[2,1-b]thiophene (S@3→b9; same skeleton as furan)
    # canonical: c1ccc2c(c1)ccc1sccc12
    "c1ccc2c(c1)ccc1sccc12": {0: 7, 1: 8, 2: 9, 3: None, 4: None, 5: 6, 6: 5, 7: 4, 8: None, 9: None, 10: 2, 11: 1, 12: None},
    # Phase 575: 9H-pyrido[2,3-b]indole (α-carboline; N@1→b8, N-H@9→b6; C positions 2-8)
    # canonical: c1ccc2c(c1)[nH]c1ncccc12
    "c1ccc2c(c1)[nH]c1ncccc12": {0: 7, 1: 6, 2: 5, 3: None, 4: None, 5: 8, 6: None, 7: None, 8: None, 9: 2, 10: 3, 11: 4, 12: None},
    # Phase 575: 9H-pyrido[3,4-b]indole (β-carboline; N@2→b9, N-H@9→b6; C positions 1,3-8)
    # canonical: c1ccc2c(c1)[nH]c1cnccc12
    "c1ccc2c(c1)[nH]c1cnccc12": {0: 7, 1: 6, 2: 5, 3: None, 4: None, 5: 8, 6: None, 7: None, 8: 1, 9: None, 10: 3, 11: 4, 12: None},
    # Phase 574: 9H-fluoren-9-one (C2v; 14 atoms; 4 unique pairs: {1,8}→b3, {2,7}→b4, {3,6}→b10, {4,5}→b9)
    # canonical: O=C1c2ccccc2-c2ccccc21
    # junctions: 2(→ring1-ring2), 7(→ring1-5ring), 8(→ring2-5ring), 13(→ring2-C9); O at 0, C9=O at 1
    "O=C1c2ccccc2-c2ccccc21": {0: None, 1: None, 2: None, 3: 1, 4: 2, 5: 6, 6: 5, 7: None, 8: None, 9: 4, 10: 3, 11: 7, 12: 8, 13: None},
    # Phase 574: 1H-phenalene (no symmetry; 13 atoms; sp3 CH2 at b12=pos1, sp2 C=C at b0=pos2,b1=pos3; aromatic 4-9)
    # canonical: C1=Cc2cccc3cccc(c23)C1
    # junctions: 2,6,10,11; sp3: 12 (pos1); sp2: 0,1 (pos2,3); aromatic CH: 3,4,5,7,8,9
    "C1=Cc2cccc3cccc(c23)C1": {0: 2, 1: 3, 2: None, 3: 4, 4: 5, 5: 6, 6: None, 7: 7, 8: 8, 9: 9, 10: None, 11: None, 12: 1},
    # Phase 573: coronene (D6h; 24 atoms; all 12 CH positions equivalent → locant 1)
    # canonical: c1cc2ccc3ccc4ccc5ccc6ccc1c1c2c3c4c5c61
    # CH: [0,1,3,4,6,7,9,10,12,13,15,16]; junctions: [2,5,8,11,14,17,18-23]
    "c1cc2ccc3ccc4ccc5ccc6ccc1c1c2c3c4c5c61": {0: 1, 1: 1, 2: None, 3: 1, 4: 1, 5: None, 6: 1, 7: 1, 8: None, 9: 1, 10: 1, 11: None, 12: 1, 13: 1, 14: None, 15: 1, 16: 1, 17: None, 18: None, 19: None, 20: None, 21: None, 22: None, 23: None},
    # Phase 572: benzo[g]isoquinoline (N@2→b7; 14 atoms; 9 unique C positions 1,3-10)
    # canonical: c1ccc2cc3cnccc3cc2c1
    # junctions: 3,5,10,12
    "c1ccc2cc3cnccc3cc2c1": {0: 7, 1: 8, 2: 9, 3: None, 4: 10, 5: None, 6: 1, 7: None, 8: 3, 9: 4, 10: None, 11: 5, 12: None, 13: 6},
    # Phase 571: acenaphthylene (C2v; 12 atoms; 4 unique pairs: {1,2}→b0, {3,8}→b9, {4,7}→b4, {5,6}→b5)
    # canonical: C1=Cc2cccc3cccc1c23
    # C2 axis through midpt(0-1) and junc 11; pairs: 0↔1, 2↔10, 9↔3, 4↔8, 5↔7; junc 6,11 on axis
    "C1=Cc2cccc3cccc1c23": {0: 1, 1: 2, 2: None, 3: 8, 4: 4, 5: 5, 6: None, 7: 6, 8: 7, 9: 3, 10: None, 11: None},
    # Phase 571: acenaphthene (C2v; 12 atoms; 4 unique pairs: {1,2}→b10, {3,8}→b1, {4,7}→b0, {5,6}→b9)
    # canonical: c1cc2c3c(cccc3c1)CC2
    # C2 axis through midpt(CH2 bond 10-11) and junc 3,8; pairs: 10↔11, 2↔4, 1↔5, 0↔6, 9↔7
    "c1cc2c3c(cccc3c1)CC2": {0: 4, 1: 3, 2: None, 3: None, 4: None, 5: 8, 6: 7, 7: 6, 8: None, 9: 5, 10: 1, 11: 2},
    # Phase 570: perylene (D2h; 20 atoms; 3 unique classes: A={1,6,7,12}→b17, B={2,5,8,11}→b0, C={3,4,9,10}→b1)
    # canonical: c1cc2cccc3c4cccc5cccc(c(c1)c23)c54
    # rings: end1=[0,1,2,16,17,18], side2=[2,3,4,5,6,18], central=[6,7,15,16,18,19], side4=[7,8,9,10,11,19], end5=[11,12,13,14,15,19]
    "c1cc2cccc3c4cccc5cccc(c(c1)c23)c54": {0: 2, 1: 3, 2: None, 3: 4, 4: 5, 5: 6, 6: None, 7: None, 8: 7, 9: 8, 10: 9, 11: None, 12: 10, 13: 11, 14: 12, 15: None, 16: None, 17: 1, 18: None, 19: None},
    # Phase 569: benzo[a]pyrene (no symmetry; 20 atoms; 12 CH at locants 1-12)
    # canonical: c1ccc2c(c1)cc1ccc3cccc4ccc2c1c34
    # junctions: 3,4,7,10,14,17,18,19
    "c1ccc2c(c1)cc1ccc3cccc4ccc2c1c34": {0: 8, 1: 9, 2: 10, 3: None, 4: None, 5: 7, 6: 6, 7: None, 8: 5, 9: 4, 10: None, 11: 3, 12: 2, 13: 1, 14: None, 15: 12, 16: 11, 17: None, 18: None, 19: None},
    # Phase 568: fluoranthene (C2v; 16 atoms; unique positions {1,6}→b7, {2,5}→b12, {3,4}→b11, {7,10}→b5, {8,9}→b0)
    # canonical: c1ccc2c(c1)-c1cccc3cccc-2c13
    # ring1=[0-5], ring2=[6,7,8,9,10,15], ring3=[10,11,12,13,14,15], 5ring=[3,4,6,14,15]
    # C2 axis through midpt(3-4) and midpt(bond 6-14 across 5-ring center 15): 0↔1, 2↔5, 3↔4, 7↔13, 8↔12, 9↔11
    "c1ccc2c(c1)-c1cccc3cccc-2c13": {0: 8, 1: 9, 2: 10, 3: None, 4: None, 5: 7, 6: None, 7: 1, 8: 5, 9: 4, 10: None, 11: 3, 12: 2, 13: 6, 14: None, 15: None},
    # Phase 418: 1H-indole-2,3-dione (isatin) — 11 atoms, two exo O
    # 0=O(C2=O), 1=C2(2), 2=N1(1,H), 3=C7a(junc), 4..7=C4..7, 8=C3a(junc), 9=C3(3), 10=O(C3=O)
    "O=C1Nc2ccccc2C1=O": {0: None, 1: 2, 2: 1, 3: None, 4: 4, 5: 5, 6: 6, 7: 7, 8: None, 9: 3, 10: None},
    # Phase 418: benzofuran-2(3H)-one — 10 atoms
    # 0=O(exo,C2=O), 1=C2(2), 2=C3(3,CH2), 3=O1(1), 4=C7a(junc), 5..8=C7..4, 9=C3a(junc)
    "O=C1COc2ccccc21": {0: None, 1: 2, 2: 3, 3: 1, 4: None, 5: 7, 6: 6, 7: 5, 8: 4, 9: None},
    # Phase 418: benzo[b]thiophen-2(3H)-one — same map but S at position 1
    "O=C1CSc2ccccc21": {0: None, 1: 2, 2: 3, 3: 1, 4: None, 5: 7, 6: 6, 7: 5, 8: 4, 9: None},
    # Phase 419: naphthalen-1(2H)-one — 11 atoms, two 6-membered rings
    # 0=O(exo), 1=C1(1,keto), 2=C2(2,CH2), 3=C3(3), 4=C4(4), 5=C4a(junc), 6..9=C5..8, 10=C8a(junc)
    "O=C1CC=Cc2ccccc21": {0: None, 1: 1, 2: 2, 3: 3, 4: 4, 5: None, 6: 5, 7: 6, 8: 7, 9: 8, 10: None},
    # Phase 419: naphthalen-2(1H)-one — 11 atoms
    # 0=O(exo), 1=C2(2,keto), 2=C3(3), 3=C4(4), 4=C4a(junc), 5..8=C5..8, 9=C8a(junc), 10=C1(1,CH2)
    "O=C1C=Cc2ccccc2C1": {0: None, 1: 2, 2: 3, 3: 4, 4: None, 5: 5, 6: 6, 7: 7, 8: 8, 9: None, 10: 1},
    # Phase 419: 1,3-benzodioxol-2-one — 10 atoms (5+6 ring, aromatic)
    # 0=O(exo), 1=C2(2), 2=O1(1), 3=C7a(junc), 4..7=C7..4, 8=C3a(junc), 9=O3(3)
    "O=c1oc2ccccc2o1": {0: None, 1: 2, 2: 1, 3: None, 4: 7, 5: 6, 6: 5, 7: 4, 8: None, 9: 3},
    # Phase 420: 2H-pyran-2-one (alpha-pyrone) — 7 atoms (monocyclic 6-membered)
    # 0=O(exo), 1=C2(2,keto), 2=C3(3), 3=C4(4), 4=C5(5), 5=C6(6), 6=O1(1,ring)
    "O=c1cccco1": {0: None, 1: 2, 2: 3, 3: 4, 4: 5, 5: 6, 6: 1},
    # Phase 420: 4H-pyran-4-one (gamma-pyrone) — 7 atoms (monocyclic 6-membered)
    # 0=O(exo), 1=C4(4,keto), 2=C3(3), 3=C2(2), 4=O1(1,ring), 5=C6(6), 6=C5(5)
    "O=c1ccocc1": {0: None, 1: 4, 2: 3, 3: 2, 4: 1, 5: 6, 6: 5},
}


def _try_fused_hetero_retained(graph: "MoleculeGraph") -> str | None:
    """
    縮合ヘテロ芳香族保留名テーブルに一致する場合、名前を返す。
    置換体にも対応: 環コア SMILES でマッチし、置換基を命名する。
    """
    if graph.rdkit_mol is None:
        return None
    from rdkit.Chem import MolToSmiles, MolFromSmiles, MolFragmentToSmiles
    canon = MolToSmiles(graph.rdkit_mol)
    # 完全一致チェック (置換なし)
    name = _FUSED_HETERO_RETAINED.get(canon)
    if name is not None:
        return name

    # 置換体チェック: 環原子のみの core SMILES を取得して照合
    all_ring_atoms = sorted({a for ring in graph.ring_atom_sets for a in ring})
    if not all_ring_atoms:
        return None
    core_smi_raw = MolFragmentToSmiles(graph.rdkit_mol, all_ring_atoms, canonical=True)
    # MolFragmentToSmiles は必ずしも MolToSmiles と同じ正準形にならないので正規化する
    _tmp_mol = MolFromSmiles(core_smi_raw)
    if _tmp_mol is not None:
        core_smi = MolToSmiles(_tmp_mol)
    else:
        # Phase 158: N-置換縮合環 (1-methylindole 等): 各 'n' を順に '[nH]' に置換して再試行
        import re as _re
        _bare_n_positions = [m.start() for m in _re.finditer(r"(?<!\[)n(?!\])", core_smi_raw)]
        _found_core = None
        for _pos in _bare_n_positions:
            _alt_raw = (core_smi_raw[:_pos] + "[nH]"
                        + core_smi_raw[_pos + 1:])
            _tmp_alt = MolFromSmiles(_alt_raw)
            if _tmp_alt is not None:
                _alt_canon = MolToSmiles(_tmp_alt)
                if _alt_canon in _FUSED_HETERO_RETAINED:
                    _found_core = _alt_canon
                    break
        core_smi = _found_core if _found_core is not None else core_smi_raw

    # Phase 407: also try extended core (ring atoms + exo C=O oxygens on ring carbons)
    # Needed for retained names whose parent includes an exo C=O (e.g., indolin-2-one, isoindolin-1-one)
    _ring_set_407 = set(all_ring_atoms)
    _exo_co_o_407: set[int] = set()
    for _ra_407 in all_ring_atoms:
        _a_407 = graph.rdkit_mol.GetAtomWithIdx(_ra_407)
        if _a_407.GetSymbol() == "C":
            for _nb_407 in _a_407.GetNeighbors():
                _ni_407 = _nb_407.GetIdx()
                if _ni_407 not in _ring_set_407:
                    _bond_407 = graph.rdkit_mol.GetBondBetweenAtoms(_ra_407, _ni_407)
                    if _nb_407.GetSymbol() == "O" and _bond_407.GetBondTypeAsDouble() == 2.0:
                        _exo_co_o_407.add(_ni_407)
    if _exo_co_o_407:
        _ext_atoms_407 = sorted(all_ring_atoms) + sorted(_exo_co_o_407)
        _ext_raw_407 = MolFragmentToSmiles(graph.rdkit_mol, _ext_atoms_407, canonical=True)
        _ext_mol_407 = MolFromSmiles(_ext_raw_407)
        if _ext_mol_407 is not None:
            _ext_smi_407 = MolToSmiles(_ext_mol_407)
            if _ext_smi_407 in _FUSED_HETERO_RETAINED and _ext_smi_407 in _FUSED_LOCANT_MAP:
                core_smi = _ext_smi_407

    base_name = _FUSED_HETERO_RETAINED.get(core_smi)
    locant_map_def = _FUSED_LOCANT_MAP.get(core_smi)

    # Phase 548: 全環原子でマッチしない場合、縮合サブクラスターを試す
    # 例: flavone (2-phenylchromone) — chromone + phenyl substituent
    if base_name is None or locant_map_def is None:
        _ring_list = list(graph.ring_atom_sets)
        _n = len(_ring_list)
        _rsets = [set(r) for r in _ring_list]
        # ring adjacency (fused = share ≥2 atoms)
        _radj: dict[int, list[int]] = {i: [] for i in range(_n)}
        for _i in range(_n):
            for _j in range(_i + 1, _n):
                if len(_rsets[_i] & _rsets[_j]) >= 2:
                    _radj[_i].append(_j)
                    _radj[_j].append(_i)
        # find connected components (fused clusters)
        _vis: set[int] = set()
        _clusters: list[frozenset[int]] = []
        for _s in range(_n):
            if _s in _vis:
                continue
            _cl: list[int] = []
            _stk = [_s]
            while _stk:
                _nd = _stk.pop()
                if _nd in _vis:
                    continue
                _vis.add(_nd)
                _cl.append(_nd)
                _stk.extend(_radj[_nd])
            _clusters.append(frozenset(_cl))
        _all_cl = frozenset(range(_n))
        for _cluster in sorted(_clusters, key=lambda c: -sum(len(_ring_list[i]) for i in c)):
            if _cluster == _all_cl:
                continue  # already tried the full ring set
            _cat: set[int] = set()
            for _ri in _cluster:
                _cat |= _rsets[_ri]
            _exo_o: set[int] = set()
            for _ra in _cat:
                _am = graph.rdkit_mol.GetAtomWithIdx(_ra)
                if _am.GetSymbol() == "C":
                    for _nb in _am.GetNeighbors():
                        _ni = _nb.GetIdx()
                        if _ni not in _cat:
                            _bd = graph.rdkit_mol.GetBondBetweenAtoms(_ra, _ni)
                            if _nb.GetSymbol() == "O" and _bd.GetBondTypeAsDouble() == 2.0:
                                _exo_o.add(_ni)
            _ext_a = sorted(_cat) + sorted(_exo_o)
            _ext_r = MolFragmentToSmiles(graph.rdkit_mol, _ext_a, canonical=True)
            _ext_m = MolFromSmiles(_ext_r)
            if _ext_m is None:
                continue
            _ext_s = MolToSmiles(_ext_m)
            _bn = _FUSED_HETERO_RETAINED.get(_ext_s)
            _lm = _FUSED_LOCANT_MAP.get(_ext_s)
            if _bn is None or _lm is None:
                continue
            # found — use this sub-cluster as base
            core_smi = _ext_s
            base_name = _bn
            locant_map_def = _lm
            all_ring_atoms = sorted(_cat)
            break

    if base_name is None or locant_map_def is None:
        return None

    # substructure match: base_mol_atom_idx → graph_rdkit_atom_idx
    base_mol = MolFromSmiles(core_smi)
    if base_mol is None:
        return None
    # Use SMARTS query when core contains [nH] AND the target mol has a ring NH,
    # to enforce N-H position and prevent tautomeric-flip mismatch.
    # Skip SMARTS when the NH is absent (N-substituted compounds like 1-methylindole).
    if "[nH]" in core_smi and any(
        a.GetSymbol() == "N" and a.GetTotalNumHs() > 0 and a.GetIsAromatic()
        for a in graph.rdkit_mol.GetAtoms()
    ):
        from rdkit.Chem import MolFromSmarts as _MolFromSmarts
        _smarts_q = _MolFromSmarts(core_smi)
        all_matches = graph.rdkit_mol.GetSubstructMatches(_smarts_q, uniquify=False) if _smarts_q is not None else graph.rdkit_mol.GetSubstructMatches(base_mol, uniquify=False)
    else:
        all_matches = graph.rdkit_mol.GetSubstructMatches(base_mol, uniquify=False)
    if not all_matches:
        return None

    # For symmetric compounds pick the match giving the minimum locant set (IUPAC rule)
    _LOCANT_PENALTY = 10000  # penalizes matches where a substituent lands on a None-locant atom

    def _sub_locants(m: tuple) -> list:
        m_set = set(m)
        locs = []
        for base_idx, rdkit_idx in enumerate(m):
            loc = locant_map_def.get(base_idx)
            for nb in graph.rdkit_mol.GetAtomWithIdx(rdkit_idx).GetNeighbors():
                if nb.GetIdx() not in m_set and nb.GetAtomicNum() != 1:
                    locs.append(loc if loc is not None else _LOCANT_PENALTY)
        return sorted(locs)

    match = min(all_matches, key=_sub_locants)
    match_set = set(match)  # atoms belonging to the base structure (not substituents)
    # match[base_idx] = rdkit_atom_idx in graph.rdkit_mol

    # rdkit atom idx → locant
    rdkit_to_locant: dict[int, int] = {}
    for base_idx, rdkit_idx in enumerate(match):
        loc = locant_map_def.get(base_idx)
        if loc is not None:
            rdkit_to_locant[rdkit_idx] = loc

    # graph atom idx → rdkit atom idx (by atom order; same since no H reordering)
    # Non-ring neighbors of ring atoms are substituents
    from .molecule_analyzer import get_atom
    from .substituent import name_substituent

    ring_set = set(all_ring_atoms)
    substituents: list[tuple[int, str]] = []

    # Phase 148: 外環主官能基（カルボン酸・アルデヒド・アミド・ニトリル）の検出
    _EXOCYCLIC_SUFFIX: dict[str, str] = {
        "carboxylic_acid": "carboxylic acid",
        "aldehyde": "carbaldehyde",
        "amide": "carboxamide",
        "nitrile": "carbonitrile",
        "aldoxime": "carbaldehyde oxime",
        "aldhydrazone": "carbaldehyde hydrazone",
    }
    from .functional_group import detect_groups, principal_group as _pgrp_fn
    fg_groups = detect_groups(graph)
    pgrp_ex = _pgrp_fn(fg_groups)
    exo_loc: int | None = None
    exo_suffix: str | None = None
    exo_anchor_c: int | None = None

    if pgrp_ex is not None and pgrp_ex.group_type in _EXOCYCLIC_SUFFIX:
        anchor_c = pgrp_ex.atom_indices[0]
        for nb_idx in graph.adjacency[anchor_c]:
            if nb_idx in ring_set:
                loc_candidate = rdkit_to_locant.get(nb_idx)
                if loc_candidate is not None:
                    exo_loc = loc_candidate
                    exo_suffix = _EXOCYCLIC_SUFFIX[pgrp_ex.group_type]
                    exo_anchor_c = anchor_c
                    break

    for base_idx, rdkit_idx in enumerate(match):
        loc = locant_map_def.get(base_idx)
        if loc is None:
            continue  # junction atom
        # Find non-ring neighbors of this ring atom in graph
        for nb_idx in graph.adjacency[rdkit_idx]:
            if nb_idx in ring_set:
                continue
            if nb_idx in match_set:
                continue  # part of base structure (e.g., exo C=O oxygen in retained name)
            nb = get_atom(graph, nb_idx)
            if nb.symbol == "H":
                continue
            if nb_idx == exo_anchor_c:
                continue  # handled as principal group suffix
            sub_name = name_substituent(graph, nb_idx, ring_set)
            substituents.append((loc, sub_name))

    if exo_suffix is not None and exo_loc is not None:
        result = f"{base_name}-{exo_loc}-{exo_suffix}"
        if substituents:
            result = _apply_hetero_suffixes(result, substituents)
        # Amide N-substituents (N-methyl, N,N-dimethyl, etc.)
        if pgrp_ex is not None and pgrp_ex.group_type == "amide":
            from .substituent import _name_carbon_substituent as _ncs_fh
            from .constants import MULTIPLIER as _MULT_fh
            from collections import Counter as _Ctr_fh
            _n_fh = next(
                (ai for ai in pgrp_ex.atom_indices if get_atom(graph, ai).symbol == "N"), None
            )
            if _n_fh is not None and exo_anchor_c is not None:
                _cn_fh = [nb for nb in graph.adjacency[_n_fh]
                          if nb != exo_anchor_c and get_atom(graph, nb).symbol == "C"]
                if _cn_fh:
                    _sn_fh = [_ncs_fh(graph, c, {_n_fh}) for c in _cn_fh]
                    _cnt_fh = _Ctr_fh(_sn_fh)
                    _pp_fh: list[str] = []
                    for _s in sorted(_cnt_fh):
                        _c = _cnt_fh[_s]
                        _ss = f"({_s})" if _s.startswith("(") else _s
                        if _c == 1:
                            _pp_fh.append(f"N-{_ss}")
                        else:
                            _pp_fh.append(f"N,N-{_MULT_fh.get(_c, str(_c))}{_ss}")
                    _sep_fh = "-" if result and result[0].isdigit() else ""
                    result = "-".join(_pp_fh) + _sep_fh + result
        return result

    if not substituents:
        return base_name
    return _apply_hetero_suffixes(base_name, substituents)


def _is_fused_hetero_aromatic(graph: "MoleculeGraph") -> bool:
    """
    複数の芳香族環が縮合（2原子以上共有）しており、かつヘテロ原子を含む系かどうかを判定する。
    ビフェニル型（単結合で連結、非縮合）は False を返す。
    テーブル未収録の縮合ヘテロ芳香族に対して誤命名を防ぐ用途。
    """
    aromatic_rings = [r for r in graph.ring_atom_sets if _is_aromatic_ring(list(r), graph)]
    if len(aromatic_rings) < 2:
        return False
    # 縮合チェック: 任意の2環が 2原子以上を共有しているか
    has_fusion = False
    for i in range(len(aromatic_rings)):
        for j in range(i + 1, len(aromatic_rings)):
            if len(set(aromatic_rings[i]) & set(aromatic_rings[j])) >= 2:
                has_fusion = True
                break
        if has_fusion:
            break
    if not has_fusion:
        return False
    # 少なくとも1つのヘテロ原子が芳香族環に含まれるか
    from .molecule_analyzer import get_atom
    ring_atoms = {idx for ring in aromatic_rings for idx in ring}
    return any(
        get_atom(graph, idx).symbol not in ("C", "H") and get_atom(graph, idx).is_aromatic
        for idx in ring_atoms
    )


def name_heterocycle(graph: "MoleculeGraph") -> str | None:
    """
    分子がヘテロ環保留名に一致する場合、その IUPAC 名を返す。
    一致しない場合は None を返す。
    """
    from .ring_handler import find_hetero_rings

    # Phase 17: 縮合ヘテロ芳香族保留名を先行チェック
    fused_name = _try_fused_hetero_retained(graph)
    if fused_name is not None:
        return fused_name

    # 未収録の縮合ヘテロ芳香族は誤命名を防ぐため None を返す
    if _is_fused_hetero_aromatic(graph):
        return None

    hetero_rings = find_hetero_rings(graph)
    if not hetero_rings:
        return None

    # 単環のみ対応
    if len(hetero_rings) != 1:
        return None

    ring = hetero_rings[0]

    # 1. Phase 272: 部分不飽和環を最初にチェック (環内二重結合あり)
    # 飽和環や芳香族環は db_pairs が空 or is_aromatic で None を返すので安全
    pu_name = _match_partial_unsat(ring, graph)
    if pu_name is not None:
        rotation = _find_best_start(ring, graph)
        base_name, is_nh = pu_name, False
    else:
        # 1b. Phase 508: 2ヘテロ原子5員部分不飽和環 (thiazoline, oxazoline 等)
        pu2_name = _match_partial_unsat_2het(ring, graph)
        if pu2_name is not None:
            rotation = _find_best_start(ring, graph)
            base_name, is_nh = pu2_name, False
        else:
            # 2. 保留名テーブルを試みる
            match = _match_retained(ring, graph)
            if match is not None:
                base_name, is_nh, rotation = match
            else:
                # 3. Hantzsch-Widman 名を試みる (単一ヘテロ原子)
                hw_name = _match_hantzsch_widman(ring, graph)
                if hw_name is not None:
                    rotation = _find_best_start(ring, graph)
                    base_name, is_nh = hw_name, False
                else:
                    # 4. 多ヘテロ原子 HW 名 (7-10員環, a-命名法)
                    multi_name = _match_multi_het_ring(ring, graph)
                    if multi_name is None:
                        return None
                    rotation = _find_best_start(ring, graph)
                    base_name, is_nh = multi_name, False

    base_name, is_nh, rotation = (base_name, is_nh, rotation)
    locant_map = _build_locant_map(rotation)

    # 完全ベース名
    nh_prefix = "1H-" if is_nh else ""
    full_base = f"{nh_prefix}{base_name}"

    # Phase 22: 外環主官能基（カルボン酸・アルデヒド・アミド・ニトリル）の検出
    _EXOCYCLIC_SUFFIX: dict[str, str] = {
        "carboxylic_acid": "carboxylic acid",
        "aldehyde": "carbaldehyde",
        "amide": "carboxamide",
        "nitrile": "carbonitrile",
        "aldoxime": "carbaldehyde oxime",
        "aldhydrazone": "carbaldehyde hydrazone",
    }
    from .functional_group import detect_groups, principal_group as _principal_group
    fg_groups = detect_groups(graph)
    pgrp_ex = _principal_group(fg_groups)
    ring_set_ex = set(ring)

    if pgrp_ex is not None and pgrp_ex.group_type in _EXOCYCLIC_SUFFIX:
        anchor_c = pgrp_ex.atom_indices[0]
        # 主官能基アンカーがこのヘテロ芳香環に隣接していない場合 → _name_cyclic に委譲
        if not any(nb in ring_set_ex for nb in graph.adjacency[anchor_c]):
            return None
        for nb_idx in graph.adjacency[anchor_c]:
            if nb_idx in ring_set_ex:
                L = locant_map[nb_idx]
                n = len(rotation)
                L_rev = n + 2 - L
                if L_rev < L:
                    rev_rotation = [rotation[0]] + list(reversed(rotation[1:]))
                    locant_map = _build_locant_map(rev_rotation)
                    L = L_rev
                fg_suffix = _EXOCYCLIC_SUFFIX[pgrp_ex.group_type]
                other_subs = _collect_hetero_substituents(
                    graph, ring, locant_map, excluded_atoms={anchor_c}
                )
                result = f"{full_base}-{L}-{fg_suffix}"
                if other_subs:
                    result = _format_substituents(result, other_subs)
                # Amide N-substituents (N-methyl, N,N-dimethyl, etc.)
                if pgrp_ex.group_type == "amide":
                    from .molecule_analyzer import get_atom as _ga_am
                    from .substituent import _name_carbon_substituent as _ncs_am
                    from .constants import MULTIPLIER as _MULT_am
                    from collections import Counter as _Ctr_am
                    _n_am = next(
                        (ai for ai in pgrp_ex.atom_indices
                         if _ga_am(graph, ai).symbol == "N"), None
                    )
                    if _n_am is not None:
                        _cn_am = [nb for nb in graph.adjacency[_n_am]
                                  if nb != anchor_c and _ga_am(graph, nb).symbol == "C"]
                        _on_am = [nb for nb in graph.adjacency[_n_am]
                                  if _ga_am(graph, nb).symbol == "O"]
                        _nh_am = ["hydroxy" for _o in _on_am
                                  if any(_ga_am(graph, nb2).symbol == "H"
                                         for nb2 in graph.adjacency[_o])]
                        if _cn_am or _nh_am:
                            _sn_am = [_ncs_am(graph, c, {_n_am}) for c in _cn_am] + _nh_am
                            _cnt_am = _Ctr_am(_sn_am)
                            _pp_am: list[str] = []
                            for _s in sorted(_cnt_am):
                                _c = _cnt_am[_s]
                                _ss = f"({_s})" if _s.startswith("(") else _s
                                if _c == 1:
                                    _pp_am.append(f"N-{_ss}")
                                else:
                                    _pp_am.append(
                                        f"N,N-{_MULT_am.get(_c, str(_c))}{_ss}"
                                    )
                            _sep_am = "-" if result and result[0].isdigit() else ""
                            result = "-".join(_pp_am) + _sep_am + result
                return result

    # Phase 25/30: ラクタム (N 含有) / ラクトン (O 含有) 検出
    # 環内 C=O (exocyclic =O) を持つ環 → "{base}-{loc}-one"
    ring_set_lactam = set(ring)
    from .molecule_analyzer import get_bond_order as _get_bo_lact, get_atom as _get_atom_lact

    def _has_exo_dbl_o(c_idx: int) -> int | None:
        """環外に =O を持つ C のその O インデックスを返す。なければ None。"""
        for nb in graph.adjacency[c_idx]:
            if nb in ring_set_lactam:
                continue
            if (_get_atom_lact(graph, nb).symbol == "O"
                    and _get_bo_lact(graph, c_idx, nb) == 2.0):
                return nb
        return None

    # Phase 152/153: 複数 exo C=O グループを持つ環 → "-dione"/"-trione"
    # Phase 58/127 より前に実行して 3+ C=O (trione) も正しく処理する
    exo_co_ring_cs_e = []
    exo_co_oxygens_e = []
    for ring_idx_e in ring:
        if _get_atom_lact(graph, ring_idx_e).symbol != "C":
            continue
        exo_o_e = _has_exo_dbl_o(ring_idx_e)
        if exo_o_e is not None:
            exo_co_ring_cs_e.append(ring_idx_e)
            exo_co_oxygens_e.append(exo_o_e)
    if len(exo_co_ring_cs_e) >= 2:
        locs_e = [locant_map[c] for c in exo_co_ring_cs_e]
        locs_sorted_e = sorted(locs_e)
        suffix_e = "dione" if len(locs_sorted_e) == 2 else "trione"
        # Phase 403: pre-check for pyrrole-imide pattern (maleimide/N-subst. maleimide)
        _is_pyrrole_imide_403 = False
        _n_403: "int | None" = None
        if base_name.endswith("pyrrole"):
            _n_403 = next((idx for idx in ring if _get_atom_lact(graph, idx).symbol == "N"), None)
            if _n_403 is not None:
                _nc_403 = [nb for nb in graph.adjacency[_n_403] if nb in ring_set_lactam]
                if (len(_nc_403) == 2 and _has_exo_dbl_o(_nc_403[0]) is not None
                        and _has_exo_dbl_o(_nc_403[1]) is not None):
                    _is_pyrrole_imide_403 = True
        # Phase 405: pre-check for phthalimide / isoindole-dione pattern
        _is_isoindole_dione_405 = False
        _n_405: "int | None" = None
        if base_name == "pyrrolidine" and len(locs_sorted_e) == 2:
            _arom_junc_405 = [
                idx for idx in ring
                if _get_atom_lact(graph, idx).is_aromatic
                and any(idx in frozenset(r) for r in graph.ring_atom_sets
                        if frozenset(r) != ring_set_lactam)
            ]
            if len(_arom_junc_405) == 2:
                _n_405 = next(
                    (idx for idx in ring if _get_atom_lact(graph, idx).symbol == "N"), None
                )
                if _n_405 is not None:
                    _is_isoindole_dione_405 = True
        if is_nh:
            # Phase 401: indicated-H dione — use canonical rotation (already IUPAC-optimal)
            _nh_atoms_e = [
                idx for idx in ring
                if _get_atom_lact(graph, idx).symbol == "N"
                and any(_get_atom_lact(graph, nb).symbol == "H"
                        for nb in graph.adjacency[idx])
            ]
            _nh_locs_e = sorted(locant_map[nh] for nh in _nh_atoms_e)
            _nh_str_e = ",".join(f"{l}H" for l in _nh_locs_e)
            loc_str_e = ",".join(str(l) for l in locs_sorted_e)
            dione_name_e = f"{base_name}-{loc_str_e}({_nh_str_e})-{suffix_e}"
        elif _is_pyrrole_imide_403:
            # Phase 403: maleimide / N-substituted maleimide → 1H-pyrrole-{loc}-dione
            loc_str_e = ",".join(str(l) for l in locs_sorted_e)
            dione_name_e = f"1H-pyrrole-{loc_str_e}-{suffix_e}"
            _excl_403 = set(exo_co_oxygens_e)
            _subs_raw_403 = _collect_hetero_substituents(
                graph, ring, locant_map, excluded_atoms=_excl_403
            )
            _n_num_403 = locant_map[_n_403] if _n_403 is not None else 1
            _subs_403 = [(_n_num_403 if loc == "N" else loc, nm) for loc, nm in _subs_raw_403]
            if not _subs_403:
                return dione_name_e
            return _format_substituents(dione_name_e, _subs_403)
        elif base_name == "hexahydropyrimidine" and len(locs_sorted_e) >= 3:
            # Phase 402: barbituric acid → pyrimidine parent + indicated H (trione only)
            _ih_locs_e = sorted(
                locant_map[idx] for idx in ring
                if any(_get_atom_lact(graph, nb).symbol == "H"
                       for nb in graph.adjacency[idx])
            )
            _ih_str_e = ",".join(f"{l}H" for l in _ih_locs_e)
            loc_str_e = ",".join(str(l) for l in locs_sorted_e)
            dione_name_e = f"pyrimidine-{loc_str_e}({_ih_str_e})-{suffix_e}"
        elif _is_isoindole_dione_405:
            # Phase 405: phthalimide / N-subst. phthalimide → isoindole-1,3(2H)-dione
            dione_name_e = "isoindole-1,3(2H)-dione"
            _excl_405 = set(exo_co_oxygens_e)
            _subs_raw_405 = _collect_hetero_substituents(
                graph, ring, locant_map, excluded_atoms=_excl_405
            )
            _subs_405 = [(2 if loc == "N" else loc, nm) for loc, nm in _subs_raw_405]
            if not _subs_405:
                return dione_name_e
            return _format_substituents(dione_name_e, _subs_405)
        else:
            rev_rotation_e = [rotation[0]] + list(reversed(rotation[1:]))
            locant_map_rev_e = _build_locant_map(rev_rotation_e)
            locs_rev_e = [locant_map_rev_e[c] for c in exo_co_ring_cs_e]
            if sorted(locs_rev_e) < sorted(locs_e):
                locant_map = locant_map_rev_e
                locs_e = locs_rev_e
            locs_sorted_e = sorted(locs_e)
            loc_str_e = ",".join(str(l) for l in locs_sorted_e)
            _eb = full_base  # -dione/-trione starts with consonant: no elision
            dione_name_e = f"{_eb}-{loc_str_e}-{suffix_e}"
        excl_e = set(exo_co_oxygens_e)
        other_subs_e = _collect_hetero_substituents(
            graph, ring, locant_map, excluded_atoms=excl_e
        )
        if not other_subs_e:
            return dione_name_e
        return _format_substituents(dione_name_e, other_subs_e)

    # Phase 58: 環状イミド検出 (N が両側の C(=O) に挟まれている)
    # O=C1NC(=O)CC1 → pyrrolidine-2,5-dione (exo C=O 1 個の場合のフォールバック)
    for n_idx_im in ring:
        if _get_atom_lact(graph, n_idx_im).symbol != "N":
            continue
        n_ring_c = [nb for nb in graph.adjacency[n_idx_im]
                    if nb in ring_set_lactam and _get_atom_lact(graph, nb).symbol == "C"]
        if len(n_ring_c) == 2:
            o1 = _has_exo_dbl_o(n_ring_c[0])
            o2 = _has_exo_dbl_o(n_ring_c[1])
            if o1 is not None and o2 is not None:
                loc1 = locant_map[n_ring_c[0]]
                loc2 = locant_map[n_ring_c[1]]
                # ロカントが小さい方から順に並べる
                locs = sorted([loc1, loc2])
                # より小さいロカント集合になる向きを選ぶ (ロカント最小化)
                n_ring = len(ring)
                locs_rev = sorted([n_ring + 2 - loc1, n_ring + 2 - loc2])
                if locs_rev < locs:
                    rev_rotation = [rotation[0]] + list(reversed(rotation[1:]))
                    locant_map = _build_locant_map(rev_rotation)
                    locs = locs_rev
                loc_str = ",".join(str(l) for l in locs)
                _imide_base = full_base  # -dione starts with consonant: no elision
                imide_name = f"{_imide_base}-{loc_str}-dione"
                excl_o = {o1, o2}
                other_subs = _collect_hetero_substituents(
                    graph, ring, locant_map, excluded_atoms=excl_o
                )
                if not other_subs:
                    return imide_name
                return _format_substituents(imide_name, other_subs)

    # Phase 127: 環状無水物 (oxa-dione): 環内 O の両隣 C が共に exo C=O → "{base}-{loc1},{loc2}-dione"
    for o_idx_ca in ring:
        if _get_atom_lact(graph, o_idx_ca).symbol != "O":
            continue
        o_ring_c = [nb for nb in graph.adjacency[o_idx_ca]
                    if nb in ring_set_lactam and _get_atom_lact(graph, nb).symbol == "C"]
        if len(o_ring_c) == 2:
            exo1 = _has_exo_dbl_o(o_ring_c[0])
            exo2 = _has_exo_dbl_o(o_ring_c[1])
            if exo1 is not None and exo2 is not None:
                loc1 = locant_map[o_ring_c[0]]
                loc2 = locant_map[o_ring_c[1]]
                locs = sorted([loc1, loc2])
                loc_str = ",".join(str(l) for l in locs)
                _anhy_base = full_base  # -dione starts with consonant: no elision
                dione_name = f"{_anhy_base}-{loc_str}-dione"
                excl_o = {exo1, exo2}
                other_subs = _collect_hetero_substituents(
                    graph, ring, locant_map, excluded_atoms=excl_o
                )
                if not other_subs:
                    return dione_name
                return _format_substituents(dione_name, other_subs)

    # Phase 152: 複数 exo C=O グループを持つ環 → "-dione"/"-trione" (例: piperazine-2,5-dione)
    exo_co_ring_cs = []
    exo_co_oxygens = []
    for ring_idx_m in ring:
        if _get_atom_lact(graph, ring_idx_m).symbol != "C":
            continue
        exo_o = _has_exo_dbl_o(ring_idx_m)
        if exo_o is not None:
            exo_co_ring_cs.append(ring_idx_m)
            exo_co_oxygens.append(exo_o)
    if len(exo_co_ring_cs) >= 2:
        # 複数 C=O → dione/trione suffix
        locs_m = [locant_map[c] for c in exo_co_ring_cs]
        n_ring_m = len(ring)
        # ロカント集合を最小化: すべての C=O ロカントの和が最小の向きを選ぶ
        rev_rotation_m = [rotation[0]] + list(reversed(rotation[1:]))
        locant_map_rev_m = _build_locant_map(rev_rotation_m)
        locs_rev_m = [locant_map_rev_m[c] for c in exo_co_ring_cs]
        if sorted(locs_rev_m) < sorted(locs_m):
            locant_map = locant_map_rev_m
            locs_m = locs_rev_m
        locs_sorted = sorted(locs_m)
        loc_str_m = ",".join(str(l) for l in locs_sorted)
        suffix_m = "dione" if len(locs_sorted) == 2 else "trione"
        _dm = full_base  # -dione/-trione starts with consonant: no elision
        dione_name_m = f"{_dm}-{loc_str_m}-{suffix_m}"
        excl_m = set(exo_co_oxygens)
        other_subs_m = _collect_hetero_substituents(
            graph, ring, locant_map, excluded_atoms=excl_m
        )
        if not other_subs_m:
            return dione_name_m
        return _format_substituents(dione_name_m, other_subs_m)

    has_ring_hetero = any(_get_atom_lact(graph, idx).symbol not in ("C", "H") for idx in ring)
    if has_ring_hetero:
        for ring_idx in ring:
            ring_atom = _get_atom_lact(graph, ring_idx)
            if ring_atom.symbol != "C":
                continue
            for nb_idx in graph.adjacency[ring_idx]:
                if nb_idx in ring_set_lactam:
                    continue
                nb = _get_atom_lact(graph, nb_idx)
                if nb.symbol == "O" and _get_bo_lact(graph, ring_idx, nb_idx) == 2.0:
                    loc = locant_map[ring_idx]
                    n_ring = len(ring)
                    loc_rev = n_ring + 2 - loc
                    if loc_rev < loc:
                        rev_rotation = [rotation[0]] + list(reversed(rotation[1:]))
                        locant_map = _build_locant_map(rev_rotation)
                        loc = locant_map[ring_idx]
                    if is_nh:
                        # Phase 400: unified IUPAC rotation — O < S < N priority
                        _all_h = [idx for idx in ring
                                  if _get_atom_lact(graph, idx).symbol in ("N", "O", "S")]
                        _nh_atom = next(
                            (idx for idx in ring
                             if _get_atom_lact(graph, idx).symbol == "N"
                             and any(_get_atom_lact(graph, nb).symbol == "H"
                                     for nb in graph.adjacency[idx])),
                            None,
                        )
                        _nr = len(ring)
                        _inf = _nr + 1
                        _bk: tuple = ([_inf] * _nr, [_inf] * _nr, [_inf] * _nr, _inf, _inf)
                        _blm: dict[int, int] = locant_map
                        _bloc: int = loc
                        _bnh: int = 1
                        for _hstart in _all_h:
                            _pos = ring.index(_hstart)
                            for _dd in (1, -1):
                                _rot = [ring[(_pos + _dd * i) % _nr] for i in range(_nr)]
                                _lm2 = {_rot[i]: i + 1 for i in range(_nr)}
                                _ok = sorted(_lm2[h] for h in _all_h
                                             if _get_atom_lact(graph, h).symbol == "O")
                                _sk = sorted(_lm2[h] for h in _all_h
                                             if _get_atom_lact(graph, h).symbol == "S")
                                _nk = sorted(_lm2[h] for h in _all_h
                                             if _get_atom_lact(graph, h).symbol == "N")
                                _co2 = _lm2[ring_idx]
                                _nh2 = _lm2[_nh_atom] if _nh_atom is not None else _inf
                                _k: tuple = (_ok, _sk, _nk, _co2, _nh2)
                                if _k < _bk:
                                    _bk = _k
                                    _blm = _lm2
                                    _bloc = _co2
                                    _bnh = _nh2
                        locant_map = _blm
                        loc = _bloc
                        _nh_loc = _bnh
                        _base_nh = (base_name[:-1] if base_name.endswith("e") else base_name)
                        lactam_name = f"{_base_nh}-{loc}({_nh_loc}H)-one"
                    else:
                        lactam_base = (full_base[:-1] if full_base.endswith("e") else full_base)
                        lactam_name = f"{lactam_base}-{loc}-one"
                    other_subs = _collect_hetero_substituents(
                        graph, ring, locant_map, excluded_atoms={nb_idx}
                    )
                    if not other_subs:
                        return lactam_name
                    return _format_substituents(lactam_name, other_subs)
                # Phase 391: thiolactam – exocyclic =S on a ring C in an N-heterocycle
                if nb.symbol == "S" and _get_bo_lact(graph, ring_idx, nb_idx) == 2.0:
                    loc = locant_map[ring_idx]
                    n_ring = len(ring)
                    loc_rev = n_ring + 2 - loc
                    if loc_rev < loc:
                        rev_rotation = [rotation[0]] + list(reversed(rotation[1:]))
                        locant_map = _build_locant_map(rev_rotation)
                        loc = locant_map[ring_idx]
                    if is_nh:
                        # Phase 400: unified IUPAC rotation — O < S < N priority
                        _all_h_tl = [idx for idx in ring
                                     if _get_atom_lact(graph, idx).symbol in ("N", "O", "S")]
                        _nh_atom_tl = next(
                            (idx for idx in ring
                             if _get_atom_lact(graph, idx).symbol == "N"
                             and any(_get_atom_lact(graph, nb2).symbol == "H"
                                     for nb2 in graph.adjacency[idx])),
                            None,
                        )
                        _nr_tl = len(ring)
                        _inf_tl = _nr_tl + 1
                        _bk_tl: tuple = (
                            [_inf_tl] * _nr_tl, [_inf_tl] * _nr_tl, [_inf_tl] * _nr_tl,
                            _inf_tl, _inf_tl,
                        )
                        _blm_tl: dict[int, int] = locant_map
                        _bloc_tl: int = loc
                        _bnh_tl: int = 1
                        for _hstart_tl in _all_h_tl:
                            _pos_tl = ring.index(_hstart_tl)
                            for _dd_tl in (1, -1):
                                _rot_tl = [ring[(_pos_tl + _dd_tl * i) % _nr_tl]
                                           for i in range(_nr_tl)]
                                _lm_tl = {_rot_tl[i]: i + 1 for i in range(_nr_tl)}
                                _ok_tl = sorted(_lm_tl[h] for h in _all_h_tl
                                                if _get_atom_lact(graph, h).symbol == "O")
                                _sk_tl = sorted(_lm_tl[h] for h in _all_h_tl
                                                if _get_atom_lact(graph, h).symbol == "S")
                                _nk_tl = sorted(_lm_tl[h] for h in _all_h_tl
                                                if _get_atom_lact(graph, h).symbol == "N")
                                _co_tl = _lm_tl[ring_idx]
                                _nh_tl2 = (_lm_tl[_nh_atom_tl]
                                           if _nh_atom_tl is not None else _inf_tl)
                                _k_tl: tuple = (_ok_tl, _sk_tl, _nk_tl, _co_tl, _nh_tl2)
                                if _k_tl < _bk_tl:
                                    _bk_tl = _k_tl
                                    _blm_tl = _lm_tl
                                    _bloc_tl = _co_tl
                                    _bnh_tl = _nh_tl2
                        locant_map = _blm_tl
                        loc = _bloc_tl
                        _nh_loc_tl = _bnh_tl
                        _base_nh_tl = (base_name[:-1]
                                       if base_name.endswith("e") else base_name)
                        thiolactam_name = f"{_base_nh_tl}-{loc}({_nh_loc_tl}H)-thione"
                    else:
                        thiolactam_name = f"{full_base}-{loc}-thione"
                    other_subs = _collect_hetero_substituents(
                        graph, ring, locant_map, excluded_atoms={nb_idx}
                    )
                    if not other_subs:
                        return thiolactam_name
                    return _format_substituents(thiolactam_name, other_subs)

    # 置換基収集: 全有効回転候補を試して最小ロカント集合 → 同点時アルファベット順
    n_ring = len(rotation)
    fwd_hetero = sorted(
        i + 1 for i, a in enumerate(rotation) if _atom_sig(graph, a) != "C"
    )

    # 各ヘテロ原子を起点として順/逆方向回転を生成し、ヘテロ原子ロカント集合が同じものを収集
    rot_candidates = [rotation]
    seen_rots: set[tuple[int, ...]] = {tuple(rotation)}
    for i in range(n_ring):
        if _atom_sig(graph, rotation[i]) == "C":
            continue
        for sign in (1, -1):
            alt = [rotation[(i + sign * j) % n_ring] for j in range(n_ring)]
            alt_key = tuple(alt)
            if alt_key in seen_rots:
                continue
            seen_rots.add(alt_key)
            alt_hetero = sorted(
                k + 1 for k, a in enumerate(alt) if _atom_sig(graph, a) != "C"
            )
            # Only swap to an alternative that starts with the same element (e.g. N↔N ok,
            # N↔O not allowed: in isoxazole O must remain at locant 1, not N).
            if alt_hetero == fwd_hetero and (
                _get_atom_lact(graph, alt[0]).symbol == _get_atom_lact(graph, rotation[0]).symbol
            ):
                rot_candidates.append(alt)

    # 各候補の置換基ロカント集合を収集
    candidate_subs = [
        _collect_hetero_substituents(graph, ring, _build_locant_map(rot))
        for rot in rot_candidates
    ]

    # 最小ロカント集合を選択し、同点時はアルファベット順で tie-break
    def _loc_sort(loc: "int | str") -> tuple:
        return (0, loc, "") if isinstance(loc, int) else (1, 0, loc)

    def _subs_key(subs: list[tuple["int | str", str]]) -> tuple:
        ordered = sorted(subs, key=lambda x: _loc_sort(x[0]))
        locs = sorted(_loc_sort(loc) for loc, _ in subs)
        return (locs, [nm for _, nm in ordered])

    candidate_subs.sort(key=_subs_key)
    substituents = candidate_subs[0]

    if not substituents:
        return full_base

    return _apply_hetero_suffixes(full_base, substituents)
