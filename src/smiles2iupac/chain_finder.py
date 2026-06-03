"""
主鎖の探索とロカント番号付け。

IUPAC 2013 Blue Book P-44 の選択規則を順に適用し、
最適な主鎖（PrincipalChain）を決定する。
"""

from __future__ import annotations

from dataclasses import dataclass

from .constants import FUNCTIONAL_GROUPS
from .functional_group import FunctionalGroup
from .molecule_analyzer import MoleculeGraph, get_atom, get_bond_order, carbon_indices


@dataclass
class PrincipalChain:
    atom_indices: list[int]        # ロカント1から始まる順序付き炭素インデックス
    locant_map: dict[int, int]     # atom_idx -> ロカント番号 (1始まり)

    @property
    def length(self) -> int:
        return len(self.atom_indices)


def find_principal_chain(
    graph: MoleculeGraph,
    principal_grp: FunctionalGroup | None,
) -> PrincipalChain:
    """
    主鎖を選択してロカントを割り当てる。

    選択規則 (IUPAC P-44):
      1. principal characteristic group の炭素を含む最長鎖
      2. 最長炭素鎖
      3. 多重結合 (C=C, C≡C) の数が最多
      4. 置換基数が最多
      5. ロカント集合が最小

    Args:
        graph: MoleculeGraph
        principal_grp: detect_groups()[0] または None (alkane)

    Returns:
        PrincipalChain
    """
    # 環原子が存在する場合は非環炭素のみを主鎖に使用
    # (フェニルアルキル化合物等: 環は substituent として扱う)
    from .molecule_analyzer import non_ring_carbon_indices
    has_ring_c = any(a.symbol == "C" and a.in_ring for a in graph.atoms)
    if has_ring_c:
        c_idxs = non_ring_carbon_indices(graph)
        if not c_idxs:
            raise NotImplementedError("Ring compounds are not supported in Phase 1-2.")
    else:
        # 全炭素インデックス
        c_idxs = carbon_indices(graph)
    if not c_idxs:
        raise ValueError("No carbon atoms found in the molecule.")

    # principal group に属する炭素インデックスのセット
    required_carbons: set[int] = set()
    if principal_grp is not None:
        for ai in principal_grp.atom_indices:
            if get_atom(graph, ai).symbol == "C":
                required_carbons.add(ai)

    # 全炭素単純経路を DFS で列挙（炭素のみ通る）
    all_paths = _enumerate_carbon_paths(graph, c_idxs)

    # required_carbons を含む経路のみに絞る
    if required_carbons:
        filtered = [p for p in all_paths if required_carbons.issubset(set(p))]
        if not filtered:
            # required を含む経路がない場合は全経路を使う（エラー回避）
            filtered = all_paths
        all_paths = filtered

    if not all_paths:
        # 炭素が1つしかない場合
        return _make_chain(graph, c_idxs[:1], principal_grp)

    # 選択規則でソート（全規則を複合キーで）
    best_path = _select_best_path(graph, all_paths)

    # ロカント方向を決定
    oriented = _orient_chain(graph, best_path, principal_grp)

    return _make_chain(graph, oriented, principal_grp)


# ─── 経路列挙 ──────────────────────────────────────────────────────────

def _enumerate_carbon_paths(
    graph: MoleculeGraph, c_idxs: list[int]
) -> list[list[int]]:
    """
    炭素原子のみを通る全単純経路をDFSで列挙する。
    端点ペア (start, end) ごとに最長経路を1つ返す形ではなく、
    全ての単純経路を返す（長さ1の単一炭素経路も含む）。
    """
    c_set = set(c_idxs)
    all_paths: list[list[int]] = []

    # 炭素同士の隣接関係
    c_adj: dict[int, list[int]] = {
        c: [nb for nb in graph.adjacency[c] if nb in c_set]
        for c in c_idxs
    }

    def dfs(current: int, visited: set[int], path: list[int]) -> None:
        path.append(current)
        visited.add(current)
        has_extension = False
        for nb in c_adj[current]:
            if nb not in visited:
                has_extension = True
                dfs(nb, visited, path)
        if not has_extension or len(c_adj[current]) == 0:
            # 末端ノードに到達 → 経路を記録
            all_paths.append(list(path))
        path.pop()
        visited.remove(current)

    # 各炭素を始点として DFS
    # 重複を避けるため (start < end) の経路のみ収集後に逆向きを除く
    seen_pairs: set[tuple[int, int]] = set()
    for start in c_idxs:
        # DFS で全末端経路を列挙
        def dfs_collect(cur: int, vis: set[int], path: list[int]) -> None:
            path.append(cur)
            vis.add(cur)
            extended = False
            for nb in c_adj[cur]:
                if nb not in vis:
                    extended = True
                    dfs_collect(nb, vis, path)
            if not extended:
                # 末端に到達
                if len(path) > 0:
                    pair = (min(path[0], path[-1]), max(path[0], path[-1]))
                    if pair not in seen_pairs:
                        seen_pairs.add(pair)
                        all_paths.append(list(path))
            path.pop()
            vis.remove(cur)

        dfs_collect(start, set(), [])

    # 単一炭素（孤立）の場合
    if not all_paths:
        all_paths = [[c] for c in c_idxs]

    return all_paths


# ─── 選択規則 ─────────────────────────────────────────────────────────

def _select_best_path(graph: MoleculeGraph, paths: list[list[int]]) -> list[int]:
    """
    IUPAC P-44 の優先規則で最良経路を選ぶ。
    """

    def sort_key(path: list[int]):
        length = len(path)
        multiple_bonds = _count_multiple_bonds(graph, path)
        substituent_count = _count_substituents(graph, path)
        # ロカント集合: 小さいほどよい → 負にして降順ソートに対応
        locant_sets = _locant_set_for_sorting(graph, path)
        return (
            length,                     # 1) 最長優先 (大きいほど良い)
            multiple_bonds,             # 2) 多重結合数 (大きいほど良い)
            substituent_count,          # 3) 置換基数 (大きいほど良い)
            [-l for l in locant_sets],  # 4) ロカント集合最小 (負にして大きいほど良い)
        )

    return max(paths, key=sort_key)


def _count_multiple_bonds(graph: MoleculeGraph, path: list[int]) -> int:
    """経路上の C=C, C≡C の数を返す。"""
    count = 0
    for i in range(len(path) - 1):
        bo = get_bond_order(graph, path[i], path[i + 1])
        if bo >= 2.0:
            count += 1
    return count


def _count_substituents(graph: MoleculeGraph, path: list[int]) -> int:
    """経路上の炭素に付く非主鎖・非H の置換基数を返す。"""
    path_set = set(path)
    count = 0
    for c_idx in path:
        for nb in graph.adjacency[c_idx]:
            nb_atom = get_atom(graph, nb)
            if nb not in path_set and nb_atom.symbol not in ("H",):
                count += 1
    return count


def _locant_set_for_sorting(graph: MoleculeGraph, path: list[int]) -> list[int]:
    """
    置換基・多重結合のロカント集合を昇順リストで返す。
    両方向で評価し、小さい方を選ぶ。
    """
    set_fwd = _compute_locant_set(graph, path)
    set_rev = _compute_locant_set(graph, list(reversed(path)))
    return min(set_fwd, set_rev)


def _compute_locant_set(graph: MoleculeGraph, path: list[int]) -> list[int]:
    """
    指定方向の経路でのロカントセット（置換基 + 多重結合）を返す。
    """
    path_set = set(path)
    locants = []
    for locant, c_idx in enumerate(path, start=1):
        # 置換基
        for nb in graph.adjacency[c_idx]:
            if nb not in path_set and get_atom(graph, nb).symbol != "H":
                locants.append(locant)
        # 多重結合（結合の小さい方のロカントを記録）
        if locant < len(path):
            next_idx = path[locant]  # locant は 1始まりなので path[locant] が次
            bo = get_bond_order(graph, c_idx, next_idx)
            if bo >= 2.0:
                locants.append(locant)
    return sorted(locants)


# ─── ロカント方向決定 ────────────────────────────────────────────────

def _orient_chain(
    graph: MoleculeGraph,
    path: list[int],
    principal_grp: FunctionalGroup | None,
) -> list[int]:
    """
    主鎖の向き（どちら端をロカント1にするか）を決定する。

    規則:
      1. carboxylic_acid / aldehyde: 官能基炭素を C1 (path[0]) にする
      2. alcohol / ketone: 官能基炭素のロカントが最小になる向き
      3. alkene / alkyne: 多重結合のロカント集合が最小になる向き
      4. alkane / 未決: 置換基のロカント集合が最小になる向き
    """
    if len(path) == 1:
        return path

    fwd = path
    rev = list(reversed(path))

    if principal_grp is None or principal_grp.group_type == "alkane":
        return _choose_by_substituent_locants(graph, fwd, rev)

    gtype = principal_grp.group_type
    spec = FUNCTIONAL_GROUPS.get(gtype)

    # alkene / alkyne: 多重結合ロカントを最小にする
    if gtype in ("alkene", "alkyne"):
        return _choose_by_multiple_bond_locants(graph, fwd, rev)

    if spec is None:
        return _choose_by_substituent_locants(graph, fwd, rev)

    grp_carbons = {
        ai for ai in principal_grp.atom_indices
        if get_atom(graph, ai).symbol == "C"
    }

    if spec.anchor_c1:
        # 官能基炭素を C1 に固定 (carboxylic_acid, aldehyde, nitrile)
        # dioic_acid / dial は両端が anchor → 多重結合ロカント→置換基ロカントで方向決定
        both_ends_are_anchors = path[0] in grp_carbons and path[-1] in grp_carbons
        if both_ends_are_anchors:
            def _mb_locs_anchor(p: list[int]) -> list[int]:
                locs = []
                for i in range(len(p) - 1):
                    if get_bond_order(graph, p[i], p[i + 1]) >= 2.0:
                        locs.append(i + 1)
                return sorted(locs)
            mb_fwd = _mb_locs_anchor(fwd)
            mb_rev = _mb_locs_anchor(rev)
            if mb_fwd < mb_rev:
                return fwd
            if mb_rev < mb_fwd:
                return rev
            return _choose_by_substituent_locants(graph, fwd, rev)
        if path[0] in grp_carbons:
            return fwd
        elif path[-1] in grp_carbons:
            return rev
        return fwd

    if spec.needs_locant:
        # 官能基炭素のロカントが最小になる向き (alcohol, ketone, amine)
        loc_fwd = min(
            (i + 1 for i, c in enumerate(fwd) if c in grp_carbons),
            default=len(fwd),
        )
        loc_rev = min(
            (i + 1 for i, c in enumerate(rev) if c in grp_carbons),
            default=len(rev),
        )
        if loc_fwd != loc_rev:
            return fwd if loc_fwd < loc_rev else rev
        # タイなら多重結合ロカントを次に試す (IUPAC P-44.1: C→D→E順)
        def _mb_locs(path: list[int]) -> list[int]:
            locs = []
            for i in range(len(path) - 1):
                if get_bond_order(graph, path[i], path[i + 1]) >= 2.0:
                    locs.append(i + 1)
            return sorted(locs)
        mb_fwd = _mb_locs(fwd)
        mb_rev = _mb_locs(rev)
        if mb_fwd != mb_rev:
            return fwd if mb_fwd < mb_rev else rev
        # それでもタイなら置換基ロカントで決める (Phase 189)
        return _choose_by_substituent_locants(graph, fwd, rev)

    return _choose_by_substituent_locants(graph, fwd, rev)


def _choose_by_substituent_locants(
    graph: MoleculeGraph, fwd: list[int], rev: list[int]
) -> list[int]:
    """置換基のロカント集合が最小になる向きを返す。"""
    path_set = set(fwd)

    def sub_locants(path: list[int]) -> list[int]:
        locs = []
        for i, c in enumerate(path, start=1):
            for nb in graph.adjacency[c]:
                if nb not in path_set and get_atom(graph, nb).symbol != "H":
                    locs.append(i)
        return sorted(locs)

    l_fwd = sub_locants(fwd)
    l_rev = sub_locants(rev)
    return fwd if l_fwd <= l_rev else rev


def _choose_by_multiple_bond_locants(
    graph: MoleculeGraph, fwd: list[int], rev: list[int]
) -> list[int]:
    """多重結合のロカント集合が最小になる向きを返す。

    IUPAC 2013 P-31.1.6.3: 組み合わせロカント集合が等しい場合は
    二重結合が低いロカントを得る向きを選ぶ。
    """

    def mb_locants(path: list[int]) -> list[int]:
        locs = []
        for i in range(len(path) - 1):
            bo = get_bond_order(graph, path[i], path[i + 1])
            if bo >= 2.0:
                locs.append(i + 1)
        return sorted(locs)

    def db_locants(path: list[int]) -> list[int]:
        locs = []
        for i in range(len(path) - 1):
            if get_bond_order(graph, path[i], path[i + 1]) == 2.0:
                locs.append(i + 1)
        return sorted(locs)

    l_fwd = mb_locants(fwd)
    l_rev = mb_locants(rev)
    if l_fwd != l_rev:
        return fwd if l_fwd < l_rev else rev
    # 全多重結合ロカント集合が等しい → 二重結合に低いロカントを与える方向を選ぶ
    db_fwd = db_locants(fwd)
    db_rev = db_locants(rev)
    if db_fwd != db_rev:
        return fwd if db_fwd < db_rev else rev
    # タイなら置換基で決める
    return _choose_by_substituent_locants(graph, fwd, rev)


# ─── PrincipalChain 生成 ────────────────────────────────────────────

def _make_chain(
    graph: MoleculeGraph,
    path: list[int],
    principal_grp: FunctionalGroup | None,
) -> PrincipalChain:
    locant_map = {atom_idx: loc for loc, atom_idx in enumerate(path, start=1)}
    return PrincipalChain(atom_indices=path, locant_map=locant_map)


def get_multiple_bond_locants(
    graph: MoleculeGraph, chain: PrincipalChain
) -> dict[str, list[int]]:
    """
    主鎖上の多重結合のロカントを返す。

    Returns:
        {'ene': [2, 4], 'yne': [6]} など
    """
    ene_locants: list[int] = []
    yne_locants: list[int] = []
    path = chain.atom_indices
    for i in range(len(path) - 1):
        bo = get_bond_order(graph, path[i], path[i + 1])
        locant = i + 1  # 低い方のロカント
        if bo == 2.0:
            ene_locants.append(locant)
        elif bo == 3.0:
            yne_locants.append(locant)
    return {"ene": sorted(ene_locants), "yne": sorted(yne_locants)}


# ─── group_namers で使用するチェーンユーティリティ ─────────────────────────

def collect_acid_chain(graph: MoleculeGraph, start_c: int, excluded_set, get_atom_fn) -> list[int]:
    """カルボニル C から excluded_set を除いた最長炭素鎖を返す (chain順)。

    芳香族環 C に到達したとき、直前の原子が非芳香族なら環への侵入を停止する。
    これにより PhCH₂-C(=O)OR の酸鎖が benzene 環を含まないようにする。
    """
    excl = set(excluded_set)
    visited: set[int] = set(excl)

    start_atom = get_atom_fn(graph, start_c)
    start_is_aromatic = start_atom.in_ring and start_atom.is_aromatic

    def _longest(idx: int, from_aromatic: bool) -> list[int]:
        visited.add(idx)
        atom = get_atom_fn(graph, idx)
        this_is_aromatic = atom.in_ring and atom.is_aromatic
        best_ext: list[int] = []
        for nb in graph.adjacency[idx]:
            if nb in visited or nb in excl:
                continue
            nb_atom = get_atom_fn(graph, nb)
            if nb_atom.symbol != "C":
                continue
            # 非芳香族領域から芳香族環に侵入しない
            nb_is_aromatic = nb_atom.in_ring and nb_atom.is_aromatic
            if nb_is_aromatic and not this_is_aromatic:
                continue
            ext = _longest(nb, this_is_aromatic)
            if len(ext) > len(best_ext):
                best_ext = ext
        visited.discard(idx)
        return [idx] + best_ext

    return _longest(start_c, start_is_aromatic)


def chain_through_pivot(
    graph: MoleculeGraph, pivot_c: int, excluded_set, get_atom_fn
) -> tuple[list[int], int]:
    """pivot_c を通る最長炭素鎖を見つけ (chain, locant) を返す。
    locant は 1始まりの pivot_c の位置。
    S, N などのヘテロ原子が pivot に付く場合に使用。
    """
    excl_for_arms = set(excluded_set) | {pivot_c}
    arms: list[list[int]] = []
    for nb in graph.adjacency[pivot_c]:
        if nb in excluded_set or get_atom_fn(graph, nb).symbol != "C":
            continue
        arm = collect_acid_chain(graph, nb, excl_for_arms, get_atom_fn)
        arms.append(arm)
    arms.sort(key=len)
    if not arms:
        return [pivot_c], 1
    if len(arms) == 1:
        return [pivot_c] + arms[0], 1
    left = list(reversed(arms[0]))
    right = arms[-1]
    chain = left + [pivot_c] + right
    locant = len(left) + 1
    return chain, locant


def chain_multiple_bonds(graph: MoleculeGraph, chain: list[int]) -> tuple[list[int], list[int]]:
    """chain 内の連続原子ペアの多重結合ロカントを返す。(ene_locs, yne_locs)"""
    ene: list[int] = []
    yne: list[int] = []
    for i in range(len(chain) - 1):
        bo = get_bond_order(graph, chain[i], chain[i + 1])
        if bo == 2.0:
            ene.append(i + 1)
        elif bo == 3.0:
            yne.append(i + 1)
    return ene, yne
