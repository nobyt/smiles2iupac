"""倍数命名法 (Multiplicative Nomenclature, IUPAC 2013 P-15.3, P-54).

対称な二官能性化合物が二価の連結基（-O-, -S-, -SO2-, -NH-, -CH2-, -CH2CH2- 等）で
結合された化合物を IUPAC Preferred IUPAC Name (PIN) に従って命名する。
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .molecule_analyzer import MoleculeGraph


def _name_multiplicative(graph: "MoleculeGraph", get_atom) -> str | None:
    """倍数命名法が適用可能な対称構造を検出して IUPAC 名を返す。"""
    from .molecule_analyzer import get_bond_order

    # 1. 2原子連結基 (-CH2-CH2-) の検出
    res_2c = _try_two_carbon_linker(graph, get_atom)
    if res_2c is not None:
        return res_2c

    # 1b. 2硫黄連結基 (-S-S-, disulfanediyl) の検出 (Phase 945)
    res_ss = _try_disulfanediyl_linker(graph, get_atom)
    if res_ss is not None:
        return res_ss

    # 2. 単一原子連結基 (-O-, -S-, -SO2-, -SO-, -NH-, -CH2-, -CO-) の検出
    for l_idx, l_atom in enumerate(graph.atoms):
        sym = l_atom.symbol
        if sym not in ("O", "S", "N", "C"):
            continue
        if l_atom.in_ring:
            continue

        # 隣接重原子
        heavy_nbs = [nb for nb in graph.adjacency[l_idx] if get_atom(graph, nb).symbol != "H"]

        # --- (A) -O- (oxy) ---
        if sym == "O" and len(heavy_nbs) == 2:
            c1, c2 = heavy_nbs
            if get_atom(graph, c1).symbol == "C" and get_atom(graph, c2).symbol == "C":
                res = _evaluate_symmetric_branches(graph, l_idx, c1, c2, "oxy", get_atom)
                if res:
                    return res

        # --- (B) -S-, -S(=O)-, -S(=O)2- (thio, sulfinyl, sulfonyl) ---
        elif sym == "S":
            o_dbl_nbs = [
                nb for nb in heavy_nbs
                if get_atom(graph, nb).symbol == "O" and get_bond_order(graph, l_idx, nb) == 2.0
            ]
            c_nbs = [nb for nb in heavy_nbs if get_atom(graph, nb).symbol == "C"]
            if len(c_nbs) == 2:
                c1, c2 = c_nbs
                if len(o_dbl_nbs) == 0 and len(heavy_nbs) == 2:
                    res = _evaluate_symmetric_branches(graph, l_idx, c1, c2, "thio", get_atom)
                    if res:
                        return res
                elif len(o_dbl_nbs) == 1 and len(heavy_nbs) == 3:
                    res = _evaluate_symmetric_branches(graph, l_idx, c1, c2, "sulfinyl", get_atom)
                    if res:
                        return res
                elif len(o_dbl_nbs) == 2 and len(heavy_nbs) == 4:
                    res = _evaluate_symmetric_branches(graph, l_idx, c1, c2, "sulfonyl", get_atom)
                    if res:
                        return res

        # --- (C) -NH- (imino / azanediyl) ---
        elif sym == "N" and len(heavy_nbs) == 2 and l_atom.formal_charge == 0:
            c1, c2 = heavy_nbs
            if get_atom(graph, c1).symbol == "C" and get_atom(graph, c2).symbol == "C":
                # カルボン酸等には imino, その他には azanediyl
                res = _evaluate_symmetric_branches(graph, l_idx, c1, c2, "imino", get_atom)
                if res:
                    return res

        # --- (D) -CH2- (methylene) or -C(=O)- (carbonyl) ---
        elif sym == "C":
            o_dbl_nbs = [
                nb for nb in heavy_nbs
                if get_atom(graph, nb).symbol == "O" and get_bond_order(graph, l_idx, nb) == 2.0
            ]
            c_nbs = [nb for nb in heavy_nbs if get_atom(graph, nb).symbol == "C"]
            if len(c_nbs) == 2:
                c1, c2 = c_nbs
                if len(o_dbl_nbs) == 0 and len(heavy_nbs) == 2:
                    res = _evaluate_symmetric_branches(graph, l_idx, c1, c2, "methylene", get_atom)
                    if res:
                        return res
                elif len(o_dbl_nbs) == 1 and len(heavy_nbs) == 3:
                    res = _evaluate_symmetric_branches(graph, l_idx, c1, c2, "carbonyl", get_atom)
                    if res:
                        return res

    return None


def _try_two_carbon_linker(graph: "MoleculeGraph", get_atom) -> str | None:
    """-CH2-CH2- 連結基 ((ethane-1,2-diyl)) の倍数命名を試みる。"""
    from .molecule_analyzer import get_bond_order

    for c1_idx, c1_atom in enumerate(graph.atoms):
        if c1_atom.symbol != "C" or c1_atom.in_ring:
            continue
        c1_heavy = [nb for nb in graph.adjacency[c1_idx] if get_atom(graph, nb).symbol != "H"]
        if len(c1_heavy) != 2:
            continue

        for c2_idx in c1_heavy:
            if c2_idx <= c1_idx:
                continue
            c2_atom = get_atom(graph, c2_idx)
            if c2_atom.symbol != "C" or c2_atom.in_ring:
                continue
            if get_bond_order(graph, c1_idx, c2_idx) != 1.0:
                continue

            c2_heavy = [nb for nb in graph.adjacency[c2_idx] if get_atom(graph, nb).symbol != "H"]
            if len(c2_heavy) != 2:
                continue

            ext1 = next(nb for nb in c1_heavy if nb != c2_idx)
            ext2 = next(nb for nb in c2_heavy if nb != c1_idx)

            if get_atom(graph, ext1).symbol != "C" or get_atom(graph, ext2).symbol != "C":
                continue

            # 連結基原子セット = {c1_idx, c2_idx}
            linker_set = {c1_idx, c2_idx}
            res = _evaluate_branches_from_set(graph, linker_set, ext1, ext2, "(ethane-1,2-diyl)", get_atom)
            if res:
                return res

    return None


def _try_disulfanediyl_linker(graph: "MoleculeGraph", get_atom) -> str | None:
    """-S-S- 連結基 (disulfanediyl) の倍数命名を試みる。"""
    from .molecule_analyzer import get_bond_order

    def _pure_disulfide_s(s_idx: int) -> tuple[int, int] | None:
        """S が純粋なジスルフィド硫黄 (二重結合Oなし・重原子隣接2つ) なら
        (相手S, C枝) を返す。"""
        s_atom = get_atom(graph, s_idx)
        if s_atom.symbol != "S" or s_atom.in_ring or s_atom.formal_charge != 0:
            return None
        heavy = [nb for nb in graph.adjacency[s_idx] if get_atom(graph, nb).symbol != "H"]
        if len(heavy) != 2:
            return None
        # 二重結合Oを持つ (sulfinyl/sulfonyl系) は対象外
        if any(get_atom(graph, nb).symbol == "O" and get_bond_order(graph, s_idx, nb) == 2.0
               for nb in heavy):
            return None
        return tuple(heavy)  # 2隣接

    for s1_idx, s1_atom in enumerate(graph.atoms):
        info1 = _pure_disulfide_s(s1_idx)
        if info1 is None:
            continue
        # s1 の隣に S があるか
        s2_idx = next((nb for nb in info1 if get_atom(graph, nb).symbol == "S"), None)
        if s2_idx is None or s2_idx <= s1_idx:
            continue
        if get_bond_order(graph, s1_idx, s2_idx) != 1.0:
            continue
        info2 = _pure_disulfide_s(s2_idx)
        if info2 is None:
            continue

        ext1 = next((nb for nb in info1 if nb != s2_idx), None)
        ext2 = next((nb for nb in info2 if nb != s1_idx), None)
        if ext1 is None or ext2 is None:
            continue
        if get_atom(graph, ext1).symbol != "C" or get_atom(graph, ext2).symbol != "C":
            continue

        res = _evaluate_branches_from_set(
            graph, {s1_idx, s2_idx}, ext1, ext2, "disulfanediyl", get_atom
        )
        if res:
            return res

    return None


def _evaluate_symmetric_branches(
    graph: "MoleculeGraph",
    linker_idx: int,
    c1: int,
    c2: int,
    linker_name: str,
    get_atom,
) -> str | None:
    return _evaluate_branches_from_set(graph, {linker_idx}, c1, c2, linker_name, get_atom)


def _evaluate_branches_from_set(
    graph: "MoleculeGraph",
    linker_set: set[int],
    c1: int,
    c2: int,
    linker_name: str,
    get_atom,
) -> str | None:
    """linker_set から伸びる c1, c2 の2つの枝が同一対称かを判定し倍数名を返す。"""
    from .molecule_analyzer import get_bond_order

    def _collect_branch(root: int) -> set[int]:
        visited: set[int] = set()
        stack = [root]
        while stack:
            curr = stack.pop()
            if curr in visited or curr in linker_set:
                continue
            visited.add(curr)
            for nb in graph.adjacency[curr]:
                if nb not in visited and nb not in linker_set:
                    stack.append(nb)
        return visited

    b1 = _collect_branch(c1)
    b2 = _collect_branch(c2)

    # 枝同士が重複している場合は非対称 / 環状
    if b1 & b2:
        return None
    if len(b1) != len(b2):
        return None

    # 1. カルボン酸鎖の判定 (-CH2-COOH or -CH2-CH2-COOH)
    res_acid = _check_acid_branch(graph, b1, b2, c1, c2, linker_name, get_atom)
    if res_acid is not None:
        return res_acid

    # 2. 芳香族官能基（4-アミノフェニル, 4-ヒドロキシフェニル, 2-ピリジル）の判定
    res_arom = _check_aromatic_branch(graph, b1, b2, c1, c2, linker_name, get_atom)
    if res_arom is not None:
        return res_arom

    return None


def _check_acid_branch(
    graph: "MoleculeGraph",
    b1: set[int],
    b2: set[int],
    c1: int,
    c2: int,
    linker_name: str,
    get_atom,
) -> str | None:
    """カルボン酸枝 (-CH2-COOH, -CH2CH2-COOH) を判定。
    
    倍数命名法 (P-54) は、ヘテロ原子二価連結基 (-O-, -S-, -SO2-, -NH-) の場合にのみ適用する。
    炭素連結基 (-CH2-, -CH2CH2-) の直鎖ジカルボン酸 (pentanedioic acid, adipic acid) を誤命名しない。
    """
    if linker_name not in ("oxy", "thio", "sulfinyl", "sulfonyl", "imino", "disulfanediyl"):
        return None

    from .molecule_analyzer import get_bond_order

    def _acid_info(branch_atoms: set[int], root_c: int) -> tuple[int, str] | None:
        # root_c からの直鎖カルボン酸を探索
        # 枝内の重原子を収集
        heavy = {a for a in branch_atoms if get_atom(graph, a).symbol != "H"}
        # 炭素数カウント
        cs = [a for a in heavy if get_atom(graph, a).symbol == "C"]
        os = [a for a in heavy if get_atom(graph, a).symbol == "O"]
        if len(os) != 2:
            return None
        # アニオン (formal charge -1) は通常のアニオン命名ハンドラに任せる
        if any(get_atom(graph, o).formal_charge != 0 for o in os):
            return None
        # COOH を持つ炭素を特定
        cooh_c = next(
            (c for c in cs if sum(1 for nb in graph.adjacency[c] if nb in os) == 2),
            None,
        )
        if cooh_c is None:
            return None

        # root_c から cooh_c までの鎖長
        if len(cs) == 2 and root_c != cooh_c:
            # -CH2-COOH -> acetic acid, locant 2
            return 2, "acetic acid"
        elif len(cs) == 3 and root_c != cooh_c:
            # -CH2-CH2-COOH -> propanoic acid, locant 3
            return 3, "propanoic acid"
        return None

    info1 = _acid_info(b1, c1)
    info2 = _acid_info(b2, c2)
    if info1 and info2 and info1 == info2:
        locant, acid_stem = info1
        # linker_name が imino の場合は 2,2'-iminodiacetic acid
        l_name = "imino" if linker_name == "imino" else linker_name
        return f"{locant},{locant}'-{l_name}di{acid_stem}"

    return None


def _check_aromatic_branch(
    graph: "MoleculeGraph",
    b1: set[int],
    b2: set[int],
    c1: int,
    c2: int,
    linker_name: str,
    get_atom,
) -> str | None:
    """芳香族官能基枝（aniline, phenol, pyridine）を判定。"""
    from .molecule_analyzer import get_bond_order

    def _aromatic_info(branch_atoms: set[int], root_c: int) -> tuple[int, str] | None:
        heavy = {a for a in branch_atoms if get_atom(graph, a).symbol != "H"}
        # 6員環を探す
        ring_six = None
        for rt in (graph.ring_atom_sets or []):
            if set(rt).issubset(heavy) and len(rt) == 6:
                ring_six = list(rt)
                break
        if ring_six is None:
            return None

        ring_set = set(ring_six)
        if root_c not in ring_set:
            return None

        # 環内原子の元素
        ring_syms = [get_atom(graph, a).symbol for a in ring_six]

        # Case 1: 全て炭素 (benzene ring)
        if all(s == "C" for s in ring_syms):
            # 環外置換基を探す (root_c 以外)
            ext_subs = []
            for ra in ring_six:
                if ra == root_c:
                    continue
                for nb in graph.adjacency[ra]:
                    if nb in heavy and nb not in ring_set:
                        ext_subs.append((ra, nb))

            if len(ext_subs) == 1:
                sub_pos_atom, sub_atom_idx = ext_subs[0]
                # root_c からの距離を測る
                # 6員環で反対側 = 最短距離 3 (4位)
                # 環の巡回順
                def _order_6(start, rset):
                    ord_ = [start]
                    prev = -1
                    curr = start
                    while len(ord_) < 6:
                        nexts = [nb for nb in graph.adjacency[curr] if nb in rset and nb != prev and nb not in ord_]
                        if not nexts:
                            break
                        prev, curr = curr, nexts[0]
                        ord_.append(curr)
                    return ord_
                ring_ord = _order_6(root_c, ring_set)
                pos_in_ring = ring_ord.index(sub_pos_atom) + 1
                locant = 4 if pos_in_ring in (4,) else (pos_in_ring if pos_in_ring <= 3 else 8 - pos_in_ring)

                sub_sym = get_atom(graph, sub_atom_idx).symbol
                # -NH2 (aniline)
                if sub_sym == "N" and heavy == (ring_set | {sub_atom_idx}):
                    return locant, "aniline"
                # -OH (phenol)
                elif sub_sym == "O" and heavy == (ring_set | {sub_atom_idx}):
                    return locant, "phenol"

        # Case 2: ピリジン環 (1個の N)
        elif ring_syms.count("N") == 1 and ring_syms.count("C") == 5:
            if heavy == ring_set:
                # root_c と N の位置関係
                n_idx = next(a for a in ring_six if get_atom(graph, a).symbol == "N")
                # N が 1位のとき root_c のロカント
                # N に隣接していれば 2位
                if root_c in graph.adjacency[n_idx]:
                    return 2, "pyridine"

        return None

    info1 = _aromatic_info(b1, c1)
    info2 = _aromatic_info(b2, c2)
    if info1 and info2 and info1 == info2:
        locant, parent_name = info1
        l_name = "azanediyl" if linker_name == "imino" else linker_name
        return f"{locant},{locant}'-{l_name}di{parent_name}"

    return None
