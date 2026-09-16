"""ホスホニウムイリド・スルホニウムイリド命名法 (Ylides / Betaines, IUPAC 2013 P-74, P-67, P-68).

ウィッティヒ試薬等のホスホニウムイリドおよびスルホニウムイリドを
IUPAC 2013 lambda 命名法 (lambda5-phosphane, lambda4-sulfane) により命名する。
"""

from __future__ import annotations

from typing import TYPE_CHECKING
from collections import defaultdict
import re

if TYPE_CHECKING:
    from .molecule_analyzer import MoleculeGraph

from .constants import MULTIPLIER, CHAIN_PREFIX
from .substituent import name_substituent


def _name_ylide(graph: "MoleculeGraph", get_atom) -> str | None:
    """ホスホニウムイリドおよびスルホニウムイリドを検出して IUPAC 名を返す。"""
    from .molecule_analyzer import get_bond_order

    for h_idx, h_atom in enumerate(graph.atoms):
        sym = h_atom.symbol
        if sym not in ("P", "S"):
            continue
        if h_atom.in_ring:
            continue

        # イリド炭素 (ylidene C) を探索
        # パターン 1: 二重結合 C
        # パターン 2: 双極性 [P+]-[C-] または [S+]-[C-]
        ylide_c: int | None = None
        other_subs: list[int] = []

        for nb_idx in graph.adjacency[h_idx]:
            nb = get_atom(graph, nb_idx)
            if nb.symbol == "H":
                continue
            bo = get_bond_order(graph, h_idx, nb_idx)
            if nb.symbol == "C":
                if bo == 2.0:
                    ylide_c = nb_idx
                elif h_atom.formal_charge == 1 and nb.formal_charge == -1:
                    ylide_c = nb_idx
                else:
                    other_subs.append(nb_idx)
            else:
                other_subs.append(nb_idx)

        if ylide_c is None:
            continue

        # ヘテロ原子の母体名
        if sym == "P":
            parent = "lambda5-phosphane"
            expected_other = 3
        else:  # S
            parent = "lambda4-sulfane"
            expected_other = 2

        if len(other_subs) != expected_other:
            continue

        # 1. イリド炭素置換基 (ylidene) の命名
        ylidene_name = _name_ylidene(graph, ylide_c, excluded={h_idx}, get_atom=get_atom)
        if ylidene_name is None:
            continue

        # 2. その他の単結合置換基の命名
        sub_names: list[str] = []
        excluded_for_subs = {h_idx, ylide_c}
        for sub_c in other_subs:
            sn = name_substituent(graph, sub_c, excluded_for_subs)
            sub_names.append(sn)

        # 3. 置換基の集約とソート
        grouped: dict[str, int] = defaultdict(int)
        for sn in sub_names:
            grouped[sn] += 1

        parts: list[tuple[str, str]] = []  # (sort_key, formatted_string)

        # ylidene 置換基
        yl_fmt = f"({ylidene_name})"
        yl_sort = re.sub(r"[^a-zA-Z]", "", ylidene_name)
        parts.append((yl_sort, yl_fmt))

        # その他の置換基
        for sn in sorted(grouped.keys()):
            cnt = grouped[sn]
            mult = MULTIPLIER.get(cnt, str(cnt)) if cnt > 1 else ""
            fmt = f"{mult}{sn}"
            sort_k = re.sub(r"[^a-zA-Z]", "", sn)
            parts.append((sort_k, fmt))

        # アルファベット順にソート
        parts.sort(key=lambda t: t[0])

        prefix_str = "".join(t[1] for t in parts)
        return f"{prefix_str}-{parent}"

    return None


def _name_ylidene(graph: "MoleculeGraph", c_idx: int, excluded: set[int], get_atom) -> str | None:
    """イリド炭素からなる ylidene 置換基（methylidene, ethylidene, propan-2-ylidene 等）を命名。"""
    heavy_nbs = [nb for nb in graph.adjacency[c_idx] if nb not in excluded and get_atom(graph, nb).symbol != "H"]

    # 炭素数 1: =CH2 -> methylidene
    if not heavy_nbs:
        return "methylidene"

    # 炭素数 2: =CH-CH3 -> ethylidene
    if len(heavy_nbs) == 1:
        c1 = heavy_nbs[0]
        c1_heavy = [nb for nb in graph.adjacency[c1] if nb != c_idx and nb not in excluded and get_atom(graph, nb).symbol != "H"]
        if not c1_heavy and get_atom(graph, c1).symbol == "C":
            return "ethylidene"
        if len(c1_heavy) == 1 and get_atom(graph, c1).symbol == "C":
            c2 = c1_heavy[0]
            c2_heavy = [nb for nb in graph.adjacency[c2] if nb != c1 and nb not in excluded and get_atom(graph, nb).symbol != "H"]
            if not c2_heavy and get_atom(graph, c2).symbol == "C":
                return "propylidene"

    # 炭素数 3 (分岐): =C(CH3)2 -> propan-2-ylidene
    if len(heavy_nbs) == 2:
        c1, c2 = heavy_nbs
        c1_heavy = [nb for nb in graph.adjacency[c1] if nb != c_idx and nb not in excluded and get_atom(graph, nb).symbol != "H"]
        c2_heavy = [nb for nb in graph.adjacency[c2] if nb != c_idx and nb not in excluded and get_atom(graph, nb).symbol != "H"]
        if not c1_heavy and not c2_heavy and get_atom(graph, c1).symbol == "C" and get_atom(graph, c2).symbol == "C":
            return "propan-2-ylidene"

    return "methylidene"
