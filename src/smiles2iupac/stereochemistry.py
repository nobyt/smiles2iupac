"""
立体化学記述子の生成。

RDKit が割り当てた CIP コード (R/S, E/Z) を読み取り、
IUPAC 命名法の立体記述子文字列リストを返す。
"""

from __future__ import annotations

from .chain_finder import PrincipalChain
from .molecule_analyzer import MoleculeGraph, get_atom, get_bond_order


def assign_stereochemistry(
    graph: MoleculeGraph,
    chain: PrincipalChain,
) -> list[str]:
    """
    主鎖上の立体中心に対して R/S・E/Z 記述子を返す。

    Returns:
        ['(2R)', '(3E)', ...] など。立体化学がなければ空リスト。
    """
    descriptors: list[str] = []

    path = chain.atom_indices
    locant_map = chain.locant_map

    # R/S: 主鎖炭素の不斉中心
    for c_idx in path:
        atom = get_atom(graph, c_idx)
        cip = atom.chiral_tag  # 'R', 'S', or None
        if cip is not None:
            locant = locant_map[c_idx]
            descriptors.append(f"({locant}{cip})")

    # E/Z: 主鎖上の C=C
    for i in range(len(path) - 1):
        bo = get_bond_order(graph, path[i], path[i + 1])
        if bo == 2.0:
            # BondInfo から stereo を取得
            stereo = _get_bond_stereo(graph, path[i], path[i + 1])
            if stereo is not None:
                locant = locant_map[path[i]]  # 低い方のロカント
                descriptors.append(f"({locant}{stereo})")

    return descriptors


def _get_bond_stereo(
    graph: MoleculeGraph, idx1: int, idx2: int
) -> str | None:
    """2原子間の結合の E/Z ステレオを返す。"""
    for bond in graph.bonds:
        if (bond.begin_idx == idx1 and bond.end_idx == idx2) or \
           (bond.begin_idx == idx2 and bond.end_idx == idx1):
            return bond.stereo
    return None
