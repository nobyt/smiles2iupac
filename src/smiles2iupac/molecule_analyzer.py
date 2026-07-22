"""
RDKit Mol オブジェクトから内部グラフ表現 MoleculeGraph を構築する。
"""

from __future__ import annotations

from dataclasses import dataclass, field

from rdkit import Chem
from rdkit.Chem import rdchem


@dataclass
class AtomInfo:
    idx: int
    symbol: str           # 'C', 'O', 'N', 'H', ...
    atomic_num: int
    is_aromatic: bool
    in_ring: bool
    num_hs: int           # implicit + explicit H (AddHs後は explicit のみ)
    chiral_tag: str | None  # 'R', 'S', or None
    formal_charge: int = 0  # 形式電荷 (Phase 146)


@dataclass
class BondInfo:
    begin_idx: int
    end_idx: int
    bond_order: float     # 1.0, 2.0, 3.0, 1.5(芳香族)
    stereo: str | None    # 'E', 'Z', or None


@dataclass
class MoleculeGraph:
    atoms: list[AtomInfo]
    bonds: list[BondInfo]
    adjacency: dict[int, list[int]]          # atom_idx -> [neighbor atom_idx]
    bond_orders: dict[tuple[int, int], float]  # (i, j) -> bond_order (i<j)
    ring_atom_sets: list[tuple[int, ...]] = field(default_factory=list)
    # RDKit SSSR (各リングを構成する重原子インデックスのタプル)
    rdkit_mol: object = None
    # AddHs() 前の RDKit Mol オブジェクト（canonical SMILES 生成用）


def build_molecule_graph(smiles: str) -> MoleculeGraph:
    """
    SMILESをパースし RDKit Mol → MoleculeGraph に変換する。

    Args:
        smiles: 有効なSMILES文字列

    Returns:
        MoleculeGraph

    Raises:
        ValueError: 無効なSMILES
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles!r}")

    # 立体化学情報を割り当て
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    # Phase 873: 疑似不斉中心 (1,4-二置換シクロヘキサン等) の小文字 r/s を含む
    # 正確な CIP ラベルを付与。legacy AssignStereochemistry は疑似不斉中心に
    # _CIPCode を付けず環立体化学が欠落していた。既存の大文字 R/S は不変
    # (rdCIPLabeler は legacy と一致、疑似不斉のみ追加)。
    try:
        from rdkit.Chem import rdCIPLabeler
        rdCIPLabeler.AssignCIPLabels(mol)
    except Exception:
        pass  # rdCIPLabeler 不在時は legacy ラベルのままフォールバック

    # SSSR (最小環系) を AddHs 前に取得（重原子インデックスが変わらない）
    ring_info = mol.GetRingInfo()
    ring_atom_sets: list[tuple[int, ...]] = list(ring_info.AtomRings())

    # AddHs() 前の mol を保存（canonical SMILES 生成用）
    mol_no_hs = mol

    # 水素を明示的に追加（官能基検出の精度向上）
    mol = Chem.AddHs(mol)

    atoms: list[AtomInfo] = []
    for atom in mol.GetAtoms():
        props = atom.GetPropsAsDict()
        chiral = props.get("_CIPCode")  # 'R', 'S', or None

        atoms.append(
            AtomInfo(
                idx=atom.GetIdx(),
                symbol=atom.GetSymbol(),
                atomic_num=atom.GetAtomicNum(),
                is_aromatic=atom.GetIsAromatic(),
                in_ring=atom.IsInRing(),
                num_hs=atom.GetTotalNumHs(),
                chiral_tag=chiral,
                formal_charge=atom.GetFormalCharge(),
            )
        )

    bonds: list[BondInfo] = []
    bond_orders: dict[tuple[int, int], float] = {}
    adjacency: dict[int, list[int]] = {a.idx: [] for a in atoms}

    for bond in mol.GetBonds():
        bi = bond.GetBeginAtomIdx()
        ei = bond.GetEndAtomIdx()
        bo = bond.GetBondTypeAsDouble()

        # E/Z stereo
        # Phase 873: rdCIPLabeler は bond stereo enum を STEREOTRANS/CIS に変換し
        # 正しい E/Z を _CIPCode に格納する。まず _CIPCode を読み、無ければ
        # legacy の STEREOE/STEREOZ enum にフォールバック。
        stereo = None
        _bond_cip = bond.GetPropsAsDict().get("_CIPCode")
        if _bond_cip in ("E", "Z"):
            stereo = _bond_cip
        else:
            s = bond.GetStereo()
            if s == rdchem.BondStereo.STEREOE:
                stereo = "E"
            elif s == rdchem.BondStereo.STEREOZ:
                stereo = "Z"

        bonds.append(BondInfo(begin_idx=bi, end_idx=ei, bond_order=bo, stereo=stereo))

        key = (min(bi, ei), max(bi, ei))
        bond_orders[key] = bo
        adjacency[bi].append(ei)
        adjacency[ei].append(bi)

    return MoleculeGraph(
        atoms=atoms,
        bonds=bonds,
        adjacency=adjacency,
        bond_orders=bond_orders,
        ring_atom_sets=ring_atom_sets,
        rdkit_mol=mol_no_hs,
    )


def get_bond_order(graph: MoleculeGraph, i: int, j: int) -> float:
    """2原子間の結合次数を返す。結合がなければ 0.0。"""
    key = (min(i, j), max(i, j))
    return graph.bond_orders.get(key, 0.0)


def get_atom(graph: MoleculeGraph, idx: int) -> AtomInfo:
    """インデックスで AtomInfo を返す。"""
    return graph.atoms[idx]


def carbon_indices(graph: MoleculeGraph) -> list[int]:
    """水素以外の炭素原子インデックス一覧。"""
    return [a.idx for a in graph.atoms if a.symbol == "C"]


def non_ring_carbon_indices(graph: MoleculeGraph) -> list[int]:
    """環外炭素のインデックス一覧。"""
    return [a.idx for a in graph.atoms if a.symbol == "C" and not a.in_ring]
