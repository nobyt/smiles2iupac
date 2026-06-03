"""
smiles2iupac: SMILES から IUPAC 名を生成するライブラリ。

IUPACの命名ロジックはスクラッチ実装。
SMILESのパース・分子グラフ構築は RDKit に委譲。

対応範囲（Phase 1-3）:
  Phase 1-2:
  - 直鎖・分岐アルカン、アルコール、ハロゲン化物
  - カルボン酸、アルデヒド、ケトン
  - アルケン・アルキン
  - 立体化学 (R/S, E/Z)
  Phase 3 (complete):
  - シクロアルカン (cyclopropane ~ cyclodecane)
  - 置換シクロアルカン (methylcyclohexane 等)
  - シクロアルカンの官能基 (cyclohexanol, cyclohexanone 等)
  - ベンゼン
  - 置換ベンゼン (chlorobenzene, 1,2-dimethylbenzene 等)
"""

from __future__ import annotations

__all__ = ["smiles_to_iupac"]

from .group_namers import (
    PGRP_DISPATCH,
    _name_thiourea_if_match,
    _name_substituted_urea_if_match,
    _is_urea,
    _is_carbonohydrazide,
    _name_hetero_n_oxide,
    _name_amine_n_oxide,
    _name_n_substituted_hydroxylamine,
    _name_nitrone,
    _name_diazo_compound,
    _name_nitrosamine,
    _name_hydrazine_compound,
    _name_azo_compound,
    _alkylidene_name,
    _name_n_substituted_imine,
    _name_secondary_tertiary_amide,
    _name_substituted_amidine,
    _name_secondary_tertiary_amine,
    _name_cyanamide,
    _name_acyl_peroxide,
    _name_sulfate_ester,
    _name_sulfite_ester,
    _name_diazonium,
)


# Phase 135/136: 保留名テーブル (IUPAC 2013 P-12.1, P-65.1.1.4, P-65.1.3)
# キー = RDKit canonical SMILES (stereo あり/なし両方)
_RETAINED_NAMES: dict[str, str] = {
    "NCC(=O)O":                          "glycine",
    # alanine
    "C[C@H](N)C(=O)O":                  "L-alanine",
    "C[C@@H](N)C(=O)O":                 "D-alanine",
    "CC(N)C(=O)O":                       "alanine",
    # valine
    "CC(C)[C@H](N)C(=O)O":              "L-valine",
    "CC(C)[C@@H](N)C(=O)O":             "D-valine",
    "CC(C)C(N)C(=O)O":                   "valine",
    # leucine
    "CC(C)C[C@H](N)C(=O)O":             "L-leucine",
    "CC(C)C[C@@H](N)C(=O)O":            "D-leucine",
    "CC(C)CC(N)C(=O)O":                  "leucine",
    # isoleucine
    "CC[C@H](C)[C@H](N)C(=O)O":         "L-isoleucine",
    # proline
    "O=C(O)[C@@H]1CCCN1":               "L-proline",
    "O=C(O)[C@H]1CCCN1":                "D-proline",
    "OC(=O)C1CCCN1":                     "proline",
    # phenylalanine
    "N[C@@H](Cc1ccccc1)C(=O)O":          "L-phenylalanine",
    "N[C@H](Cc1ccccc1)C(=O)O":          "D-phenylalanine",
    "NC(Cc1ccccc1)C(=O)O":               "phenylalanine",
    # tryptophan
    "N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O": "L-tryptophan",
    "N[C@H](Cc1c[nH]c2ccccc12)C(=O)O":  "D-tryptophan",
    # methionine
    "CSCC[C@H](N)C(=O)O":               "L-methionine",
    "CSCC[C@@H](N)C(=O)O":              "D-methionine",
    "CSCCC(N)C(=O)O":                    "methionine",
    # serine
    "N[C@@H](CO)C(=O)O":                 "L-serine",
    "N[C@H](CO)C(=O)O":                 "D-serine",
    "NC(CO)C(=O)O":                      "serine",
    # threonine
    "C[C@H](O)[C@H](N)C(=O)O":          "L-threonine",
    "C[C@@H](O)[C@@H](N)C(=O)O":        "D-threonine",
    # cysteine
    "N[C@@H](CS)C(=O)O":                 "L-cysteine",
    "N[C@H](CS)C(=O)O":                 "D-cysteine",
    "NC(CS)C(=O)O":                      "cysteine",
    # tyrosine
    "N[C@@H](Cc1ccc(O)cc1)C(=O)O":       "L-tyrosine",
    "N[C@H](Cc1ccc(O)cc1)C(=O)O":       "D-tyrosine",
    # asparagine
    "NC(=O)C[C@H](N)C(=O)O":            "L-asparagine",
    "NC(=O)C[C@@H](N)C(=O)O":           "D-asparagine",
    "NC(=O)CC(N)C(=O)O":                 "asparagine",
    # glutamine
    "NC(=O)CC[C@H](N)C(=O)O":           "L-glutamine",
    "NC(=O)CC[C@@H](N)C(=O)O":          "D-glutamine",
    "NC(=O)CCC(N)C(=O)O":                "glutamine",
    # lysine
    "NCCCC[C@H](N)C(=O)O":              "L-lysine",
    "NCCCC[C@@H](N)C(=O)O":             "D-lysine",
    "NCCCCC(N)C(=O)O":                   "lysine",
    # arginine
    "N=C(N)NCCC[C@H](N)C(=O)O":         "L-arginine",
    "N=C(N)NCCC[C@@H](N)C(=O)O":        "D-arginine",
    # histidine
    "N[C@@H](Cc1cnc[nH]1)C(=O)O":        "L-histidine",
    "N[C@H](Cc1cnc[nH]1)C(=O)O":        "D-histidine",
    "NC(Cc1cnc[nH]1)C(=O)O":             "histidine",
    # aspartic acid
    "N[C@@H](CC(=O)O)C(=O)O":            "L-aspartic acid",
    "N[C@H](CC(=O)O)C(=O)O":            "D-aspartic acid",
    "NC(CC(=O)O)C(=O)O":                 "aspartic acid",
    # glutamic acid
    "N[C@@H](CCC(=O)O)C(=O)O":           "L-glutamic acid",
    "N[C@H](CCC(=O)O)C(=O)O":           "D-glutamic acid",
    "NC(CCC(=O)O)C(=O)O":                "glutamic acid",

    # ── Phase 136: 二塩基酸・ヒドロキシ酸 保留名 (P-65.1.1.4, P-65.1.3) ──
    # aliphatic diacids
    "O=C(O)C(=O)O":                      "oxalic acid",
    "O=C(O)CC(=O)O":                     "malonic acid",
    "O=C(O)CCC(=O)O":                    "succinic acid",
    "O=C(O)CCCC(=O)O":                   "glutaric acid",
    "O=C(O)CCCCC(=O)O":                  "adipic acid",
    "O=C(O)CCCCCC(=O)O":                 "pimelic acid",
    "O=C(O)CCCCCCC(=O)O":               "suberic acid",
    "O=C(O)CCCCCCCC(=O)O":              "azelaic acid",
    "O=C(O)CCCCCCCCC(=O)O":             "sebacic acid",
    # unsaturated diacids
    "O=C(O)/C=C\\C(=O)O":               "maleic acid",
    "O=C(O)/C=C/C(=O)O":                "fumaric acid",
    # cinnamic acid
    "O=C(O)/C=C/c1ccccc1":              "(E)-cinnamic acid",
    "O=C(O)/C=C\\c1ccccc1":             "(Z)-cinnamic acid",
    # aromatic diacids (phthalic, isophthalic, terephthalic) P-65.1.1.4
    "O=C(O)c1ccccc1C(=O)O":             "phthalic acid",
    "O=C(O)c1cccc(C(=O)O)c1":           "isophthalic acid",
    "O=C(O)c1ccc(C(=O)O)cc1":           "terephthalic acid",
    # ── Phase 270: benzene tri/tetracarboxylic acids (IUPAC 2013 P-65.1.2.3) ──
    # tricarboxylic (3 isomers)
    "O=C(O)c1cccc(C(=O)O)c1C(=O)O":    "benzene-1,2,3-tricarboxylic acid",
    "O=C(O)c1ccc(C(=O)O)c(C(=O)O)c1":  "benzene-1,2,4-tricarboxylic acid",
    "O=C(O)c1cc(C(=O)O)cc(C(=O)O)c1":  "benzene-1,3,5-tricarboxylic acid",
    # tetracarboxylic (3 isomers)
    "O=C(O)c1ccc(C(=O)O)c(C(=O)O)c1C(=O)O": "benzene-1,2,3,4-tetracarboxylic acid",
    "O=C(O)c1cc(C(=O)O)c(C(=O)O)c(C(=O)O)c1": "benzene-1,2,3,5-tetracarboxylic acid",
    "O=C(O)c1cc(C(=O)O)c(C(=O)O)cc1C(=O)O": "benzene-1,2,4,5-tetracarboxylic acid",
    # pentacarboxylic
    "O=C(O)c1cc(C(=O)O)c(C(=O)O)c(C(=O)O)c1C(=O)O": "benzene-1,2,3,4,5-pentacarboxylic acid",
    # hexacarboxylic (mellitic acid — retained name, IUPAC 2013 P-65.1.2.3)
    "O=C(O)c1c(C(=O)O)c(C(=O)O)c(C(=O)O)c(C(=O)O)c1C(=O)O": "mellitic acid",
    # lactic acid P-65.1.1.4
    "C[C@@H](O)C(=O)O":                 "L-lactic acid",
    "C[C@H](O)C(=O)O":                  "D-lactic acid",
    "CC(O)C(=O)O":                       "lactic acid",
    # malic acid
    "O=C(O)C[C@H](O)C(=O)O":           "L-malic acid",
    "O=C(O)C[C@@H](O)C(=O)O":          "D-malic acid",
    "O=C(O)CC(O)C(=O)O":               "malic acid",
    # tartaric acid
    "O=C(O)[C@@H](O)[C@H](O)C(=O)O":   "L-tartaric acid",
    "O=C(O)[C@H](O)[C@@H](O)C(=O)O":   "D-tartaric acid",
    "O=C(O)[C@@H](O)[C@@H](O)C(=O)O":  "meso-tartaric acid",
    "O=C(O)C(O)C(O)C(=O)O":            "tartaric acid",
    # citric acid
    "O=C(O)CC(O)(CC(=O)O)C(=O)O":      "citric acid",

    # ── Phase 137: 核酸塩基 保留名 (IUPAC 2013 P-14.5) ──
    "Nc1ncnc2[nH]cnc12":               "adenine",
    "Nc1nc2[nH]cnc2c(=O)[nH]1":        "guanine",
    "Nc1ccnc(=O)[nH]1":                "cytosine",
    "O=c1cc[nH]c(=O)[nH]1":           "uracil",
    "Cc1c[nH]c(=O)[nH]c1=O":          "thymine",
    "O=c1[nH]cnc2[nH]cnc12":          "hypoxanthine",
    "O=c1[nH]c(=O)c2[nH]cnc2[nH]1":  "xanthine",
    "c1ncc2[nH]cnc2n1":               "9H-purine",

    # ── Phase 146: 無機オキソ酸・オキソ酸アニオン 保留名 ──
    "O=C([O-])O":                      "bicarbonate",
    "O=C([O-])[O-]":                   "carbonate",
    "O=S(=O)(O)O":                     "sulfuric acid",
    "O=S(O)O":                         "sulfurous acid",
    "O=[N+]([O-])O":                   "nitric acid",
    "O=NO":                            "nitrous acid",
    "Cl":                              "hydrochloric acid",
    "OB(O)O":                          "boric acid",
    "O=C(O)O":                         "carbonic acid",
    "O=[N+]([O-])[O-]":               "nitrate",
    "O=S(=O)([O-])[O-]":              "sulfate",

    # ── Phase 144: 無機イオン・ヒドロキシルアミン 保留名 ──
    "NO":                              "hydroxylamine",
    "[NH4+]":                          "ammonium",
    "[OH-]":                           "hydroxide",
    "[F-]":                            "fluoride",
    "[Cl-]":                           "chloride",
    "[Br-]":                           "bromide",
    "[I-]":                            "iodide",
    "[Na+]":                           "sodium",
    "[K+]":                            "potassium",
    "[Li+]":                           "lithium",
    "[Ca+2]":                          "calcium",
    "[Mg+2]":                          "magnesium",
    "[Zn+2]":                          "zinc",
    "[Fe+2]":                          "iron(2+)",
    "[Fe+3]":                          "iron(3+)",
    "[Cu+]":                           "copper(1+)",
    "[Cu+2]":                          "copper(2+)",
    "[Al+3]":                          "aluminium",

    # ── Phase 143: 無機リン・ホウ素・ケイ素 保留名 ──
    "P":                               "phosphane",
    "PP":                              "diphosphane",
    "O=P(O)(O)O":                      "phosphoric acid",
    "OP(O)O":                          "phosphorous acid",
    "O=[PH](O)O":                      "hypophosphorous acid",
    "B":                               "borane",
    "BB":                              "diborane",
    "[SiH4]":                          "silane",
    "[SiH3][SiH3]":                    "disilane",

    # ── Phase 237: 追加無機化合物保留名 (IUPAC 2013) ──
    "O=[Se](O)O":                      "selenious acid",
    "O=[Se](=O)(O)O":                  "selenic acid",
    "O=P(O)(O)OP(=O)(O)O":            "pyrophosphoric acid",
    "Cl[Si](Cl)(Cl)Cl":               "tetrachlorosilane",
    "Cl[SiH](Cl)Cl":                  "trichlorosilane",
    "Cl[SiH2]Cl":                     "dichlorosilane",
    "[SiH3]Cl":                        "chlorosilane",
    "F[SiH2]F":                        "difluorosilane",

    # ── Phase 252: 無機リン/イオウハロゲン化物 (IUPAC 2013 P-68) ──
    "ClP(Cl)Cl":                       "phosphorus trichloride",
    "FP(F)F":                          "phosphorus trifluoride",
    "BrP(Br)Br":                       "phosphorus tribromide",
    "IP(I)I":                          "phosphorus triiodide",
    "O=P(Cl)(Cl)Cl":                   "phosphoryl trichloride",
    "O=P(F)(F)F":                      "phosphoryl trifluoride",
    "S=P(Cl)(Cl)Cl":                   "phosphorothioyl trichloride",
    "O=S(Cl)Cl":                       "thionyl chloride",
    "O=S(F)F":                         "thionyl fluoride",
    "O=S(=O)(Cl)Cl":                   "sulfuryl chloride",
    "O=S(=O)(F)F":                     "sulfuryl fluoride",
    # Phase 255: 無機リン窒素化合物
    "NP(N)(N)=O":                      "phosphoric triamide",
    "NP(N)N":                          "phosphorous triamide",

    # ── Phase 258: ジチオカルバミン酸・関連保留名 (IUPAC 2013 P-65.1.1) ──
    "NC(=S)S":                         "carbamodithioic acid",
    "OC(=S)S":                         "carbonodithioic O-acid",

    # ── Phase 163: 炭素・炭素化合物特殊保留名 ──
    "O=C=O":                           "carbon dioxide",
    "[C-]#[O+]":                       "carbon monoxide",
    "S=C=S":                           "carbon disulfide",
    "O=C=S":                           "carbonyl sulfide",
    "ClC(Cl)(Cl)Cl":                   "carbon tetrachloride",  # IUPAC 2013 retained
    "N":                               "ammonia",
    "N#N":                             "dinitrogen",
    "O=O":                             "dioxygen",
    "ClCl":                            "dichlorine",
    "BrBr":                            "dibromine",
    "II":                              "diiodine",
    "FF":                              "difluorine",

    # ── Phase 170: その他単純無機化合物 ──
    "OO":                              "hydrogen peroxide",
    "[H][H]":                          "dihydrogen",
    "S":                               "hydrogen sulfide",
    "O":                               "water",

    # ── Phase 246: デカリン保留名 (IUPAC 2013 P-31.1.3.4) ──
    "C1CCC2CCCCC2C1":                  "decahydronaphthalene",
    "C1CCC2CCCCC2CC1":                 "decahydronaphthalene",  # alt canonical

    # ── Phase 245: 炭素を含まない無機水素化物 (IUPAC 2013 P-68) ──
    "[AsH3]":                          "arsane",
    "[SbH3]":                          "stibane",
    "[BiH3]":                          "bismuthane",
    "[GeH4]":                          "germane",
    "[SnH4]":                          "stannane",
    "[PbH4]":                          "plumbane",

    # ── Phase 176: カルボジイミド・セミカルバジド保留名 ──
    "N=C=N":                           "methanediimine",
    "NNC(N)=O":                        "semicarbazide",
    "NNC(N)=S":                        "thiosemicarbazide",

    # ── Phase 208: シアナミド保留名 (IUPAC 2013 P-66.4.1.1) ──
    "N#CN":                            "cyanamide",

    # ── Phase 228: ホルムアルデヒド保留名 (IUPAC 2013 P-66.6.3.1) ──
    "C=O":                             "formaldehyde",

    # ── Phase 180: アセトアルデヒド・アセトン保留名 (IUPAC 2013 P-31.1.3.4) ──
    "CC=O":                            "acetaldehyde",
    "CC(C)=O":                         "acetone",

    # ── Phase 229: トルエン・酢酸無水物保留名 (IUPAC 2013 P-31.1.3.4, P-65.1.1.3.2) ──
    "Cc1ccccc1":                       "toluene",
    "CC(=O)OC(C)=O":                   "acetic anhydride",

    # ── Phase 268: ベンゾフェノン保留名 (IUPAC 2013 P-31.1.3.4) ──
    "O=C(c1ccccc1)c1ccccc1":           "benzophenone",

    # ── Phase 268: サリチル酸メチル (methyl salicylate, IUPAC 2013 P-65.1.2) ──
    "COC(=O)c1ccccc1O":                "methyl salicylate",

    # ── Phase 240: アセトニトリル・アクリロニトリル保留名 (IUPAC 2013 P-66.6.1.1.1) ──
    "CC#N":                            "acetonitrile",
    "C=CC#N":                          "acrylonitrile",

    # ── Phase 236: 環状イミド保留名 (IUPAC 2013 P-66.8.3) ──
    # succinimide/glutarimide は PIN でないため体系名を使用 (Phase 257)
    "O=C1NC(=O)c2ccccc21":             "phthalimide",

    # ── Phase 230: スチレン保留名 (IUPAC 2013 P-31.1.3.4) ──
    "C=Cc1ccccc1":                     "styrene",

    # ── Phase 181: アクリル酸系・グリコール酸保留名 (IUPAC 2013 P-65.1.1.4) ──
    "C=CC(=O)O":                       "acrylic acid",
    "C=C(C)C(=O)O":                    "methacrylic acid",
    "O=C(O)CO":                        "glycolic acid",

    # ── Phase 182: ケト酸保留名 (IUPAC 2013 P-65.1.1.4) ──
    "CC(=O)C(=O)O":                    "pyruvic acid",
    "CC(=O)CCC(=O)O":                  "levulinic acid",
    "CC(=O)CC(=O)O":                   "acetoacetic acid",

    # ── Phase 179: ハロ酢酸・フェニル酢酸保留名 (IUPAC 2013 P-65.1.1.4) ──
    "O=C(O)CCl":                       "chloroacetic acid",
    "O=C(O)CBr":                       "bromoacetic acid",
    "O=C(O)CF":                        "fluoroacetic acid",
    "O=C(O)CI":                        "iodoacetic acid",
    "O=C(O)C(Cl)(Cl)Cl":              "trichloroacetic acid",
    "O=C(O)C(Cl)Cl":                   "dichloroacetic acid",
    "O=C(O)C(F)(F)F":                  "trifluoroacetic acid",
    "O=C(O)C(F)F":                     "difluoroacetic acid",
    "O=C(O)Cc1ccccc1":                 "phenylacetic acid",

    # ── Phase 259: 芳香族ヒドロキシ酸・トリチオ炭酸 保留名 (IUPAC 2013 P-65.1.1) ──
    # 注: nicotinic acid / anthranilic acid / mandelic acid は PIN でないため除外
    "O=C(O)c1ccccc1O":                "salicylic acid",
    "S=C(S)S":                        "trithiocarbonic acid",

    # ── Phase 259: アミノ酸エステル保留名 (IUPAC 2013 P-65.1.2.4) ──
    # glycinate
    "COC(=O)CN":                      "methyl glycinate",
    "CCOC(=O)CN":                     "ethyl glycinate",
    # alaninate (racemic)
    "COC(=O)C(C)N":                   "methyl alaninate",
    "CCOC(=O)C(C)N":                  "ethyl alaninate",
    # L-alaninate
    "COC(=O)[C@@H](C)N":              "methyl L-alaninate",
    "CCOC(=O)[C@@H](C)N":             "ethyl L-alaninate",
    # D-alaninate
    "COC(=O)[C@H](C)N":               "methyl D-alaninate",
    "CCOC(=O)[C@H](C)N":              "ethyl D-alaninate",
}


def _try_retained_name(graph) -> str | None:
    """canonical SMILES で保留名テーブルを照合する。"""
    if graph.rdkit_mol is None:
        return None
    from rdkit.Chem import MolToSmiles
    canon = MolToSmiles(graph.rdkit_mol)
    return _RETAINED_NAMES.get(canon)


def _name_multicomponent(smiles: str) -> str:
    """複数成分 SMILES (ドット区切り) の命名 — 塩形式に対応  (Phase 147)"""
    from rdkit import Chem
    from collections import Counter
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles!r}")
    frags = Chem.GetMolFrags(mol, asMols=True)
    # 各成分を名前付きリストに
    cat_names: list[str] = []
    ani_names: list[str] = []
    neu_names: list[str] = []
    for frag in frags:
        frag_smi = Chem.MolToSmiles(frag)
        name = smiles_to_iupac(frag_smi)
        charge = sum(a.GetFormalCharge() for a in frag.GetAtoms())
        if charge > 0:
            cat_names.append(name)
        elif charge < 0:
            ani_names.append(name)
        else:
            neu_names.append(name)

    def _compact(names: list[str]) -> str:
        counts = Counter(names)
        from .constants import MULTIPLIER
        parts = []
        for nm in sorted(set(names), key=lambda n: names.index(n)):
            cnt = counts[nm]
            if cnt > 1:
                mult = MULTIPLIER.get(cnt, f"{cnt}")
                parts.append(f"{mult}{nm}")
            else:
                parts.append(nm)
        return " ".join(parts)

    parts = []
    if cat_names:
        parts.append(_compact(cat_names))
    if ani_names:
        parts.append(_compact(ani_names))
    if neu_names:
        parts.append(_compact(neu_names))
    return " ".join(parts)


def _is_isocyanic_acid(graph) -> bool:
    """H-N=C=O パターン (isocyanic acid) を検出する。"""
    from .molecule_analyzer import get_atom, get_bond_order
    for c_idx in range(len(graph.atoms)):
        if get_atom(graph, c_idx).symbol != "C":
            continue
        has_dbl_n_h = False
        has_dbl_o = False
        for nb_idx in graph.adjacency[c_idx]:
            nb = get_atom(graph, nb_idx)
            bo = get_bond_order(graph, c_idx, nb_idx)
            if nb.symbol == "N" and bo == 2.0:
                n_non_ch = [n for n in graph.adjacency[nb_idx]
                            if n != c_idx and get_atom(graph, n).symbol not in ("H",)]
                n_h = sum(1 for n in graph.adjacency[nb_idx]
                          if get_atom(graph, n).symbol == "H")
                if not n_non_ch and n_h >= 1:
                    has_dbl_n_h = True
            elif nb.symbol == "O" and bo == 2.0:
                o_non_ch = [n for n in graph.adjacency[nb_idx]
                            if n != c_idx and get_atom(graph, n).symbol not in ("H",)]
                if not o_non_ch:
                    has_dbl_o = True
        if has_dbl_n_h and has_dbl_o:
            heavy_other = [nb for nb in graph.adjacency[c_idx]
                           if get_atom(graph, nb).symbol not in ("H", "N", "O")]
            if not heavy_other:
                return True
    return False


def smiles_to_iupac(smiles: str) -> str:
    """
    SMILES 文字列から IUPAC 系統名を生成する。

    Args:
        smiles: 有効な SMILES 文字列

    Returns:
        IUPAC 系統名文字列

    Raises:
        ValueError: 無効な SMILES
        NotImplementedError: 未対応構造（縮合多環系など）
    """
    from .name_assembler import fix_enclosing_marks as _fem
    return _fem(_smiles_to_iupac_raw(smiles))


def _smiles_to_iupac_raw(smiles: str) -> str:
    """Internal implementation without enclosing-mark post-processing."""
    # 遅延インポート（rdkit が不要なモジュールを単独テスト可能にする）
    from .molecule_analyzer import build_molecule_graph, get_atom
    from .functional_group import detect_groups, principal_group, FunctionalGroup
    from .chain_finder import find_principal_chain, get_multiple_bond_locants
    from .substituent import collect_substituents
    from .stereochemistry import assign_stereochemistry
    from .name_assembler import assemble_name
    from .ring_handler import (
        has_ring, find_rings, find_principal_ring,
        collect_ring_substituents, assemble_ring_name,
        has_hetero_ring,
    )
    from .heterocycle_handler import name_heterocycle, _try_fused_hetero_retained
    from .polycyclic_handler import name_polycyclic

    # ─── Phase 147: 複数成分 SMILES (塩形式) ─────────────────────────────
    if "." in smiles:
        return _name_multicomponent(smiles)

    # ─── 1. 分子グラフ構築 ───────────────────────────────────────────
    graph = build_molecule_graph(smiles)

    # Phase 135/136: 保留名チェック (早期チェック)
    _aa = _try_retained_name(graph)
    if _aa is not None:
        return _aa

    # ─── 2. 官能基を先行検出（ester/acid_halide は専用パス）───────────
    # Phase 128: イソシアン酸 H-N=C=O (isocyanic acid) の早期検出
    if _is_isocyanic_acid(graph):
        return "isocyanic acid"

    from .functional_group import detect_groups as _detect, principal_group as _pg
    _groups = _detect(graph)
    _pgrp = _pg(_groups)

    # Phase 87: チオ尿素 NC(=S)N (thioamide より先に確認)
    _thiourea = _name_thiourea_if_match(graph, get_atom)
    if _thiourea is not None:
        return _thiourea

    # Phase 297: carbonohydrazide NNC(=O)NN — hydrazide PGRP_DISPATCH より先に確認
    if _is_carbonohydrazide(graph, get_atom):
        return "carbonohydrazide"

    # ─── 2b. PGRP_DISPATCH: group_type → handler ──────────────────────────────
    if _pgrp is not None:
        _handler = PGRP_DISPATCH.get(_pgrp.group_type)
        if _handler is not None:
            _result = _handler(graph, _pgrp, get_atom)
            if _result is not None:
                return _result

    # Special cases not in PGRP_DISPATCH
    if _pgrp is not None and _pgrp.group_type == "sulfamide":
        from .group_namers import _name_n_substituted_sulfamide
        return _name_n_substituted_sulfamide(graph, _pgrp, get_atom)

    # Phase 49: urea 特別ケース (NC(=O)N)
    if _is_urea(graph, get_atom):
        return "urea"

    # Phase 83: N-置換尿素 (CNC(=O)N → N-methylurea 等)
    _sub_urea = _name_substituted_urea_if_match(graph, get_atom)
    if _sub_urea is not None:
        return _sub_urea

    # Phase 165: ヘテロ芳香環 N-オキシド (pyridine 1-oxide 等)
    _het_oxide = _name_hetero_n_oxide(graph, get_atom)
    if _het_oxide is not None:
        return _het_oxide

    # Phase 125: アミン N-オキシド (R₃N→O)
    _n_oxide_name = _name_amine_n_oxide(graph, get_atom)
    if _n_oxide_name is not None:
        return _n_oxide_name

    # Phase 215: ジアゾニウム塩 (C[N+]#N → methanediazonium)
    _diazonium = _name_diazonium(graph, get_atom)
    if _diazonium is not None:
        return _diazonium

    # Phase 225: 亜硫酸エステル (COS(=O)OC → dimethyl sulfite)
    _sulfite_ester = _name_sulfite_ester(graph, get_atom)
    if _sulfite_ester is not None:
        return _sulfite_ester

    # Phase 211: 硫酸エステル・スルファミン酸エステル (COS(=O)(=O)O → methyl hydrogen sulfate)
    _sulfate_ester = _name_sulfate_ester(graph, get_atom)
    if _sulfate_ester is not None:
        return _sulfate_ester

    # Phase 209: アシルペルオキシド (CC(=O)OOC(=O)C → diethanoyl peroxide)
    _acyl_perox = _name_acyl_peroxide(graph, get_atom)
    if _acyl_perox is not None:
        return _acyl_perox

    # Phase 208: N-置換シアナミド (CNC#N → N-methylcyanamide)
    _cyanamide_name = _name_cyanamide(graph, get_atom)
    if _cyanamide_name is not None:
        return _cyanamide_name

    # Phase 202: N-置換ヒドロキシルアミン (CNO → N-methylhydroxylamine)
    _n_hya_name = _name_n_substituted_hydroxylamine(graph, get_atom)
    if _n_hya_name is not None:
        return _n_hya_name

    # Phase 346: O-置換ヒドロキシルアミン (NOC → O-methylhydroxylamine)
    from .group_namers import _name_o_substituted_hydroxylamine
    _o_hya_name = _name_o_substituted_hydroxylamine(graph, get_atom)
    if _o_hya_name is not None:
        return _o_hya_name

    # Phase 347: O-置換オキシム (CC=NOC → O-methylethanal oxime)
    from .group_namers import _name_o_substituted_oxime
    _o_oxime_name = _name_o_substituted_oxime(graph, get_atom)
    if _o_oxime_name is not None:
        return _o_oxime_name

    # Phase 204: ニトロン / イミン N-オキシド (C=[N+]([O-])C → N-methylmethanimine N-oxide)
    _nitrone_name = _name_nitrone(graph, get_atom)
    if _nitrone_name is not None:
        return _nitrone_name

    # Phase 123: ジアゾ化合物 (C=N=N, diazo prefix)
    _diazo_name = _name_diazo_compound(graph, get_atom)
    if _diazo_name is not None:
        return _diazo_name

    # Phase 164: ニトロソアミン (R₂N-N=O)
    _nitrosamine_name = _name_nitrosamine(graph, get_atom)
    if _nitrosamine_name is not None:
        return _nitrosamine_name

    # Phase 113: ヒドラジン (N-N 単結合、環外、カルボニル隣接なし)
    _hydrazine_name = _name_hydrazine_compound(graph, get_atom)
    if _hydrazine_name is not None:
        return _hydrazine_name

    # Phase 115: アゾ化合物 (N=N 二重結合)
    _azo_name = _name_azo_compound(graph, get_atom)
    if _azo_name is not None:
        return _azo_name

    # ─── 3. 環状 or 非環状の分岐 ─────────────────────────────────────
    if has_ring(graph):
        # 3a. スピロ・架橋二環（Phase 15 / Phase 273）
        # acyclic ルート判定より先に実行: 置換スピロ/架橋化合物が acyclic path に
        # 流れてしまうのを防ぐ。非スピロ/非架橋なら None が返るので影響なし。
        poly_name = name_polycyclic(graph)
        if poly_name is not None:
            return poly_name

        # Phase 28/32: アンカー C が環外にある場合 → acyclic path
        # Phase 28: anchor_c1=False (alcohol/amine/ketone 等) でアンカー C が環外
        # Phase 32: anchor_c1=True (carboxylic_acid 等) でアンカー C が環に隣接しない
        if _pgrp is not None and _pgrp.atom_indices:
            from .constants import FUNCTIONAL_GROUPS as _FGS
            _spec = _FGS.get(_pgrp.group_type)
            _anchor_idx = _pgrp.atom_indices[0]
            _anchor_in_ring = get_atom(graph, _anchor_idx).in_ring
            if _spec is not None and not _anchor_in_ring:
                _anchor_adj_ring = any(
                    get_atom(graph, nb).in_ring
                    for nb in graph.adjacency[_anchor_idx]
                )
                # Route to acyclic when:
                # - anchor_c1=False and anchor not in ring (Phase 28)
                # - anchor_c1=True and anchor not adjacent to ring (Phase 32)
                # Phase 40: alkene/alkyne で環の方が大きい場合は ring path を優先
                if not _spec.anchor_c1 or not _anchor_adj_ring:
                    if _pgrp.group_type in ("alkene", "alkyne"):
                        from .molecule_analyzer import non_ring_carbon_indices
                        _nr_c = len(non_ring_carbon_indices(graph))
                        _rings_tmp = find_rings(graph)
                        _max_ring = max((len(r) for r in _rings_tmp), default=0)
                        # ヘテロ芳香族環も考慮 (find_rings は純炭素環のみ)
                        if _max_ring == 0 and graph.ring_atom_sets:
                            _max_ring = max((len(r) for r in graph.ring_atom_sets), default=0)
                        if _max_ring > _nr_c:
                            pass  # ring path へ
                        else:
                            return _name_acyclic(graph, detect_groups, principal_group,
                                                 find_principal_chain, get_multiple_bond_locants,
                                                 collect_substituents, assign_stereochemistry,
                                                 assemble_name, get_atom)
                    else:
                        return _name_acyclic(graph, detect_groups, principal_group,
                                             find_principal_chain, get_multiple_bond_locants,
                                             collect_substituents, assign_stereochemistry,
                                             assemble_name, get_atom)

        # 3b. 縮合環保留名テーブル (Phase 133: 全環系対象、ヘテロ・炭素両方)
        fused_name = _try_fused_hetero_retained(graph)
        if fused_name is not None:
            return fused_name

        # 3c. ヘテロ環保留名チェック（Phase 12）
        if has_hetero_ring(graph):
            hetero_name = name_heterocycle(graph)
            if hetero_name is not None:
                return hetero_name

        return _name_cyclic(graph, find_rings, find_principal_ring,
                            collect_ring_substituents, assemble_ring_name,
                            detect_groups, principal_group, assign_stereochemistry,
                            get_atom)
    else:
        return _name_acyclic(graph, detect_groups, principal_group,
                             find_principal_chain, get_multiple_bond_locants,
                             collect_substituents, assign_stereochemistry,
                             assemble_name, get_atom)



def _name_cyclic(graph, find_rings, find_principal_ring,
                 collect_ring_substituents, assemble_ring_name,
                 detect_groups, principal_group, assign_stereochemistry,
                 get_atom) -> str:
    """環状化合物の命名フロー。"""
    # 官能基検出（環状分子用）
    groups = detect_groups(graph)
    pgrp = principal_group(groups)
    pgrp_type = pgrp.group_type if pgrp is not None else "alkane"

    # 環の検出
    rings = find_rings(graph)
    if not rings:
        # フォールバック: 純炭素環がない場合 (ヘテロ環のみ) は ring_atom_sets を使う
        if graph.ring_atom_sets:
            ring_list = sorted(graph.ring_atom_sets, key=len, reverse=True)
            rings = [list(ring_list[0])]
        else:
            raise ValueError("Ring detection failed despite in_ring=True atoms.")

    # Phase 40: alkene/alkyne で官能基アンカーが環外の場合 → alkane 扱いにして
    # 環外の C=C/C≡C を置換基として収集できるようにする
    if pgrp is not None and pgrp_type in ("alkene", "alkyne"):
        ring_atom_set = set()
        for r in rings:
            ring_atom_set.update(r)
        if not any(ai in ring_atom_set for ai in pgrp.atom_indices):
            pgrp = None
            pgrp_type = "alkane"

    # 主環の選択とロカント割り当て（主官能基アンカーを C1 にする）
    pgrp_atoms_list = pgrp.atom_indices if pgrp is not None else []
    ring_chain = find_principal_ring(graph, rings, pgrp_type, pgrp_atoms_list)

    # suffix ロカント (cyclohexanol, cyclohexanone, cyclohexanamine 等)
    from .constants import FUNCTIONAL_GROUPS as _FG
    suffix_locant: int | None = None
    if pgrp is not None and _FG.get(pgrp.group_type, _FG["alkane"]).needs_locant:
        for ai in pgrp.atom_indices:
            atom = get_atom(graph, ai)
            if atom.symbol == "C" and ai in ring_chain.locant_map:
                suffix_locant = ring_chain.locant_map[ai]
                break
    # 縮合多環 (naphthalene 等) の環外官能基ロカント (naphthalene-2-carboxylic acid 等)
    if (suffix_locant is None and pgrp is not None
            and ring_chain.base_name is not None):
        _spec_ex = _FG.get(pgrp.group_type)
        if _spec_ex is not None and _spec_ex.anchor_c1:
            _anchor_idx = pgrp.atom_indices[0] if pgrp.atom_indices else None
            if _anchor_idx is not None and _anchor_idx not in ring_chain.locant_map:
                for _nb in graph.adjacency[_anchor_idx]:
                    if _nb in ring_chain.locant_map:
                        suffix_locant = ring_chain.locant_map[_nb]
                        break

    # 置換基収集
    pgrp_atoms = pgrp.atom_indices if pgrp is not None else []
    substituents = collect_ring_substituents(graph, ring_chain, pgrp_atoms)

    # 立体化学
    from .chain_finder import PrincipalChain
    dummy_chain = PrincipalChain(
        atom_indices=ring_chain.ring_atoms,
        locant_map=ring_chain.locant_map,
    )
    stereo = assign_stereochemistry(graph, dummy_chain)

    # Phase 27/33: ポリオール / ジオン on ring (benzene 芳香環 or シクロアルカン)
    if pgrp is not None and pgrp_type in ("diol", "triol", "dione", "diamine", "triamine", "dithiol"):
        fg_c_atoms = {ai for ai in pgrp.atom_indices
                      if get_atom(graph, ai).symbol == "C"
                      and ai in ring_chain.locant_map}
        if fg_c_atoms:
            # 全回転 × 2方向で最小ロカントセットを探す
            ring_atoms = ring_chain.ring_atoms
            n_ring = len(ring_atoms)
            best_locs: tuple[int, ...] | None = None
            for start in range(n_ring):
                for direction in (1, -1):
                    if direction == 1:
                        order = [ring_atoms[(start + i) % n_ring] for i in range(n_ring)]
                    else:
                        order = [ring_atoms[(start - i) % n_ring] for i in range(n_ring)]
                    loc_map = {idx: i + 1 for i, idx in enumerate(order)}
                    locs = tuple(sorted(loc_map[a] for a in fg_c_atoms))
                    if best_locs is None or locs < best_locs:
                        best_locs = locs
            loc_str = ",".join(str(l) for l in (best_locs or ()))
            suffix_word = {
                "diol": "diol", "triol": "triol", "dione": "dione",
                "diamine": "diamine", "triamine": "triamine",
                "dithiol": "dithiol",
            }[pgrp_type]
            from .name_assembler import _build_prefix as _bp
            from .constants import CHAIN_PREFIX as _CP
            if ring_chain.is_aromatic and ring_chain.ring_size == 6:
                ring_base = "benzene"
            elif ring_chain.base_name is not None:
                ring_base = ring_chain.base_name
            else:
                _stem = _CP.get(ring_chain.ring_size, f"C{ring_chain.ring_size}")
                _db_locs = ring_chain.double_bond_locants
                if _db_locs:
                    # 環内二重結合を含む: cyclohexa-2,5-diene など
                    _db_str = ",".join(str(l) for l in sorted(_db_locs))
                    _num_db = len(_db_locs)
                    if _num_db == 1:
                        ring_base = f"cyclo{_stem}-{_db_str}-ene"
                    elif _num_db == 2:
                        ring_base = f"cyclo{_stem}a-{_db_str}-diene"
                    else:
                        _mult = {3: "triene", 4: "tetraene"}.get(_num_db, f"{_num_db}ene")
                        ring_base = f"cyclo{_stem}a-{_db_str}-{_mult}"
                else:
                    ring_base = f"cyclo{_stem}ane"
            prefix_part = _bp(substituents) if substituents else ""
            base_name = f"{ring_base}-{loc_str}-{suffix_word}"
            # 立体記述子 (Phase 66: ring diol/dione 立体化学)
            if stereo:
                combined = ",".join(d.strip("()") for d in stereo)
                stereo_str = f"({combined})-"
            else:
                stereo_str = ""
            full = f"{prefix_part}{base_name}" if prefix_part else base_name
            return f"{stereo_str}{full}"

    # N-置換アミド on 環 (benzamide → N-methylbenzamide など) の処理
    if pgrp is not None and pgrp_type == "amide":
        n_atom_idx = next(
            (ai for ai in pgrp.atom_indices if get_atom(graph, ai).symbol == "N"), None
        )
        if n_atom_idx is not None:
            h_count = sum(
                1 for nb in graph.adjacency[n_atom_idx]
                if get_atom(graph, nb).symbol == "H"
            )
            if h_count < 2:
                from .constants import MULTIPLIER
                from .substituent import _name_carbon_substituent
                from collections import Counter as _Counter
                carbonyl_c = next(
                    (ai for ai in pgrp.atom_indices if get_atom(graph, ai).symbol == "C"), None
                )
                base_amide_name = assemble_ring_name(
                    ring_chain=ring_chain,
                    substituents=substituents,
                    principal_grp_type=pgrp_type,
                    suffix_locant=suffix_locant,
                    stereo_descriptors=stereo,
                )
                excluded_n = {n_atom_idx}
                if carbonyl_c is not None:
                    excluded_n.add(carbonyl_c)
                c_nbrs_n = [
                    nb for nb in graph.adjacency[n_atom_idx]
                    if nb not in excluded_n and get_atom(graph, nb).symbol == "C"
                ]
                o_nbrs_n = [
                    nb for nb in graph.adjacency[n_atom_idx]
                    if get_atom(graph, nb).symbol == "O"
                ]
                n_hydroxy = [
                    "hydroxy" for o in o_nbrs_n
                    if any(get_atom(graph, nb2).symbol == "H"
                           for nb2 in graph.adjacency[o])
                ]
                if not c_nbrs_n and not n_hydroxy:
                    return base_amide_name
                from .name_assembler import _needs_bis_tris as _nbp2
                n_subs = [_name_carbon_substituent(graph, c, {n_atom_idx}) for c in c_nbrs_n]
                n_subs += n_hydroxy
                sub_counts = _Counter(n_subs)
                prefix_parts = []
                for sub in sorted(sub_counts):
                    cnt = sub_counts[sub]
                    sub_str = f"({sub})" if _nbp2(sub) else sub
                    if cnt == 1:
                        prefix_parts.append(f"N-{sub_str}")
                    else:
                        mult = MULTIPLIER.get(cnt, f"{cnt}")
                        prefix_parts.append(f"N,N-{mult}{sub_str}")
                prefix = "-".join(prefix_parts)
                return f"{prefix}{base_amide_name}"

    # 名前組み立て
    return assemble_ring_name(
        ring_chain=ring_chain,
        substituents=substituents,
        principal_grp_type=pgrp_type,
        suffix_locant=suffix_locant,
        stereo_descriptors=stereo,
    )



def _name_acyclic(graph, detect_groups, principal_group,
                  find_principal_chain, get_multiple_bond_locants,
                  collect_substituents, assign_stereochemistry,
                  assemble_name, get_atom) -> str:
    """非環状化合物の命名フロー（Phase 1-2）。"""
    from .functional_group import FunctionalGroup

    # 官能基検出
    groups = detect_groups(graph)
    pgrp: FunctionalGroup | None = principal_group(groups)

    # 二級・三級アミドは専用パスで処理
    if pgrp is not None and pgrp.group_type == "amide":
        n_atom_idx = next((ai for ai in pgrp.atom_indices if get_atom(graph, ai).symbol == "N"), None)
        if n_atom_idx is not None:
            h_count = sum(1 for nb in graph.adjacency[n_atom_idx]
                          if get_atom(graph, nb).symbol == "H")
            if h_count < 2:
                carbonyl_c = pgrp.atom_indices[0]
                return _name_secondary_tertiary_amide(graph, carbonyl_c, n_atom_idx, get_atom)

    # 二級・三級アミンは専用パスで処理
    if pgrp is not None and pgrp.group_type == "amine":
        if get_atom(graph, pgrp.atom_indices[0]).symbol == "N":
            n_idx = pgrp.atom_indices[0]
            c_nbrs = pgrp.atom_indices[1:]
            return _name_secondary_tertiary_amine(graph, n_idx, c_nbrs, get_atom)

    # アミジン (N-置換含む) は専用パスで処理 (Phase 185)
    if pgrp is not None and pgrp.group_type == "amidine":
        from .functional_group import _get_amidine_nitrogens
        _amidine_c = pgrp.atom_indices[0]
        _n_imine, _n_amine = _get_amidine_nitrogens(graph, _amidine_c)
        return _name_substituted_amidine(graph, _amidine_c, _n_imine, _n_amine, get_atom)

    # N-置換イミン (Phase 64)
    if pgrp is not None and pgrp.group_type == "imine":
        n_atom_idx = next((ai for ai in pgrp.atom_indices if get_atom(graph, ai).symbol == "N"), None)
        if n_atom_idx is not None:
            imine_c = pgrp.atom_indices[0]
            n_c_subs = [nb for nb in graph.adjacency[n_atom_idx]
                        if nb != imine_c and get_atom(graph, nb).symbol == "C"]
            if n_c_subs:
                return _name_n_substituted_imine(
                    graph, imine_c, n_atom_idx, n_c_subs,
                    find_principal_chain, get_multiple_bond_locants,
                    collect_substituents, assemble_name, get_atom,
                )

    # 主鎖探索・ロカント割り当て
    chain = find_principal_chain(graph, pgrp)

    # 多重結合ロカント
    mb = get_multiple_bond_locants(graph, chain)

    # suffix ロカント（複数官能基対応）
    from .constants import FUNCTIONAL_GROUPS as _FG
    suffix_locants_list: list[int] = []
    if pgrp is not None and _FG.get(pgrp.group_type, _FG["alkane"]).needs_locant:
        if pgrp.group_type in ("diol", "triol"):
            # ジェミナルジオール対応: 各 O から親 C のロカントを収集 (重複許可)
            for ai in pgrp.atom_indices:
                atom = get_atom(graph, ai)
                if atom.symbol == "O":
                    for nb in graph.adjacency[ai]:
                        if get_atom(graph, nb).symbol == "C" and nb in chain.locant_map:
                            suffix_locants_list.append(chain.locant_map[nb])
                            break
        elif pgrp.group_type == "dithiol":
            # ジチオール: 各 S から親 C のロカントを収集
            for ai in pgrp.atom_indices:
                atom = get_atom(graph, ai)
                if atom.symbol == "S":
                    for nb in graph.adjacency[ai]:
                        if get_atom(graph, nb).symbol == "C" and nb in chain.locant_map:
                            suffix_locants_list.append(chain.locant_map[nb])
                            break
        else:
            for ai in pgrp.atom_indices:
                atom = get_atom(graph, ai)
                if atom.symbol == "C" and ai in chain.locant_map:
                    loc = chain.locant_map[ai]
                    if loc not in suffix_locants_list:
                        suffix_locants_list.append(loc)
        suffix_locants_list.sort()

    suffix_locant = suffix_locants_list[0] if suffix_locants_list else None
    suffix_locants = suffix_locants_list if len(suffix_locants_list) > 1 else None

    # 置換基収集
    pgrp_atom_indices = pgrp.atom_indices if pgrp is not None else []
    substituents = collect_substituents(
        graph, chain.atom_indices, chain.locant_map, pgrp_atom_indices,
    )

    # 立体化学
    stereo = assign_stereochemistry(graph, chain)

    # Phase 349: C=N 結合の E/Z 記述子 (oxime, imine, hydrazone, semicarbazone 等)
    _cn_stereo_groups = (
        "aldoxime", "ketoxime",
        "aldhydrazone", "kethydrazone",
        "aldsemicarbazone", "semicarbazone",
        "aldthiosemicarbazone", "thiosemicarbazone",
        "imine", "diimine",
    )
    if pgrp is not None and pgrp.group_type in _cn_stereo_groups:
        from .stereochemistry import _get_bond_stereo
        _cn_c = pgrp.atom_indices[0]
        _cn_ns = [ai for ai in pgrp.atom_indices[1:]
                  if get_atom(graph, ai).symbol == "N"]
        if _cn_ns:
            _cn_stereo = _get_bond_stereo(graph, _cn_c, _cn_ns[0])
            if _cn_stereo is not None:
                stereo = list(stereo) + [f"({_cn_stereo})"]

    # 名前組み立て
    pgrp_type = pgrp.group_type if pgrp is not None else "alkane"
    return assemble_name(
        chain_length=chain.length,
        principal_group_type=pgrp_type,
        multiple_bonds=mb,
        substituents=substituents,
        stereo_descriptors=stereo,
        suffix_locant=suffix_locant,
        suffix_locants=suffix_locants,
    )
