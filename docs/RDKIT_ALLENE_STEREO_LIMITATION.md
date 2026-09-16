# RDKit におけるアレン軸性キラリティ (IUPAC P-92.4) の SMILES 規格適合性調査・バグ性検証・修正設計レポート

**文書作成日**: 2026-08-26  
**対象規格**: 
- **OpenSMILES Specification §3.8.4 (Tetrahedral Allene-like Systems)**
- **Daylight SMILES Theory Manual §Stereochemistry (Chiral Class: AL)**
- **IUPAC 2013 Blue Book Rule P-92.4 (Axial Chirality)**
**対象ライブラリ**: RDKit (`rdkit.Chem` / C++ Core)

---

## 1. エグゼクティブサマリー (Executive Summary)

### 1.1 「この SMILES 記法は正当な仕様か？」に対する結論
**結論: 100% 正当な仕様（Standard-compliant）です。**

アレン中心炭素に `@` / `@@` を付与する記法（例: `NC(Br)=[C@]=C(O)C` や `C/C=[C@]=C/C`）は、**Daylight SMILES（元祖規格）および OpenSMILES（オープンコミュニティ標準規格）の両方の公式仕様書に明記されている正式な構文** です。独自解釈や非標準記法ではありません。

### 1.2 RDKit の挙動評価
- **規格適合性の観点**: **OpenSMILES 規格に対するパーサーおよびデータモデルの欠陥（Bug / Incomplete Implementation）** です。
- **RDKit チームへの PR の妥当性**: 非常に正当であり、化学情報学コミュニティにとっても長年望まれている重要な機能追加です。ただし、C++ コアの広範なモジュール（パーサー、キラリティ判定器、シリアライザー）にまたがるため、**「OpenSMILES 規格準拠に向けた拡張四面体立体サポート」** として Issue を起票し、設計をすり合わせた上で PR を提出するのが最も建設的です。

---

## 2. 公的規格におけるアレン SMILES の厳密な仕様定義

### 2.1 OpenSMILES 規格 (§3.8.4)
> **OpenSMILES Specification §3.8.4: Tetrahedral Allene-like Systems**  
> *"Extended tetrahedral configurations can be specified for conjugated allenes with an even number of double bonds. The normal tetrahedral rules using `'@'` and `'@@'` apply, but the "neighbor" atoms to which the chirality refers are at the ends of the allene system. For example: `NC(Br)=[C@]=C(O)C`"*  
> *"To determine the correct clockwise or anticlockwise specification, the allene is conceptually 'collapsed' into a single tetrahedral chiral center, and the resulting chirality is marked as a property of the center atom of the extended allene system."*

#### 仮想四面体凝縮モデル (Collapsed Tetrahedron Model)
OpenSMILES では、アレン $A(B)\text{C}=\text{C}=\text{C}(D)E$ を次のように解釈します：
1. 2つの二重結合を仮想的に長さゼロに縮約（Collapse）し、両端の 4 つの置換基 $A, B, D, E$ が中心炭素に直接結合した「仮想的な四面体中心」とみなす。
2. SMILES 文字列の走査順で第1の隣接置換基から中心原子を見通したとき、残りの 3 置換基の出現順序が：
   - **反時計回り (Anticlockwise)**: `[C@]`
   - **時計回り (Clockwise)**: `[C@@]`

```
   (手前)                    (仮想四面体)
     A                          A
      \                        / \
       C = C = C      ===>    / C \
      /        \             B     D
     B          D             \   /
                 \              E
                  E           (奥)
```

---

### 2.2 Daylight SMILES 規格 (Theory Manual §Stereochemistry)
Daylight のオリジナル仕様でも、アレンは正式なキラルクラス **`AL` (Allene)** として定義されています。

- 完全形: `OC(Cl)=[C@AL1]=C(C)F` （反時計回り） / `OC(Cl)=[C@AL2]=C(C)F` （時計回り）
- 短縮形: `OC(Cl)=[C@]=C(C)F` （`@AL1` のデフォルト省略形として `@` が認められている）

したがって、`[C@]` / `[C@@]` をアレン中心に配置する記法は、30年以上の歴史を持つ SMILES の確立された正当な文法です。

---

## 3. RDKit 内部で情報が消失するメカニズムの特定

RDKit 2026.03 の C++ ソースコードを追跡した結果、以下のフローで立体情報が破棄されていることが判明しました：

```
[SMILES入力: "CC=[C@]=CC"]
       │
       ▼
①【Bison パーサー】(Code/GraphMol/SmilesParse/smiles.yy)
   ・"[C@]" トークンを受理し、Atom オブジェクトに ChiralTag = CHI_TETRAHEDRAL_CCW を設定。
   ・（ここでは正常に保持されている）
       │
       ▼
②【サニタイズ処理】(Code/GraphMol/Sanitizer.cpp -> cleanUp())
   ・各原子の価数、芳香属性、混成軌道（sp/sp2/sp3）を計算。
       │
       ▼
③【立体化学割り当て】(Code/GraphMol/Chirality.cpp -> assignStereochemistry())
   ・各原子の配位数（Degree）と混成軌道を検査。
   ・sp 炭素（重原子配位数 = 2）は 4 配位四面体ではないため、
     「無効な ChiralTag」と判定されて CHI_UNSPECIFIED にリセットされる。
       │
       ▼
[最終 Mol オブジェクト] -> Atom(2).GetChiralTag() == CHI_UNSPECIFIED (消失)
```

### RDKit 内の既存コードの痕跡
興味深いことに、RDKit の `Code/GraphMol/Atom.h` の `ChiralType` enum には以下のような定義が存在します：
```cpp
enum ChiralType {
  CHI_UNSPECIFIED = 0,
  CHI_TETRAHEDRAL_CW,
  CHI_TETRAHEDRAL_CCW,
  CHI_OTHER,
  CHI_TETRAHEDRAL,
  CHI_ALLENE,          // ← 開発初期に定義されたが、実装が未完のまま放置されている
  CHI_SQUAREPLANAR,
  ...
};
```
つまり、RDKit 開発チーム自身もアレン対応の必要性を認識して enum を用意していたものの、四面体不斉（`CHI_TETRAHEDRAL`）の実装完了後、アレンの割り当てロジック（`assignStereochemistry`）が未実装のまま今日に至っているという経緯が読み取れます。

---

## 4. IUPAC 2013 P-92.4 の CIP (Ra)/(Sa) への決定論的変換

SMILES の `@`/`@@` から IUPAC 2013 推奨名 (PIN) の $(Ra)/(Sa)$ への変換は、完全に数学的・決定論的に定義可能です。

### 変換アルゴリズム
1. **アレン骨格の同定**: $C_L = C_C = C_R$
2. **端点置換基の CIP 優先順位**:
   - 手前の端点 $C_L$ の置換基: 上位 $L_1$、下位 $L_2$（CIP 優先度 $L_1 > L_2$）
   - 奥の端点 $C_R$ の置換基: 上位 $R_1$、下位 $R_2$（CIP 優先度 $R_1 > R_2$）
   - IUPAC P-92.4 の軸性 CIP 順位: **$L_1 (1) > L_2 (2) > R_1 (3) > R_2 (4)$**（手前 > 奥）
3. **SMILES `@`/`@@` の置換基走査順序とのパリティ計算**:
   - SMILES の出現順序で並べた置換基列 $(S_1, S_2, S_3, S_4)$ に対し、CIP 順位列 $(L_1, L_2, R_1, R_2)$ への置換の偶奇（Permutation Sign）を計算。
   - `@`（反時計回り）× 偶置換 $\implies$ **$(Ra)$** (または $(P)$)
   - `@@`（時計回り）× 偶置換 $\implies$ **$(Sa)$** (または $(M)$)

この変換式は普遍的であり、RDKit 本体の C++ 改修でも、`smiles2iupac` 内の Python モジュールでも全く同一のロジックで実装できます。

---

## 5. RDKit への Pull Request / Issue 起票の設計提案

RDKit チームに迷惑をかけず、最も歓迎される形でコントリビューションを行うためのロードマップです。

### 5.1 推奨される Issue 起票内容 (Draft)
- **Title**: `[Feature Request / Standard Compliance] Support OpenSMILES §3.8.4 extended tetrahedral (allene) chirality`
- **Summary**:
  - OpenSMILES §3.8.4 defines allene chirality on central sp carbon using `@`/`@@` (e.g. `NC(Br)=[C@]=C(O)C`).
  - Currently, `MolFromSmiles` silently clears `ChiralTag` on 2-coordinate carbons in `assignStereochemistry`.
  - Propose activating `CHI_ALLENE` or preserving extended tetrahedral stereo for cumulated double bonds ($C=C=C$).
- **Reference**: Link to OpenSMILES §3.8.4 and Daylight Theory Manual §Stereochemistry.

### 5.2 C++ 実装のパッチ対象
1. **`Code/GraphMol/Chirality.cpp`**:
   - `assignStereochemistry()`: 2 配位炭素であっても、両隣が二重結合（`C=C=C`）である場合は `CHI_ALLENE`（または `CHI_TETRAHEDRAL`）を保持し、両端の 4 置換基の相対立体配置を計算。
2. **`Code/GraphMol/CIPLabeler/`**:
   - `CIPLabeler.cpp`: 軸性キラリティ（Extended Tetrahedral CIP rule, IUPAC P-92.4）の `Ra` / `Sa` ラベル付けメソッドの追加。
3. **`Code/GraphMol/SmilesParse/SmilesWrite.cpp`**:
   - アレン中心炭素への `[C@]` / `[C@@]` 出力のサポート。

---

## 6. smiles2iupac における当面の対応方針

1. **RDKit への Issue/PR 提案**: 上記の設計に基づき、RDKit コミュニティへ正式に Issue を起票・提案可能。
2. **smiles2iupac での自立的先行実装**:
   - ユーザーが `smiles_to_iupac("C/C(Cl)=[C@]=C(\\F)C")` を入力した際、RDKit の古いバージョンでも動作するよう、`smiles2iupac` 内部の `stereo_perception.py` で SMILES プリパース＋独自 CIP 判定を行い、即座に `(3Ra)-4-chloro-2-fluoropenta-2,3-diene` を出力できる自立構造を整える。
