# ADR-007: IUPAC 未実装規則の実装ロードマップ

**Status:** Proposed  
**Date:** 2026-05-22  
**Deciders:** smiles2iupac 開発チーム  
**Reference:** IUPAC 2013 Blue Book (Nomenclature of Organic Chemistry), IUPAC Brief Guide to Nomenclature (2020)

---

## Context

smiles2iupac は Phase 1–6 を経て、炭素系化合物（アルカン〜三環縮合芳香族）の命名を実装した。
しかし IUPAC 2013 Blue Book が定義する有機化合物命名規則のうち、窒素・酸素・硫黄ヘテロ原子を含む官能基、エステル・酸ハライド、複数官能基の同時処理、ヘテロ環、架橋多環系など多くの規則が未実装のまま残っている。

本 ADR は IUPAC Brief Guide の seniority order（Table 3）と Blue Book の各章を照合して不足規則を列挙し、実装優先度と工数見積を示すロードマップを定める。

---

## 現在の実装状況

| Phase | 対象 | 状態 |
|-------|------|------|
| 1 | 直鎖・分岐アルカン、アルコール、ハロゲン化物、アルケン・アルキン | ✅ 完了 |
| 2 | カルボン酸、アルデヒド、ケトン、立体化学 (R/S, E/Z) | ✅ 完了 |
| 3 | シクロアルカン、ベンゼン、置換ベンゼン | ✅ 完了 |
| 4 | ベンゼン官能基誘導体 (phenol, benzoic acid, benzaldehyde) | ✅ 完了 |
| 5 | ナフタレン / 置換ナフタレン | ✅ 完了 |
| 6 | 三環縮合多環芳香族 (anthracene, phenanthrene) | ✅ 完了 |

---

## ギャップ分析

IUPAC Brief Guide Table 3 の seniority order（優先順位降順）と現実装を対比する。

| 優先順位 | 官能基クラス | 構造式 | Suffix | Prefix | 実装状況 |
|----------|------------|--------|--------|--------|----------|
| 1 | Carboxylic acids | –COOH | -oic acid / -carboxylic acid | carboxy | ✅ 一部 |
| 2 | Esters | –COOR | -(R)...oate | (R)oxycarbonyl | ❌ 未実装 |
| 3 | Acid halides | –COX | -oyl halide | halocarbonyl | ❌ 未実装 |
| 4 | Amides | –CONH₂ | -amide / -carboxamide | carbamoyl | ❌ 未実装 |
| 5 | Nitriles | –C≡N | -nitrile / -carbonitrile | cyano | ❌ 未実装 |
| 6 | Aldehydes | –CHO | -al / -carbaldehyde | formyl / oxo | ✅ 一部※1 |
| 7 | Ketones | =O | -one | oxo | ✅ |
| 8 | Alcohols | –OH | -ol | hydroxy | ✅ |
| 9 | Thiols | –SH | -thiol | sulfanyl | ❌ 未実装 |
| 10 | Amines | –NH₂ | -amine | amino | ❌ 未実装 |
| 11 | Imines | =NH | -imine | imino | ❌ 未実装 |
| — | Ethers | –OR | — | (R)oxy | ❌ 未実装 |
| — | Nitro | –NO₂ | — | nitro | ❌ 未実装 |

※1 アルデヒドは `propan-1-al` 形式は実装済み。環付加 (–CHO が環外) の `cyclohexanecarbaldehyde` 形式は未対応。

### 未実装の構造型

| 種別 | 例 | Blue Book 参照 | 実装状況 |
|------|----|--------------|----------|
| 二酸 (dioic acid) | propanedioic acid, butanedioic acid | P-65.1.1 | ❌ |
| ポリオール (diol 等) | ethane-1,2-diol | P-63.1.3 | ❌ |
| 複合置換基の bis/tris | bis(2-chloroethyl) → bis でなく di を使用中 | P-14.5 | ❌ (部分) |
| ヘテロ環 — 保留名 | pyridine, furan, thiophene, pyrrole | P-31.1.3 | ❌ |
| ヘテロ環 — Hantzsch-Widman | aziridine, oxirane, azetidine | P-31.1.3.4 | ❌ |
| 架橋多環 (von Baeyer) | bicyclo[2.2.1]heptane | P-31.1.2 | ❌ |
| スピロ化合物 | spiro[4.5]decane | P-31.1.1 | ❌ |
| 4 環以上の縮合多環 | pyrene, chrysene | P-31.1.3 | ❌ |
| シクロアルケン | cyclohex-1-ene | P-31.1.3 | ❌ |
| 部分飽和環の indicated hydrogen | 1H-indene | P-31.1.2 | ❌ |

---

## 提案する実装フェーズ

### Phase 7 — 含窒素官能基（優先度: 高）

**対象規則:** IUPAC P-62 (amines), P-66.3 (nitriles)

**実装内容:**
- Primary amines: `–NH₂` → suffix `-amine`（例: methanamine, ethanamine, propan-2-amine）
- Nitriles: `–C≡N` → suffix `-nitrile` (鎖末端) / `-carbonitrile`（環付加）（例: propanenitrile, benzonitrile）
- Nitro: `–NO₂` → prefix `nitro`（接頭辞のみ、例: nitrobenzene, 1-nitropropane）

**変更ファイル:**
- `functional_group.py`: `detect_groups()` に `amine`, `nitrile`, `nitro` の検出を追加。seniority: amine=40, nitrile=80
- `constants.py`: `FUNCTIONAL_GROUP_PRIORITY`, `SUFFIX_MAP` に追加
- `substituent.py`: `name_substituent()` に NH₂ → `amino`, N≡ → `cyano`, NO₂ → `nitro` の処理を追加
- `name_assembler.py`: `_build_name_body()` に `-amine` / `-nitrile` の suffix 処理を追加

**代表テストケース:**

```
CC(N)C          → propan-2-amine
CCN             → ethanamine
CCC#N           → propanenitrile
c1ccc(N)cc1     → aniline  (保留名; PIN = benzenamine)
c1ccc([N+](=O)[O-])cc1 → nitrobenzene
```

**工数見積:** S（2–3日）  
**依存関係:** なし（Phase 1–2 の infrastructure を再利用）

---

### Phase 8 — チオール・エーテル（優先度: 高）

**対象規則:** IUPAC P-63.6 (thiols), P-63.4 (ethers)

**実装内容:**
- Thiols: `–SH` → suffix `-thiol`（例: ethanethiol, propane-1-thiol）
- Ethers: `–OR` → prefix `(alkyl)oxy`（例: methoxybenzene, 1-ethoxypropane）
- Sulfides: `–SR` → prefix `(alkyl)sulfanyl`（例: methylsulfanylbenzene）

**変更ファイル:**
- `functional_group.py`: thiol 検出（O → S に類比）。seniority: thiol=45
- `substituent.py`: S に H → `sulfanyl`, S に C → `(R)sulfanyl`, O に C（非 OH）→ `(R)oxy`
- `name_assembler.py`: `-thiol` suffix を追加

**代表テストケース:**

```
CCS             → ethanethiol
CCOC            → methoxyethane   (ethyl methyl ether の PIN)
c1ccccc1OC      → methoxybenzene  (anisole の PIN = methoxybenzene)
```

**工数見積:** S（2–3日）  
**依存関係:** なし

---

### Phase 9 — 複数官能基・ポリオール・二酸（優先度: 高）

**対象規則:** IUPAC P-14.5 (multiplying prefixes), P-63.1.3 (diols), P-65.1.1 (dioic acids)

**実装内容:**

**(a) bis/tris 接頭辞** — 複合置換基（ロカント付きまたは括弧必要な名前）には `di/tri` でなく `bis/tris` を使う。

```
判定規則: 置換基名に数字・括弧・倍数詞が含まれる → bis/tris
例: 2 × (1-methylethyl) → bis(1-methylethyl)  [not di(1-methylethyl)]
```

**(b) 二酸 (dioic acid)** — 鎖両端が `–COOH` → suffix `-dioic acid`。

```
OC(=O)CC(=O)O   → propanedioic acid
OC(=O)CCC(=O)O  → butanedioic acid
```

**(c) ポリオール / ジケトン** — `–OH` / `C=O` が複数 → suffix `-diol` / `-dione`。

```
OCC O            → ethane-1,2-diol
CC(=O)CC(=O)C    → pentane-2,4-dione
```

**変更ファイル:**
- `functional_group.py`: 複数の同一官能基を group_type に `_count` 付きで返す、または groups リストの後処理で多重検出
- `name_assembler.py`: `_build_name_body()` に `-dioic acid`, `-diol`, `-dione` の suffix を追加; `_build_prefix()` に bis/tris 判定ロジックを追加
- `constants.py`: `MULTIPLIER` に bis/tris 用の `BIS_MULTIPLIER` を追加

**工数見積:** M（4–5日）  
**依存関係:** Phase 7, 8 の後でも前でも実装可能

---

### Phase 10 — エステル・酸ハライド（優先度: 中）

**対象規則:** IUPAC P-65.1.2 (esters), P-65.1.3 (acid halides); functional class names 優先

**実装内容:**
- Esters: `–COOR` → functional class name `(R) (stem)oate`（例: ethyl ethanoate, methyl propanoate）
- Acid halides: `–COX` → functional class name `(stem)oyl (halide)`（例: ethanoyl chloride）

**SMILES 例:**
```
CCOC(=O)C   → ethyl ethanoate      (C に =O と OCC)
CC(=O)Cl    → ethanoyl chloride
CC(=O)F     → ethanoyl fluoride
```

**設計の注意:**
- エステルの SMILES は `O=C(OCC)C` または `CCOC(=O)C`。R 基（アルキルオキシ側）の命名に `name_acyclic_chain()` を再利用する。
- IUPAC 2013 は `ethyl ethanoate` を PIN、`ethyl acetate` を retained として許容。ここでは PIN を優先する。

**変更ファイル:**
- `functional_group.py`: ester 検出（`–C(=O)O–C`）, acid_halide 検出（`–C(=O)X`）; seniority: ester=90, acid_halide=85
- `__init__.py`: ester/acid_halide は functional class name なので `_name_acyclic` のロジックを一部変更（R 基 + 親酸名 の組み立て）

**工数見積:** M（4–5日）  
**依存関係:** Phase 7（amine seniority との相互作用確認）

---

### Phase 11 — アミド・イミン（優先度: 中）

**対象規則:** IUPAC P-66.1 (amides), P-66.1.5 (imines)

**実装内容:**
- Primary amides: `–CONH₂` → suffix `-amide` (鎖末端) / `-carboxamide`（環付加）
- Imines: `=NH` → suffix `-imine`, prefix `imino`

```
CC(=O)N    → acetamide  (保留名; PIN = ethanamide)
CC(N)=O    → ethanamide
```

**工数見積:** S（2–3日）  
**依存関係:** Phase 7（amine との seniority 判定が必要）

---

### Phase 12 — ヘテロ環 保留名（優先度: 中）

**対象規則:** IUPAC P-31.1.3, Table 2 (retained names)

最も使用頻度が高い保留名ヘテロ環を lookup table で先行実装する。

**実装する保留名（優先 20 種）:**

| SMILES パターン | 名称 | 環サイズ | ヘテロ原子 |
|---------------|------|---------|-----------|
| `c1ccncc1` | pyridine | 6 | N |
| `c1ccoc1` | furan | 5 | O |
| `c1ccsc1` | thiophene | 5 | S |
| `c1cc[nH]c1` | pyrrole | 5 | NH |
| `c1cnc[nH]1` | imidazole | 5 | N, NH |
| `c1ccnn1` | pyrazole | 5 | N, NH |
| `c1ccnc1` | pyridine (alias) | — | — |
| `c1ccncc1` | pyridine | 6 | N |
| `c1ccncn1` | pyrimidine | 6 | 2N |
| `c1cnccn1` | pyrazine | 6 | 2N |
| `c1cnncc1` | pyridazine | 6 | 2N |
| `C1CCNCC1` | piperidine | 6 | N |
| `C1CCOCC1` | oxane (tetrahydropyran) | 6 | O |
| `C1CCNCO1` | morpholine | 6 | N, O |
| `c1ccc2ncccc2c1` | quinoline | 10 | N |
| `c1ccc2cnccc2c1` | isoquinoline | 10 | N |
| `c1cc2ccccc2[nH]1` | indole | 10 | NH |
| `c1csc2ccccc12` | benzothiophene | 10 | S |
| `c1cc2cccc c2o1` | benzofuran | 10 | O |
| `C1CCCC1` | (cyclopentane; 既存) | — | — |

**変更ファイル:**
- `ring_handler.py`: `find_rings()` を拡張してヘテロ原子含有環を許容（現在は炭素環のみ）
- 新規ファイル `heterocycle_handler.py`: ヘテロ環の検出・保留名照合・置換基命名を担当
- `molecule_analyzer.py`: `find_rings` の炭素限定フィルタを調整

**工数見積:** L（7–10日）  
**依存関係:** Phase 7, 8（heteroatom 検出基盤）

---

### Phase 13 — Hantzsch-Widman ヘテロ環（優先度: 低〜中）

**対象規則:** IUPAC P-31.1.3.4, Table 1 (Hantzsch-Widman names)

保留名にない小環ヘテロ環を Hantzsch-Widman 規則で自動命名する。

| 環サイズ | O | S | N |
|---------|---|---|---|
| 3 | oxirane | thiirane | aziridine |
| 4 | oxetane | thietane | azetidine |
| 5 (飽和) | oxolane (THF) | thiolane | pyrrolidine |
| 5 (不飽和) | — | — | 3H-pyrrole 等 |
| 6 (飽和) | oxane | thiane | piperidine |

**構成:** ヘテロ原子前置詞（oxa, thia, aza）＋ 環サイズ語幹（ir, et, ol, an）＋ ine/ane

**工数見積:** M（4–5日）  
**依存関係:** Phase 12

---

### Phase 14 — シクロアルケン（優先度: 中）

**対象規則:** IUPAC P-31.1.3 (cycloolefins)

**実装内容:**
- 環内の C=C 二重結合を検出し、ロカントを付与する。
- 例: `C1=CCCCC1` → `cyclohex-1-ene`
- 複数二重結合: `C1=CC=CCC1` → `cyclohexa-1,3-diene`

**変更ファイル:**
- `ring_handler.py`: `_assign_ring_locants()` 内でアルケン結合を検出し `ring_size` に対して二重結合ロカントを計算
- `functional_group.py`: `detect_groups()` のアルケン検出フィルタ（現在 `atom.in_ring` のものは除外）を緩和

**工数見積:** S（2–3日）  
**依存関係:** Phase 3 実装（cycloalkane）に隣接

---

### Phase 15 — 架橋多環系・スピロ化合物（優先度: 低）

**対象規則:** IUPAC P-31.1.1 (spiro), P-31.1.2 (bridged bicyclic)

**実装内容:**

**(a) スピロ化合物** — 1つの原子（スピロ原子）で2環が結合。  
命名: `spiro[小環サイズ.大環サイズ]アルカン`  
例: `C1CCC2(CC1)CCCC2` → `spiro[4.5]decane`

**(b) 架橋二環 (von Baeyer)** — 2つの橋頭原子を持つ二環。  
命名: `bicyclo[橋1.橋2.橋3]アルカン`（橋の炭素数を降順に）  
例: `C1CC2CCC1CC2` → `bicyclo[2.2.2]octane`  
例: `C1CC2CCC1C2` → `bicyclo[2.2.1]heptane` (norbornane)

**変更ファイル:**
- `ring_handler.py`: `has_ring()` の分岐の前に spiro / bridged bicyclic を検出するルーティングを追加
- 新規ファイル `polycyclic_handler.py`: von Baeyer 命名・スピロ命名を実装

**工数見積:** L（8–12日）  
**依存関係:** Phase 12 のヘテロ環検出基盤（共通のグラフ解析が必要）

---

### Phase 16 — 4 環以上の縮合多環（優先度: 低）

**対象規則:** IUPAC P-31.1.3, retained names for polycyclic hydrocarbons

**実装する保留名:**

| 名称 | 環数 | 炭素数 |
|------|------|--------|
| pyrene | 4 | 16 |
| chrysene | 4 | 18 |
| triphenylene | 4 | 18 |
| benzo[a]pyrene | 5 | 20 |
| coronene | 7 | 24 |
| perylene | 5 | 20 |

**アーキテクチャ:** 三環検出（Phase 6）の延長として、4環以上の縮合 aromatic cluster を SSSR 数と全炭素数でフィンガープリント照合する方式を採用。fusion nomenclature（furo[2,3-b]pyridine 形式）は別途 Phase 17 として切り出す。

**工数見積:** M（5–7日）  
**依存関係:** Phase 6（三環縮合）

---

## 実装優先度サマリー

```
高優先度（ユーザーへのインパクト大・実装コスト小）
  Phase 7  — 含窒素官能基 (amines, nitriles, nitro)      S: 2–3日
  Phase 8  — チオール・エーテル                          S: 2–3日
  Phase 14 — シクロアルケン                              S: 2–3日

中優先度（IUPAC 網羅性・中程度コスト）
  Phase 9  — 複数官能基・二酸・ポリオール・bis/tris       M: 4–5日
  Phase 10 — エステル・酸ハライド                        M: 4–5日
  Phase 11 — アミド・イミン                              S: 2–3日
  Phase 12 — ヘテロ環 保留名                             L: 7–10日
  Phase 13 — Hantzsch-Widman ヘテロ環                   M: 4–5日

低優先度（特殊構造・高実装コスト）
  Phase 15 — 架橋多環・スピロ                           L: 8–12日
  Phase 16 — 4 環以上の縮合多環                         M: 5–7日
```

**推奨実装順序:**

```
Phase 7 → Phase 8 → Phase 14 → Phase 9 → Phase 11
    → Phase 10 → Phase 12 → Phase 13 → Phase 16 → Phase 15
```

Phase 7, 8, 14 は独立しており並行実装可能。Phase 12 は基盤変更を伴うため Phase 7–11 の後が安全。

---

## 既存実装への改善事項

実装と並行して対処すべき known issues:

| # | 問題 | 影響 | 修正難度 |
|---|------|------|---------|
| 1 | 環外 `–CHO` の命名が `cycloalkanal` → 正しくは `cycloalkanecarbaldehyde` | `ring_handler.py` の `assemble_ring_name()` | S |
| 2 | bis/tris が必要な複合置換基に `di/tri` を使用 | `name_assembler._build_prefix()` | S |
| 3 | ベンゼン系の `–COOH` が環外付加の場合に `benzoic acid` の保留名を使用するが、置換体（例: 4-chlorobenzoic acid）の正確なテストが不足 | `ring_handler.assemble_ring_name()` | S |
| 4 | R/S 立体化学の cyclic compound への適用が不完全 | `stereochemistry.py` | M |

---

## Decision

以下の順で実装を進める:

1. **即時着手**: Phase 7（アミン・ニトリル・ニトロ）→ Phase 8（チオール・エーテル）→ 既存改善 #1, #2
2. **第2波**: Phase 14（シクロアルケン）→ Phase 9（複数官能基）→ Phase 11（アミド）→ Phase 10（エステル）
3. **第3波**: Phase 12（ヘテロ環保留名）→ Phase 13（Hantzsch-Widman）
4. **第4波**: Phase 16（4 環以上多環）→ Phase 15（架橋・スピロ）

各 Phase の開始前にテストファイル `tests/test_phaseN.py` を先行して作成し、TDD 方式で実装する。

---

## Consequences

- **易しくなること**: 有機化合物の広いクラスを系統名で表現できる。特に Phase 7–9 で医薬品中間体・生体分子に多い構造を網羅できる。
- **難しくなること**: ヘテロ原子の導入により `functional_group.py` の seniority 判定が複雑になる。複数官能基の並立処理（principal characteristic group の選択）は現在の単一官能基前提のアーキテクチャに修正を要する。
- **要再検討**: Phase 12 のヘテロ環対応時に `find_rings()` の炭素専用フィルタを緩和するが、これが既存の芳香族命名ロジックに予期しない副作用を与えないよう回帰テストを徹底する。

---

## Action Items

- [ ] Phase 7: `functional_group.py` に amine / nitrile / nitro 検出を追加
- [ ] Phase 7: `tests/test_phase7.py` を TDD で作成（20 ケース以上）
- [ ] Phase 8: thiol / ether suffix・prefix を実装
- [ ] 改善 #1: `ring_handler.assemble_ring_name()` の carbaldehyde 形式修正
- [ ] 改善 #2: `name_assembler._build_prefix()` に bis/tris 判定を追加
- [ ] Phase 14: `ring_handler._assign_ring_locants()` にシクロアルケン対応を追加
- [ ] Phase 9: `detect_groups()` を複数同一官能基に対応させる
- [ ] Phase 12: `heterocycle_handler.py` を新規作成、保留名 lookup table を実装
- [ ] SETUP.md / README を Phase 7 完了後に更新
