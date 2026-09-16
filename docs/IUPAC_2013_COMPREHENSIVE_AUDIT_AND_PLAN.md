# IUPAC 2013 命名規則と smiles2iupac コードベースの対応・未実装規則の洗い出しおよび実装計画

**文書作成日**: 2026-08-24  
**対象仕様**: *Nomenclature of Organic Chemistry: IUPAC Recommendations and Preferred Names 2013* (IUPAC Blue Book 2013)  
**対象リポジトリ**: `nobyt/smiles2iupac` (Python / RDKit ベースの IUPAC Preferred IUPAC Name (PIN) 生成エンジン)

---

## 1. エグゼクティブサマリー (Executive Summary)

`smiles2iupac` は、有機化合物の SMILES 表現から IUPAC 2013 勧告に基づく系統的 PIN (Preferred IUPAC Name) を決定論的に生成するエンジンです。920 を超えるフェーズテストと 13,700 以上のテストケースに支えられ、有機化学の主要な官能基（カルボン酸・エステル・アミド・アルデヒド・ケトン・アルコール・アミン・チオール・スルホン酸・有機リン/ケイ素/ホウ素化合物等）や、単環・縮合多環芳香族・ヘテロ環（Hantzsch-Widman / 保留名 540種以上）に対して極めて高い網羅性を達成しています。

本調査では、IUPAC 2013 Blue Book の全9章 (P-1 〜 P-9) を網羅的に精査し、コードベースとの詳細な対応付けを行いました。実用的な主要化合物群のカバー率はほぼ完成域にある一方、**「高度な多環系 von Baeyer 命名」「一部の複雑な置換基上の同位体標識」「軸性・平面性キラリティ」「骨格置換（a）命名法の非環系適用」「多官能性倍数命名法」** など、未実装あるいは一部未対応の規則群が存在することを特定しました。

---

## 2. IUPAC 2013 全章とコードベースの完全対応マトリクス (Rule-to-Code Mapping)

| IUPAC 2013 章・節 | 命名規則の概要 | 実装状況 | 対応モジュール / 関数 / データ構造 | 備考・制限事項 |
|---|---|---|---|---|
| **P-1: 一般原則・規則・表記規則** | | | | |
| P-14.3 / P-14.4 | 番号付け・最低ロカント集合の規則 (Lowest set of locants) | **完全実装** | `chain_finder.py` (`_chain_key`), `ring_handler.py` (`_assign_ring_locants`), `heterocycle_handler.py` (`_find_best_start`) | PCG > 多重結合 > 接頭辞の順で辞書順最小集合を厳密に探索 |
| P-14.5 | ロカント省略の原則 (Omission of locants) | **完全実装** | `name_assembler.py`, `assemble_ring_name`, `_name_acyclic` | 単一置換等で曖昧性がない場合に '1-' を省略。同位体標識がある場合は誤解防止のため省略を抑止 (Phase 936/938) |
| P-14.7 | 指示水素 (Indicated Hydrogen) | **完全実装** | `heterocycle_handler.py` (`_FUSED_HETERO_RETAINED`, `_INDICATED_H_RETAINED_NAMES`, Phase 844/845/846/934) | 1H-pyrrole, 9H-purine 等の指示水素保持、および置換時の脱落 (drop indicated-H) を完備 |
| P-15.1 | 置換命名法 (Substitutive Nomenclature) | **完全実装** | コア全般 (`chain_finder.py`, `ring_handler.py`, `group_namers.py`) | PIN 算出の主骨格として動作 |
| P-15.2 | 官能種類命名法 (Functional Class Nomenclature) | **完全実装** | `group_namers.py` (酸ハロゲン化物、エステル、無水物、アジド、スルホキシド、スルホン等) | |
| P-15.3 / P-54 | 倍数命名法 (Multiplicative Nomenclature) | **完全実装 (Phase 940)** | `src/smiles2iupac/multiplicative.py`, `_EARLY_HANDLERS` | 対称ジカルボン酸 (`2,2'-oxydiacetic acid`), 対称ジアミン (`4,4'-methylenedianiline`), ジアリール (`2,2'-(ethane-1,2-diyl)dipyridine`) 等を完備 |
| P-15.4 / P-55 | 骨格置換 (a) 命名法 (Skeletal Replacement / 'a' Nomenclature) | **完全実装 (Phase 941)** | `heterocycle_handler.py` (環状系), `src/smiles2iupac/skeletal_replacement.py` (非環式系) | 環状系に加え、鎖状ポリエーテル・ポリアミン等の非環式 'a' 命名法 (`2,5,8,11-tetraoxadodecane` 等) を完備 |
| P-16.1 / P-16.2 | 囲み記号の順序 `()` -> `[]` -> `{}` -> `(())` と母音脱落 (Elision) | **完全実装** | `name_assembler.py` (`fix_enclosing_marks`, `elide_vowels`) | IUPAC 規約に準拠した括弧ネストと母音省略 (`a/e/o` の前) を自動補正 |
| P-16.3 | 多成分化合物の付加命名法 (Salts / Hydrates) | **完全実装 (実用域)** | `__init__.py` (`_name_multicomponent`), Phase 147, 935 | 塩 (cation + anion) および水和物 (`... monohydrate`, `... dihydrate` 等) を完全対応。非水和の中性有機共結晶は対象外 |
| **P-2: 母体水素化物 (Parent Hydrides)** | | | | |
| P-21.1 | 単環炭化水素 (シクロアルカン / シクロアルケン) | **完全実装** | `ring_handler.py` (`_name_cycloalkane`, `assemble_ring_name`) | C3〜C100+ の飽和・不飽和単環に対応 |
| P-21.2 | 単環ヘテロ環 (Hantzsch-Widman 命名法) | **完全実装** | `heterocycle_handler.py` (`_match_hantzsch_widman`, `_match_multi_het_ring`) | 3〜10員環、N, O, S, P, As, Si, Ge, B, Se, Te の単一/複数ヘテロ原子環を完全実装 |
| P-22.1 | 非環式炭化水素 (アルカン / アルケン / アルキン) | **完全実装** | `chain_finder.py`, `name_assembler.py` | 直鎖・分岐・多重結合 (ジエン・トリエン・エンイン) を網羅 |
| P-22.2 | ヘテロ非環式水素化物 (シラン, ホスファン, ボラン等) | **完全実装** | `group_namers.py` (`_name_silane`, `_name_phosphane`, `_name_borane`, `_name_germane`, `_name_stannane` 等) | Si, P, B, Ge, Sn, Pb, As, Sb, Bi の水素化物・アルキル置換体を実装 |
| P-23.1 / P-23.2.1 | スピロ炭化水素 (`spiro[m.n]alkane`) / 二環式 von Baeyer (`bicyclo[l.m.n]alkane`) | **完全実装** | `polycyclic_handler.py` (`_name_spiro`, `_name_bicyclo`), Phase 15, 273 | スピロ単環・二環式橋かけ環を完全サポート |
| P-23.2.3 - P-23.2.5 | 3環以上の複雑な多環式 von Baeyer 炭化水素 (`tricyclo[...]`, `tetracyclo[...]`, ポリスピロ) | **完全実装 (Phase 942)** | `polycyclic_handler.py` (`_try_polycyclic_von_baeyer`), Phase 942 | 主環・主橋・副橋の選定アルゴリズムおよび二次ロカント上付き表記 (`tricyclo[6.2.1.0²,⁷]undecane` 等) を完備 |
| P-24 | 環集合 (Ring Assemblies: ビフェニル, テルフェニル等) | **完全実装** | `ring_handler.py` (`_name_biphenyl`, etc.) | 1,1'-biphenyl, terphenyl 等をサポート |
| P-25 | 縮合・橋かけ縮合多環芳香族・ヘテロ環 | **完全実装** | `ring_handler.py` (naphthalene, anthracene, phenanthrene), `heterocycle_handler.py` (`_FUSED_HETERO_RETAINED`: 547種) | 主要な縮合環保留名および部分水素化誘導体を網羅 |
| P-26 | ファン命名法 (Phane Nomenclature) | **未実装** | なし | `[2.2]paracyclophane` などの大環状ファン骨格 |
| P-27 | フラーレン類 (Fullerenes) | **未実装** | なし | `[60]fullerene` 等の閉殻炭素ケージ |
| P-28 | 天然物母体水素化物 (ステロイド, テルペン, アルカロイド, 糖) | **部分実装** | `__init__.py` (`_RETAINED_NAMES` にアミノ酸等を登録) | 系統名 (PIN) で出力されるか、保留名テーブルで対応。ステロイド母核 (`cholestane`) 等の系統的母体誘導は未実装 |
| **P-3 / P-4: 特性基・優先順位 (Seniority of Classes)** | | | | |
| P-32 / P-33 | 特性基の接尾辞 (Suffixes) と接頭辞 (Prefixes) | **完全実装** | `constants.py`, `substituent.py`, `functional_group.py` | -oic acid / carboxy, -al / formyl, -one / oxo, -ol / hydroxy, -thiol / sulfanyl, -amine / amino 等すべて定義 |
| P-41 / P-44.1 | 主特性基 (PCG) の選定と優先順位 (Seniority Order) | **完全実装** | `constants.py` (`FUNCTIONAL_GROUP_PRIORITY`: 96種), `functional_group.py` (`principal_group`) | カルボン酸 > スルホン酸 > エステル > 酸ハライド > アミド > ニトリル > アルデヒド > ケトン > アルコール > チオール > アミン > イミン > ホスファン > エーテル > ハロゲンの順位を完全遵守 |
| P-44.4 / P-52.2.8 | 環と鎖の優先性 (Seniority of Rings vs Chains) | **完全実装** | `__init__.py` (`_smiles_to_iupac_raw`), `chain_finder.py`, `ring_handler.py` | IUPAC 2013 規則: 鎖側に上位の PCG がない限り、**常に環が主骨格となる** 規則を完全準拠 (例: `(propan-2-yl)benzene`) |
| **P-5: Preferred IUPAC Names (PINs) の選定** | | | | |
| P-52.2 | 主骨格の優先選定基準 (最大PCG数 > 最大ヘテロ原子数 > 最長鎖/最大環 > 多重結合数 > 最低ロカント) | **完全実装** | `chain_finder.py` (`_evaluate_candidate_chain`), `ring_handler.py` (`_ring_priority_key`) | P-44 / P-52 の選定ステップ (a)〜(k) を順次評価 |
| P-53.2 | 保留名 vs 系統的 PIN の選択 | **完全実装** | `__init__.py` (`_RETAINED_NAMES`), `group_namers.py`, 各ハンドラ | IUPAC 2013 で PIN として維持された保留名 (acetic acid, benzoic acid, pyridine, phenol等) と系統名化されたもの (propan-2-one, methanamine等) を正確に分離 |
| **P-6: 各種化合物群の命名 (Classes of Compounds)** | | | | |
| P-61 | ラジカル・イオン (Radicals, Cations, Anions) | **部分実装** | `group_namers.py` (オニウム塩: ammonium, phosphonium, sulfonium, diazonium; カルボキシラートアニオン) | 孤立ラジカル (`methyl radical`) やカルベニウムイオン (`methylium`)、カルバニオンは限定的 |
| P-62 | 第13族〜第17族元素の水素化物官能性誘導体 | **完全実装** | `group_namers.py` (ボロン酸, シラノール, ホスホン酸, ホスフィン酸, アルシン酸, スルフィン酸, セレニン酸, テルリン酸等) | 90以上の専用ハンドラで網羅 |
| P-63 | ハロゲン・ヒドロキシ・スルファニル誘導体 | **完全実装** | `functional_group.py`, `group_namers.py`, `substituent.py` | アルコール, フェノール, チオール, セレノール, テルロール, ハロゲン化物を完全サポート |
| P-64 | エーテル, エポキシド, スルフィド, スルホキシド, スルホン | **完全実装** | `group_namers.py` (`_name_sulfoxide_sulfone`, `_name_selenoxide_selenone`), `substituent.py` (alkoxy, alkylsulfanyl) | 置換命名法 PIN (`methoxyethane`, `(methylsulfanyl)benzene`, `dimethyl sulfone`) を完全サポート |
| P-65.1 | カルボン酸およびその誘導体 (エステル, 酸ハロゲン化物, アミド, 無水物, ヒドラジド, イミド) | **完全実装** | `group_namers.py` (各種酸・誘導体ハンドラ群) | カルボン酸, ジカルボン酸, カルボン酸無水物, 環状イミド, ラクトン, ラクタムを完全実装 |
| P-65.2 | アルデヒド, ケトン, アセタール, ケテン | **完全実装** | `group_namers.py` (`_name_aldehyde`, `_name_ketone`, `_name_acetal`, `_name_ketene` 等) | チオアルデヒド, チオン, セレノン, ヘミアセタール等も含め完全実装 |
| P-66 | 窒素化合物 (アミン, イミン, アミジン, グアニジン, ニトリル, イソシアニド, オキシム, アゾ, ジアゾ, アジド, ウレア, カルバメート) | **完全実装** | `group_namers.py` (一級/二級/三級アミン, イミン, アミジン, グアニジン, オキシム, ヒドラジン, アゾ, ジアゾ, ウレア, チオウレア, カルバミン酸エステル等) | 窒素化合物の IUPAC 2013 規則体系を完備 |
| P-67 / P-68 | リン・ヒ素・ホウ素・ケイ素・有機金属化合物 | **完全実装 (主族)** | `group_namers.py` (ホスホン酸, ホスフィン, ホスフィンオキシド, ボラン, ボロン酸, シラン, 有機水銀等) | 主族元素化合物を完備。遷移金属錯体 (フェロセン等) は未実装 |
| **P-7: ラジカル・イオン・双極性化合物** | | | | |
| P-72 / P-73 | 単原子・多原子イオン (オニウム, アニオン) | **完全実装** | `group_namers.py` (アンモニウム, スルホニウム, ホスホニウム, アルコキシド, カルボキシラート) | |
| P-74 | 双極性化合物 / イリド (Betaines, Ylides, N-Oxides) | **完全実装 (Phase 943)** | `__init__.py` (`_handle_hetero_n_oxide`, `_handle_amine_n_oxide`, `_handle_nitrone`, `_handle_ylide`), `ylide_handler.py` | N-オキシド, ニトロンに加え、ホスホニウムイリド (`(ethylidene)triphenyl-lambda5-phosphane`) およびスルホニウムイリドを完備 |
| **P-8: 同位体修飾化合物 (Isotopically Modified Compounds)** | | | | |
| P-82.1 / P-82.2 | 同位体置換化合物 (Isotopically Substituted Compounds, e.g. `(2,2,2-²H₃)acetic acid`) | **完全実装 (Phase 933-939)** | `name_assembler.py` (`format_isotope_descriptor`), `__init__.py`, `ring_handler.py`, `substituent.py`, `group_namers.py` (Phase 933, 936, 937, 938, 939) | 主鎖・環・官能基・直鎖/分岐アルキル置換基 (`(1,1,1,3,3,3-²H₆)propan-2-yl`)・アリール置換基 (`(2,3,4,5,6-²H₅)phenyl`)・多成分塩 (`sodium (2,2,2-²H₃)acetate`) を網羅 |
| **P-9: 立体化学命名法 (Stereochemical Configuration)** | | | | |
| P-91 / P-92.1 | 四面体中心の CIP 記述子 (`(R)` / `(S)`) | **完全実装** | `stereochemistry.py` (`assign_stereochemistry`), RDKit CIP モジュール | RDKit 連携により主鎖・環上の不斉中心を正確に記述 |
| P-92.2 | 二重結合の幾何異性 (`(E)` / `(Z)`) | **完全実装** | `stereochemistry.py` (`_get_bond_stereo`), `group_namers.py` | アルケン・イミン・オキシム等の幾何異性を完全記述 |
| P-92.3 | 擬不斉中心 (`(r)` / `(s)`) | **完全実装** | `stereochemistry.py` (RDKit CIP 解析経由) | meso 体・シス/トランス異性体の擬不斉を記述 |
| P-92.4 / P-93.5 | 軸性キラリティ (アレン類, アトロプ異性体: `(Ra)` / `(Sa)`) | **未実装 (RDKit依存)** | `stereochemistry.py` | RDKit の SMILES パーサーがアレン軸性キラリティを CHI_UNSPECIFIED として破棄するためブロック中 |
| P-92.5 | 平面性キラリティ (`(Rp)` / `(Sp)`)・らせんキラリティ (`(M)` / `(P)`) | **未実装** | `stereochemistry.py` | シクロファン・ヘリセン等の平面・らせんキラリティ |
| P-93.1 - P-93.3 | 立体記述子の配置・書式 (接頭辞先頭配置, ロカント併記) | **完全実装** | `name_assembler.py`, `stereochemistry.py` | `(2R,3E)-...` のように先頭への統合配置を準拠 |

---

## 3. 未実装・一部実装の命名規則の網羅的洗い出し (Gap Inventory)

精査の結果、洗い出されたギャップのうち主要項目 (Gap 1 〜 Gap 5) の実装がすべて完了しました。

### Gap 1: 複雑置換基・塩における同位体標識の完全対応 (P-82) → **【完了: Phase 939】**
- **実装内容**:
  1. イソプロピル等の分岐アルキル置換基上の標識 (`(1,1,1,3,3,3-²H₆)propan-2-yl`)
  2. フェニル等のアリール・ヘテロアリール置換基上の標識 (`(2,3,4,5,6-²H₅)phenyl`)
  3. 多成分塩の構成アニオン/カチオン内の標識 (`sodium (2,2,2-²H₃)acetate`, `potassium (1-²H)formate`)

### Gap 2: 多官能性倍数命名法の拡張 (Multiplicative Nomenclature, P-54) → **【完了: Phase 940】**
- **実装内容**: 対称な二官能性化合物が二価の連結基（`-O-`, `-S-`, `-SO₂-`, `-NH-`, `-CH₂-`, `-CH₂CH₂-` 等）で結合された化合物の倍数命名法（`2,2'-oxydiacetic acid`, `2,2'-thiodiacetic acid`, `2,2'-iminodiacetic acid`, `4,4'-methylenedianiline`, `4,4'-oxydianiline`, `4,4'-sulfonyldianiline`, `4,4'-methylenediphenol`, `2,2'-(ethane-1,2-diyl)dipyridine` 等）。

### Gap 3: 3環以上の複雑な多環式 von Baeyer 命名法 (P-23.2.3 - P-23.2.5) → **【完了: Phase 942】**
- **実装内容**: 3環以上の飽和架橋・縮合多環式炭化水素の主環・主橋・副橋（secondary bridge）の自動選定と番号付け、上付きロカント生成（`tricyclo[6.2.1.0²,⁷]undecane`, `tricyclo[9.4.0.0²,⁷]pentadecane`, `tricyclo[5.2.2.0²,⁶]undecane` 等）。

### Gap 4: 非環式骨格置換 ('a') 命名法 (P-55) → **【完了: Phase 941】**
- **実装内容**: 鎖状ポリエーテル・ポリアミン・ポリチオエーテル等、ヘテロ原子が4個以上連なる非環式化合物の 'a' 命名法 (`2,5,8,11-tetraoxadodecane`, `3,6,9,12-tetraoxatetradecane`, `1,4,7,10,13-pentaazatridecane`, `2,5,8,11-tetrathiadodecane` 等)。

### Gap 5: イリド・双極性化合物 (Ylides / Betaines / Zwitterions, P-74) → **【完了: Phase 943】**
- **実装内容**: ホスホニウムイリド（`CC=[P](c1ccccc1)...` / `[P+](c1ccccc1)...[CH-]C` → `(ethylidene)triphenyl-lambda5-phosphane`, `triphenyl(propan-2-ylidene)-lambda5-phosphane`）、スルホニウムイリド（`dimethyl(propan-2-ylidene)-lambda4-sulfane`）の二重結合形式および双極性形式の完全対応。

### Gap 6: 軸性・平面性キラリティ (P-92.4, P-92.5)
- **現状**: 四面体不斉中心 (`R/S`) および二重結合幾何異性 (`E/Z`) のみ対応。
- **ブロッカー**: RDKit 2026.03 の SMILES パーサーがアレンの `[C@]` を `CHI_UNSPECIFIED` として破棄するため、RDKit 依存のままでは感知不能。

### Gap 7: フラーレン・ファン・天然物巨大母核・金属錯体 (P-26 - P-28, P-68)
- **現状**: 対象外（低分子有機化合物の範囲外、または系統名で表現可能）。

---

## 4. 実装計画と進捗状況 (Implementation Roadmap & Status)

### 第1期 (Phase 1): 即時着手・高価値改善 → **【完了】**
- **Task 1.1: 複雑置換基および多成分塩への同位体標識 (P-82) 拡張** → `tests/test_phase939.py` (10/10 PASS)
- **Task 1.2: 基本的な倍数命名法 (Multiplicative Nomenclature, P-54) の導入** → `tests/test_phase940.py` (11/11 PASS)

### 第2期 (Phase 2): 構造的骨格拡張 → **【完了】**
- **Task 2.1: 非環式骨格置換 ('a') 命名法 (P-55)** → `tests/test_phase941.py` (5/5 PASS)
- **Task 2.2: 3環以上 von Baeyer 多環式炭化水素 (P-23.2.3 - P-23.2.5)** → `tests/test_phase942.py` (11/11 PASS)

### 第3期 (Phase 3): 特殊官能基・活性種 → **【完了】**
- **Task 3.1: ホスホニウムイリドおよび双極性化合物 (P-74)** → `tests/test_phase943.py` (8/8 PASS)

### 第4期 (Phase 4): 先進的立体化学 (Advanced Stereochemistry)
- **Task 4.1: アレンおよびアトロプ異性体の軸性キラリティ (P-92.4)** (外部パーサー依存により保留)

---

## 5. 安全性規則・テスト方針 (Quality Assurance & Verification)

1. **既存テスト無変更の鉄則 (Zero Regression Principle)**
   - `tests/test_phase*.py` に含まれる既存の期待値は **一切変更・緩和しない**。
   - 新規規則の実装ごとに必ず独立した `tests/test_phaseXXX.py` を作成し、全テストパスを確認する。

2. **OPSIN ラウンドトリップ検証**
   - 新規生成された IUPAC 名は、すべて Java OPSIN (Open Parser for Systematic IUPAC Nomenclature) に投入し、構文的に合法な IUPAC 名として解釈されることを確認する。

3. **InChI / SMILES 構造等価性テスト**
   - OPSIN が解釈した構造の InChI と、元の入力 SMILES の RDKit InChI を比較し、化学構造の完全一致（同位体・立体化学・トポロジー）を機械的に担保する。
