# 開発計画: 未知アンプリコン向け Nanopore リード処理ツール

> ブランチ: `refless`  
> 作成日: 2026-07-03  
> 関連: `tool_concept`（初期アイデア）、現行 `miaoseq`（`R/functions.R`）

---

## 1. ツールの位置づけ

### 1.1 何をするツールか

本ツールは **16S 特化ツールではない**。参照ゲノム配列や既知のアンプリコン配列を前提としたリード処理ツールでもない。

**前提:**

- **既知:** インデックス配列、PCR プライマー配列
- **未知:** プライマーより内側に増幅されるアンプリコン配列（増幅対象領域）

**目的:**

インデックス付き Nanopore PCR リードから、増幅対象領域の真の配列を以下の流れで推定する。

```text
FASTQ
  → Demultiplex（サンプル割当）
  → Read clustering（サンプル内の配列多様性の分離）
  → Consensus generation（クラスタごとの代表配列）
  → Confidence evaluation（信頼度評価）
  → 推定アンプリコン配列（クラスタごと）
```

**下流解析はユーザー次第:**

- Taxonomy 解析（推定コンセンサスを分類 DB に照合）
- ゲノム編集パターン解析（推定コンセンサスを参照ゲノムと比較）

コアパイプラインは参照配列に依存せず動作し、下流モジュールがオプションとして接続される設計とする。

### 1.2 現行 miaoseq との関係

| 領域 | 方針 |
|------|------|
| Demultiplex | **大幅改修**（現行 BLAST + index のみは置き換え） |
| Clustering / Consensus / 信頼度 | **新規構築** |
| 参照ゲノムベースの編集パターン同定（下記モジュール全体） | **コアから分離し、統合・最適化して残す** |
| 複数コンセンサスからの編集代表選択 | **新規構築**（上記モジュール内） |

現行アルゴリズムは大幅に削除・改修してよい。ただし、参照ゲノム配列と比較してゲノム編集パターンを同定する機能は、下記 **Edit-calling Module** として統合的に維持・最適化する。

#### 現行コードの整理: 別関数だが同一目的

現行 miaoseq では以下が **別関数** として実装されているが、いずれも「リードを参照ゲノム上の編集標的に位置づけ、編集パターンを検出する」という **同一の目的** を担う。

| 現行関数 | 役割 |
|---------|------|
| `prepAmpliconDB()` | 参照ゲノム + プライマーから amplicon 座標・配列を事前構築 |
| `doAlign()` | 各リードを amplicon DB に BLAST アラインメント |
| `doEditcall()` | アラインメント結果から編集型を集計・判定 |

これらは demultiplex 以前の段階（`prepAmpliconDB`）と demultiplex 以降（`doAlign` → `doEditcall`）に分かれているが、論理的には **1 つの Edit-calling Module** である。新ツールではこれらを統合再設計し、コアパイプライン（demux → cluster → consensus）から分離する。

#### 現行 vs 新ツールの根本的な違い

```text
【現行 miaoseq】
prepAmpliconDB（参照ゲノム + プライマー → amplicon DB）
    ↓
demultiplex
    ↓
doAlign: 全リード × amplicon を BLAST（遅い）
    ↓
doEditcall: リード単位の編集型集計

【新ツール】
Core: demux → cluster → consensus（参照ゲノム不要）
    ↓
Edit-calling Module（オプション）:
  prepRefTargets（参照ゲノム + プライマー + gRNA → 標的座標）
    ↓
  alignConsensusToRef: コンセンサス配列のみを参照ゲノムにアラインメント（高速）
    ↓
  callEditsFromConsensus: コンセンサス単位の編集パターン判定
    ↓
  selectEditCluster: 複数コンセンサスから代表を選択
```

**要点:** 生リードではなく、demux + clustering で得た **コンセンサス配列（サンプルあたり複数あり得る）** を参照ゲノムにアラインメントする。BLAST による全リードアラインメントは廃止する。

---

## 2. 入力仕様

### 2.1 必須入力（コアパイプライン）

| 入力 | 説明 |
|------|------|
| FASTQ | Basecalled リード（必須） |
| インデックス情報（CSV 等） | サンプル名、Forward index、Reverse index |
| PCR プライマー配列 | Forward primer、Reverse primer |

### 2.2 オプション入力

| 入力 | 用途 |
|------|------|
| POD5 | Basecalling（Dorado 等） |
| BAM | Basecalled リードの代替入力 |

### 2.3 下流解析用入力（モジュール別）

| モジュール | 追加入力 |
|-----------|---------|
| Taxonomy | 分類 DB（SILVA, NCBI nt 等）、分類器設定 |
| ゲノム編集解析 | **参照ゲノム配列（必須）**、**ガイド RNA 配列情報（必須）**、PAM 情報（必要に応じて） |

ゲノム編集パターン同定では、ガイド RNA 配列情報を **必須入力** とする。複数コンセンサスが得られた場合にどれを編集判定の対象とするかは、ガイド RNA 周辺の編集シグナル等を用いた新規アルゴリズムで決定する。

---

## 3. 全体アーキテクチャ

```text
                    ┌─────────────────────────────────────┐
                    │         Core Pipeline               │
                    │  (参照配列・分類 DB 不要)            │
                    └─────────────────────────────────────┘
FASTQ ──► ① Demultiplex ──► ② Clustering ──► ③ Consensus ──► ④ Confidence
                                      │
                                      ▼
                         推定アンプリコン（クラスタごと）
                                      │
                    ┌─────────────────┴─────────────────┐
                    ▼                                   ▼
           ┌────────────────┐              ┌──────────────────────────────┐
           │ Taxonomy Module │              │ Edit-calling Module          │
           │ (オプション)     │              │ (オプション・統合再設計)      │
           │                 │              │ ・prepRefTargets             │
           │                 │              │ ・alignConsensusToRef        │
           │                 │              │ ・callEditsFromConsensus     │
           │                 │              │ ・selectEditCluster          │
           │                 │              │ 入力: 参照ゲノム + gRNA 必須  │
           └────────────────┘              └──────────────────────────────┘
```

---

## 4. コアパイプライン詳細

### 4.1 ① Demultiplex

#### 設計思想

検索対象は **Index のみ** ではなく **Index + Primer** の複合配列とする。

```text
Forward end:  IndexF + PrimerF
Reverse end:  PrimerR + IndexR  （reverse complement を考慮）
```

両端を独立に探索し、以下を区別する。

- Forward のみ一致
- Reverse のみ一致
- 両方一致

#### スコアリング

各端について算出:

- alignment score
- identity
- mismatch 数
- gap 数

さらに **best hit と second-best hit のスコア差** を利用する。

```text
Sample1: score 58
Sample2: score 42  → 十分区別可能 → High confidence

Sample1: score 58
Sample2: score 57  → 差が小さい → Ambiguous
```

#### 距離関数: edit distance（Levenshtein）

Nanopore リードは indel エラーが多いため、Hamming 距離ではなく **edit distance** を推奨する。

- mismatch / insertion / deletion を unit cost でカウント
- 正規化: `identity = 1 - edit_distance / query_length`
- 閾値はエラー率ベースで設定（例: 最大 15–20% エラー許容）

**補足:** Barbell（bioRxiv 2025）では、フランク位置特定に edit distance、バーコード識別に subsequence-based scoring を併用している。単純な edit distance だけでなく、**位置制約付き部分一致スコア** の併用も検討する。

#### 出力分類

各 read を以下の 4 段階に分類する。

| 分類 | 条件（案） |
|------|-----------|
| **High confidence** | 両端一致、スコア ≥ T_high、best − second ≥ δ |
| **Low confidence** | 片端一致、またはスコア ≥ T_low |
| **Ambiguous** | best − second < δ |
| **Unclassified** | スコア < T_low |

#### 現行実装との差分

現行 `doDemultiplex()` の問題点:

| 問題 | 詳細 |
|------|------|
| Index のみ検索 | Primer を demultiplex に未使用 |
| `complete_match` が厳しすぎる | mismatch=0 & gap=0 を要求 → ONT では未分類が増加 |
| 信頼度分類なし | High / Low / Ambiguous / Unclassified がない |
| best vs second-best 判定なし | スコア差による曖昧判定がない |
| BLAST 依存 | read 末端の短い local alignment には過剰 |

**改修方針:** `doDemultiplex2()` を新規実装し、末端ウィンドウに対する semi-global alignment（`pwalign` または Edlib）に置き換える。

---

### 4.2 ② Read clustering

各サンプルごとに、単一コンセンサスを作るのではなく **まずクラスタリング** する。

```text
sample1
  cluster1  92%   ← 主アンプリコン
  cluster2   7%   ← コンタミ候補
  cluster3   1%   ← chimera / 低頻度変異体
```

**目的:**

- コンタミネーションの検出
- PCR chimera の分離
- 複数テンプレート・複数編集パターンの共存検出

#### 方針: 3 パイプライン切替可能

クラスタリング手法は **ユーザーが切り替え可能** とする。いずれも共通の後段処理（リファインメント → コンセンサス）に接続する。

| パイプライン ID | 手法 | 位置づけ |
|--------------|------|---------|
| `"meshclust"` | meshclust のみ | **本番推奨**。ONT アンプリコンでの精度実績（LACA） |
| `"mmseqs2"` | MMseqs2 のみ | **高速プロトタイプ代替**。開発初期・大規模データの速度検証 |
| `"umap_meshclust"` | UMAPclust → meshclust | **高精度候補**。LACA でノイジー ONT に最良の組み合わせ |

```text
                    ┌─────────────────────────────────────────┐
                    │  cluster_method パラメータで切替         │
                    └─────────────────────────────────────────┘
         meshclust only │  mmseqs2 only  │  UMAPclust → meshclust
                        ▼                ▼                ▼
                    ┌─────────────────────────────────────────┐
                    │  Stage 2: minimap2 リファインメント      │
                    │  （全パイプライン共通・オプション）       │
                    └─────────────────────────────────────────┘
                                        ▼
                    ┌─────────────────────────────────────────┐
                    │  Stage 3: pileup consensus              │
                    │  （+ オプション Racon/Medaka）            │
                    └─────────────────────────────────────────┘
```

**API 案:**

```r
doCluster(
  reads_fasta,
  cluster_method = c("meshclust", "mmseqs2", "umap_meshclust"),
  refine = TRUE,          # Stage 2 リファインメント
  meshclust_id = 0.90,    # meshclust --id
  mmseqs_min_seq_id = 0.85,
  umap_min_samples = 50,  # UMAPclust パラメータ
  ...
)
```

---

#### パイプライン 1: meshclust のみ（本番推奨）

**根拠:**

- [LACA](https://doi.org/10.1080/19490976.2025.2516703) ベンチマークで、ノイジー ONT リードのクラスタリングに isONclust より優位
- k-mer ベース（alignment-free）+ オプション alignment で identity 閾値をデータから推定可能
- 未知アンプリコン内の **異なるテンプレート分離** に適する（isONclust は同一遺伝子ファミリーの統合向けで本用途には不向き）

**コマンド例:**

```bash
meshclust reads.fasta --id 0.90 --align --output clusters.clstr
```

**注意点:**

- `--id` の設定が精度に直結（厳しすぎると過分割、緩すぎると異種混在）
- ストランド情報（正負鎖）で過分割する場合あり → Stage 2 リファインメントで統合

---

#### パイプライン 2: MMseqs2 のみ（高速プロトタイプ代替）

**根拠:**

- 大規模データで meshclust より **数倍〜数十倍高速**
- 実装・デバッグが容易（CLI が安定、依存が少ない）
- パイプライン全体のプロトタイプ構築・スループット検証に最適

**コマンド例:**

```bash
mmseqs easy-linclust reads.fasta clusterRes tmp \
  --min-seq-id 0.85 \
  -c 0.8 \
  --cov-mode 0 \
  --alignment-mode 3
```

**注意点:**

- ONT 専用エラーモデルなし。indel 多発領域で過分割の可能性
- 精度ベンチマークでは meshclust / UMAPclust+meshclust に劣る可能性 — 本番運用前に `"meshclust"` または `"umap_meshclust"` への切替を推奨

**用途:** Phase 2 初期開発、大規模サンプル数の速度検証、パラメータ探索

---

#### パイプライン 3: UMAPclust → meshclust（高精度候補）

**根拠:**

- LACA 論文で **UMAPclust が純度・NMI 最高**、ノイジー ONT では **UMAPclust + meshclust の組み合わせが最良**
- UMAPclust: 5-mer frequency による alignment-free 粗クラスタリング
- meshclust: UMAPclust の出力を統合・identity 閾値で精緻化

```text
reads.fasta
  → UMAPclust（5-mer frequency、粗クラスタ）
  → meshclust（クラスタ統合・identity 精緻化）
  → Stage 2 リファインメント
```

**注意点:**

- UMAPclust は正負鎖で過分割しやすい → meshclust 統合で改善（LACA Fig. 2）
- isONclust + UMAPclust + meshclust の 3 連続は過分割を招く — **UMAPclust → meshclust の 2 段のみ** を採用
- 計算コストの低い順: MMseqs2 < meshclust のみ < UMAPclust → meshclust

**用途:** 精度重視の本番解析、ノイジー ONT（R9.4.1 等）、近縁テンプレートの分離が重要な場合

---

#### 手法比較サマリー

| | meshclust | MMseqs2 | UMAPclust → meshclust |
|--|-----------|---------|----------------------|
| 速度 | ○ | ◎ | △ |
| ONT 精度（文献） | ○ | △ | ◎ |
| 実装容易性 | ○ | ◎ | △ |
| 本番推奨 | **◎** | △（プロトタイプ用） | ○（高精度が必要な場合） |
| 異種テンプレート分離 | ○ | ○ | ◎ |
| 近縁配列分離（編集パターン差） | △ 要 Stage 2 | △ 要 Stage 2 | ○ 要 Stage 2 |

---

#### 共通後段処理（全パイプライン）

ONT indel の問題を考慮し、いずれのパイプラインでも以下を共通適用する。

```text
Stage 2（クラスタ内リファインメント）: minimap2 アラインメント + 距離再計算
  → indel を考慮した真のクラスタ境界を確定
  → ストランド統合、過分割クラスタのマージ
  → chimera 検出（UCHIME / 参照なし chimera フィルタ）

Stage 3: pileup consensus（+ オプション Racon/Medaka）
```

`refine = FALSE` で Stage 2 をスキップ可能（ベンチマーク比較用）。

---

#### Phase 2 ベンチマーク計画

3 パイプラインを模擬データ（既知テンプレート mix + ONT エラー注入）で比較し、デフォルトを決定する。

| 評価指標 | 内容 |
|---------|------|
| クラスタ純度（purity） | 正解ラベルとの一致 |
| NMI | クラスタリング品質の総合指標 |
| コンセンサスエラー率 | Stage 3 後の配列精度 |
| 処理時間 | サンプルあたり・リード数あたり |
| 過分割率 | 真のクラスタ数に対する過剰クラスタ数 |

**検証項目:**

1. 3 パイプライン × Stage 2 有無 の 6 条件比較
2. meshclust `--id` / MMseqs2 `--min-seq-id` の感度解析
3. UMAPclust パラメータ（`min_samples`, `min_bin_size`）の最適化
4. 本番デフォルト: `"meshclust"`、高精度モード: `"umap_meshclust"`、開発用: `"mmseqs2"`

**isONclust は本用途（異種テンプレート分離）に不向きのためベンチマーク対象外。**

**追加検討:** chimera 検出・除去（UCHIME 等）。

---

### 4.3 ③ Consensus generation

各クラスタについて:

```text
Reads → Alignment → Consensus
```

**基本:** minimap2 / MAFFT 等でアラインメント → pileup による majority consensus

**高精度化（オプション）:**

- Racon（iterative polishing）
- Medaka（ONT 向け polish）

**出力:** クラスタごとのコンセンサス配列（プライマー間の増幅領域）

---

### 4.4 ④ 信頼度評価

各コンセンサスについて以下を算出する。

| メトリクス | 例 |
|-----------|-----|
| クラスタサイズ | 1,523 reads |
| クラスタ割合 | 95.3% |
| 平均 identity | 99.6% |
| 各塩基支持率 | position 234: A 98%, G 2% |
| 平均 QV | 可能なら算出 |
| 配列長 | 異常長の検出 |

#### 最終出力例

```text
Sample01
──────────────────
Cluster1
  Support:    1,452 reads
  Fraction:   97.1%
  Consensus:  1,458 bp
  Mean identity: 99.7%

Cluster2
  Support:    31 reads
  Fraction:   2.1%
  Consensus:  1,457 bp
```

Taxonomy や編集パターンは、この出力を下流モジュールに渡して解析する。

---

## 5. 下流モジュール

### 5.1 Taxonomy Module（オプション）

- 入力: クラスタごとのコンセンサス配列
- 手法: BLAST（16S / 汎用 nt）、SINTAX、Kraken 等（ユーザー選択）
- 出力: クラスタごとの分類結果

コアパイプラインは分類 DB に依存しない。

### 5.2 Edit-calling Module（オプション・統合再設計）

現行の `prepAmpliconDB` + `doAlign` + `doEditcall` を **1 つのモジュール** として統合再設計する。コアパイプラインの出力（クラスタごとのコンセンサス）を入力とし、**コンセンサス配列のみ** を参照ゲノムにアラインメントして編集パターンを同定する。

#### 必須入力

- 参照ゲノム配列
- PCR プライマー配列（標的領域の座標特定に使用）
- **ガイド RNA 配列情報**（必須）
- PAM 情報（必要に応じて）
- コアパイプライン出力: サンプル × クラスタごとのコンセンサス + 信頼度メトリクス

#### 現行方式の問題点

| 問題 | 詳細 |
|------|------|
| 全リード BLAST | `doAlign` は amplicon を query、全 basecall リードを DB に BLAST。リード数に比例して遅い |
| amplicon DB 経由の間接比較 | `prepAmpliconDB` で抽出した amplicon を介してリードを位置づけ。コンセンサスを使えば不要 |
| リード単位のノイズ | 個々の ONT リードのエラーが編集型集計に混入 |
| 関数分割による見通しの悪さ | 準備・アラインメント・編集判定が別関数だが、実質は一連の処理 |
| 複数コンセンサス未対応 | 1 サンプル 1 編集型の前提に近い |

#### 新モジュール構成（4 ステップ）

```text
入力: consensus（サンプル × クラスタ）+ 参照ゲノム + プライマー + gRNA

  ⑤ prepRefTargets
       参照ゲノム上にプライマー位置・gRNA 切断サイト座標を定義
           ↓
  ⑥ alignConsensusToRef
       各コンセンサスを参照ゲノムの標的領域にアラインメント（リードは使わない）
           ↓
  ⑦ callEditsFromConsensus
       アラインメント結果から編集パターンをコンセンサス単位で同定
           ↓
  ⑧ selectEditCluster
       複数コンセンサスがあるサンプルで、編集報告の代表クラスタを選択
```

---

#### ⑤ `prepRefTargets()` — 現行 `prepAmpliconDB()` の再設計

**目的:** 参照ゲノム上の編集標的領域を定義する。demultiplex には不要であり、Edit-calling Module 専用の準備ステップ。

**処理:**

1. プライマー配列を参照ゲノムにマッピング（minimap2 `map-ont` または `blastn-short`、初回のみ・リード不要）
2. プライマー間の amplicon 座標を特定（現行 `prepAmpliconDB` と同様のロジック）
3. gRNA 配列情報から標的遺伝子・**切断期待位置（cut site）** を特定
4. 各標的について **野生型（intact）参照配列** を amplicon 全体（またはプライマー間領域）として抽出

**`check_window` は廃止する（後述）。**

**出力（amplicon FASTA ではなく座標・参照配列セット）:**

```text
ref_targets/
  targets.bed          # 標的領域座標（amplicon 範囲）
  intact_sequences.fa  # 各標的の野生型 amplicon 配列
  target_metadata.csv  # gRNA ID, 遺伝子名, PAM 位置, cut site 座標
```

**現行との差分:**

- amplicon DB を「リードアラインメント用クエリ」として作る必要がなくなる
- 出力は座標 + 野生型参照配列。コンセンサスとの比較に直接使う
- プライマーマッピングは参照ゲノムに対する **1 回限りの処理**（リード数に非依存）
- 現行の `check_window`（±N bp の評価ウィンドウ）は採用しない

#### `check_window` 廃止の理由

現行 `doAlign()` における `check_window` は、**生リード** を amplicon に BLAST アラインメントした際に、切断サイト周辺の狭い領域だけを切り出して編集判定することで、**アンプリコン上の他部位のシーケンスエラーを編集変異と誤認しない** ためのパラメータであった。

```r
# 現行: 切断位置 ± check_window bp のみを intact 参照として抽出
edit_site_gr <- GRanges(...,
  ranges = IRanges(start = edit_site$V3 - check_window,
                   end = edit_site$V3 + check_window))
```

コンセンサスベースでは:

- シーケンスエラーはクラスタ内で多数決により大幅に低減されている
- 比較対象はリード 1 本ではなく **コンセンサス 1 本 vs 野生型参照 1 本**
- global/semi-global alignment で diff を直接取得できる

したがって、**ノイズ抑制のための `check_window` は不要** である。

**代わりに必要なもの:**

| 入力 | 役割 |
|------|------|
| **gRNA 切断サイト座標（必須）** | 編集を評価する生物学的アンカー。「どこで編集が起きたか」の報告位置 |
| 野生型 amplicon 参照配列 | コンセンサスとの diff の基準 |

編集判定は **コンセンサス vs 野生型の full alignment** で行い、diff のうち **切断サイトに近接する変異**（例: cut site ± 数 bp 以内の indel/substitution）を編集イベントとして報告する。cut site 近傍の判定閾値は生物学的に固定した小さな値（例: ±3 bp）とし、ユーザーがチューニングする広い `check_window` パラメータとは区別する。

---

#### ⑥ `alignConsensusToRef()` — 現行 `doAlign()` の置き換え

**目的:** クラスタごとのコンセンサス配列を参照ゲノム上の標的領域にアラインメントする。

**入力:** コンセンサス FASTA（サンプル × クラスタ）、`prepRefTargets` の出力

**方針:**

| 項目 | 現行 `doAlign` | 新 `alignConsensusToRef` |
|------|---------------|-------------------------|
| アラインメント対象 | 全 basecall リード | コンセンサス配列のみ |
| 件数 | 数万〜数百万リード | サンプル数 × クラスタ数（通常 < 100） |
| アライナ | BLAST（amplicon → read DB） | minimap2 または `pwalign`（consensus → ref region） |
| 速度 | ボトルネック | コンセンサス数に比例、実質瞬時 |

**推奨アルゴリズム:**

1. **Primary:** `pwalign::pairwiseAlignment`（global/semi-global）でコンセンサス vs 野生型参照配列
   - コンセンサスは既にエラー補正済み → 高精度 global alignment が可能
   - R パッケージ内で完結、依存が少ない
2. **Alternative:** minimap2 でコンセンサスを参照ゲノム全体にマッピング → 標的領域でフィルタ
   - 複数標的・染色体レベルの配置が必要な場合に有効
   - 既に `makeMMI()` で minimap2 インフラあり

**各コンセンサスについて記録するメトリクス:**

- mapping quality / alignment score
- identity, CIGAR
- 切断サイト近傍における indel/substitution の位置と塩基変化
- アラインメント失敗フラグ（コンセンサスが参照と大きく乖離 → コンタミ候補）

---

#### ⑦ `callEditsFromConsensus()` — 現行 `doEditcall()` の再設計

**目的:** アラインメント済みコンセンサスから編集パターンを判定する。

**現行 `doEditcall` から継承するロジック:**

- gRNA 切断サイト近傍の indel/substitution 検出
- 野生型（intact）配列との比較による編集型分類
- 編集型のカテゴリ定義（現行の編集評価レベル）

**現行からの変更点:**

| 項目 | 現行 | 新 |
|------|------|-----|
| 入力単位 | リード単位のアラインメント集計 | コンセンサス 1 本あたり 1 判定 |
| 評価範囲 | `check_window` で狭い領域を切り出し | full alignment の diff から cut site 近傍を抽出 |
| `check_window` | ユーザー指定（デフォルト 10 bp） | **廃止**（ノイズ抑制目的は不要） |
| 集計 | `table(read_seq)` でリード数をカウント | 不要（コンセンサス自体が代表配列） |
| 信頼度 | リード count ベース | クラスタの support / fraction / per-base support を利用 |
| 出力 | サンプル × 遺伝子 × 編集型 | サンプル × クラスタ × 遺伝子 × 編集型 |

**出力例:**

```text
Sample01 / Cluster1
  Target: GeneA
  Edit pattern: -3bp deletion at gRNA site
  vs intact: edited
  Alignment identity: 99.2%
  Cluster support: 1,452 reads (97.1%)

Sample01 / Cluster2
  Target: GeneA
  Edit pattern: intact (no edit detected)
  vs intact: wild-type
  Cluster support: 31 reads (2.1%)  ← コンタミ候補
```

---

#### ⑧ `selectEditCluster()` — 新規

1 サンプルに複数クラスタ（複数コンセンサス）が得られた場合、**サンプルレポート用の代表編集パターン** を選択する。

**スコアリング関数（案）:**

```text
cluster_score =
    w1 * log(cluster_fraction)
  + w2 * consensus_confidence
  + w3 * edit_signal_strength(gRNA_window)
  + w4 * alignment_quality
  - w5 * penalty_if_no_edit_signal_and_small_cluster
```

**考慮要素:**

| 要素 | 重み付けの考え方 |
|------|----------------|
| クラスタサイズ / 割合 | 主クラスタを優先しつつ、小クラスタも棄却しない |
| コンセンサス信頼度 | 低信頼クラスタは編集判定から除外 |
| ガイド RNA 周辺の編集シグナル | 編集が期待される位置に indel/substitution があるクラスタを優先 |
| 参照ゲノムとの整合性 | アラインメント品質が高いクラスタを優先 |
| 野生型との一致 | 編集実験では「編集あり」クラスタを、コントロールでは「野生型」クラスタを優先 |

**出力:**

```text
Sample01
  Primary cluster for edit-calling: Cluster1 (score: 0.92, confidence: high)
    Reason: largest fraction + clear -3bp deletion at gRNA target
  Secondary clusters:
    Cluster2: intact, 2.1% — likely contaminant
  Sample-level edit call: -3bp deletion (based on Cluster1)
```

全クラスタの編集判定結果は保持し、代表選択はレポート・サマリー用とする。

---

## 6. 既存ツール調査と比較

### 6.1 Demultiplex 関連

| ツール | アルゴリズム | 検索対象 | 両端 | 曖昧性処理 | 自作 index |
|--------|------------|---------|------|-----------|-----------|
| **Dorado / Guppy** | 改変 Needleman-Wunsch | ONT キット adapter+barcode | キット依存 | 最高スコア + 閾値 | ❌ |
| **Porechop** | 末端 alignment、% identity | ONT キット barcode | ✅ | threshold + diff | ❌ |
| **cutadapt** | semi-global edit distance | ユーザー定義 | `-g`/`-G` | 最初のマッチ | ✅（短 index に弱い） |
| **Barbell** | edit distance + subsequence scoring | flank + barcode | ✅ | パターン認識 | ✅ |
| **TagGen** | Levenshtein | 設計済み barcode | end モード | 距離閾値 | ✅ |
| **MysteryMaster** | Cola (SW)、位置ウィンドウ | ONT barcode | 位置制約 | SW スコア cutoff | キット向け |
| **本ツール（計画）** | edit distance + Index+Primer | index + primer | ✅ | best/second + 4 段階 | ✅ |

### 6.2 本ツールの差別化

既存 demultiplex ツールは **Demultiplex →（場合により trim）** で終わることが多い。

本ツールは:

1. **未知アンプリコン** の配列推定まで一貫して行う
2. **参照配列不要** のコアパイプライン
3. クラスタごとの **信頼度付きコンセンサス** を出力
4. Taxonomy / 編集解析を **オプションモジュール** として接続
5. 複数コンセンサス環境での **編集代表クラスタ選択**（新規）

### 6.3 既存ツールから取り入れる要素

- **Porechop:** `barcode_diff`（best vs second-best マージン）
- **Barbell:** flank + barcode パターン、subsequence scoring
- **TagGen:** index 設計時の Levenshtein 距離評価
- **cutadapt:** edit distance 閾値設計
- **MiniBar:** Index + Primer 複合探索、edit distance による demux（DuBA.flow / ONTbarcoder でも使用）

### 6.4 関連ツールとの比較（参考）

開発計画策定時に調査した、自作インデックス付き Nanopore アンプリコン解析に関連するツールとの比較。いずれも本ツールと **完全に同一のものではない**。

#### 調査対象ツールの概要

| ツール | URL | 本質 | 主な用途 |
|--------|-----|------|---------|
| **OpenPortablePipeline** | [GitHub](https://github.com/c2997108/OpenPortablePipeline) | NGS 解析の GUI ラッパー／スクリプト集 | メタゲノム、RNA-seq、Nanopore 各種（汎用プラットフォーム） |
| **ONTbarcoder** | [GitHub](https://github.com/asrivathsan/ONTbarcoder) | デュアルタグ amplicon の demux + バーコードコンセンサス | バイオダイバーシティ・DNA バーコーディング（COI 等） |
| **DuBA.flow** | [GitHub](https://github.com/RGSchindler/DuBA.flow) | デュアルバーコード amplicon の end-to-end ワークフロー | 合成 DNA 構築物のバリデーション |

#### 機能比較表

| 機能 | 本ツール（計画） | ONTbarcoder | DuBA.flow | OpenPortablePipeline |
|------|------------------|-------------|-----------|---------------------|
| 自作インデックス demux | ✅ | ✅ | ✅ | △（スクリプト個別） |
| Index + Primer 複合探索 | ✅ | ✅ | ✅（MiniBar） | △ |
| edit distance | ✅ | △ | ✅（MiniBar） | 不明 |
| best/second-best 曖昧判定 | ✅（4 段階） | △ | ❌ | ❌ |
| **未知アンプリコン**（参照不要） | ✅ **コア** | ❌ | ❌ | ❌ |
| サンプル内マルチクラスタリング | ✅ | ❌ | ❌ | △ |
| クラスタごとコンセンサス + 信頼度 | ✅ | △ | △ | △ |
| コンタミ・複数テンプレート検出 | ✅ | ❌ | ❌ | ❌ |
| 参照ゲノム必須 | ❌（コア） | ❌ | **✅ 必須** | スクリプト依存 |
| ゲノム編集パターン解析 | ✅（オプション） | ❌ | ❌ | ❌ |
| 複数コンセンサスから編集代表選択 | ✅（新規） | ❌ | ❌ | ❌ |
| Taxonomy（オプション） | ✅ | ✅（COI 特化） | ❌ | △ |

#### 各ツールの詳細

**DuBA.flow — demux 思想が最も近い**

- [ACS Synthetic Biology 2023](https://doi.org/10.1021/acssynbio.3c00522) で発表。Docker ベースの統合パイプライン
- demultiplex に **MiniBar** を使用。入力形式は Index + Primer の 5 列 TSV:

```text
SampleID | FwIndex | FwPrimer | RvIndex | RvPrimer
```

- MiniBar は read 末端で Index + Primer を **edit distance** で探索（本ツールの demux 案と同型）
- demux 後: minimap2 で **サンプルごとの参照配列** にマッピング → 1 サンプル 1 コンセンサス
- Sanger トレース風の HTML レポートを出力

**限界（本ツールとの差）:** 参照配列必須、単一コンセンサス、未知アンプリコン・コンタミ分離は想定外。

**ONTbarcoder — demux + バーコード特化、用途が異なる**

- [BMC Biology 2021](https://doi.org/10.1186/s12915-021-01141-x)、[bioRxiv 2023](https://doi.org/10.1101/2023.06.26.546538)
- 3 モジュール: (1) demultiplex、(2) barcode calling（コンセンサス）、(3) バーコード比較
- GUI、リアルタイムバーコーディング対応
- COI 等の **既知長・タンパーク質コード遺伝子** 向けに最適化。4 つの QC 基準（長さ、翻訳可能性等）
- 論文内で demux に MiniBar とも比較（ONTbarcoder の方が高精度と報告）

**限界:** バイオダイバーシティ／標本バーコーディングが主目的。1 標本 → 1 バーコード。コンタミ・複数テンプレート分離・ゲノム編集解析は非対応。

**OpenPortablePipeline — 直接競合ではない**

- Java GUI + Docker/Singularity 上で各種 NGS スクリプトを実行する **プラットフォーム**
- Nanopore 関連スクリプト: `nanopore~split-barcode`, `nanopore~get-consensus`, `mapping-nanopore~minimap2` 等が個別に存在
- 統合された unknown-amplicon パイプラインではない

#### ポジショニング

```text
                    参照配列
                    不要 ←────────────────────────→ 必須
                     │                              │
  単一コンセンサス   │  DuBA.flow                   │
       ↑             │  (合成DNA validation)        │
       │             │                              │
       │   ONTbarcoder                              │
       │   (COI barcoding)                          │
       │             │                              │
       │             │    ★ 本ツール                │
       │             │    ・未知アンプリコン         │
       │             │    ・マルチクラスタ          │
       │             │    ・信頼度評価              │
       │             │    ・編集解析オプション       │
       ↓             │                              │
  マルチコンセンサス │                              │
  + 信頼度           │                              │
```

**推奨ポジショニング:**

> 自作インデックス付き Nanopore PCR データから、**参照配列なしに**増幅領域の真の配列を推定し、コンタミ・複数テンプレートを分離したうえで、Taxonomy またはゲノム編集解析につなげる統合パイプライン

- DuBA.flow = 参照ありの合成 DNA 検証
- ONTbarcoder = COI バーコーディング
- 本ツール = **未知アンプリコンの de novo 配列推定 + 下流解析**

### 6.5 優位性・新規性の評価（参考）

#### 明確にある優位性

1. **参照配列不要のコアパイプライン** — 既知はインデックス + プライマーのみ、増幅領域は未知
2. **サンプル内マルチクラスタリング + 信頼度評価** — 主アンプリコン / コンタミ / chimera / 複数編集パターンの共存検出
3. **下流モジュールの柔軟性** — Taxonomy またはゲノム編集解析（ユーザー選択）、`selectEditCluster()` による複数コンセンサスからの代表選択
4. **demux の信頼度分類** — High / Low / Ambiguous / Unclassified、best vs second-best スコア差（MiniBar には標準でない）

#### 新規性が限定的な部分（正直な評価）

| 要素 | 評価 |
|------|------|
| Index + Primer 複合 demux | MiniBar（DuBA.flow, ONTbarcoder でも使用）で **既報** |
| edit distance による demux | MiniBar の `-e` / `-E` で **既報** |
| 両端インデックス探索 | MiniBar Method 2、ONTbarcoder で **既報** |
| コンセンサス生成 | 各ツール・各所で **既報** |

→ **demux アルゴリズム単体では論文レベルの新規性は出にくい**。差別化の中心は demux 以降のパイプライン全体にある。

#### 総合評価

| 観点 | 評価 |
|------|------|
| demux 単体の新規性 | 低〜中（MiniBar と重複。best/second-best 分類が差分） |
| パイプライン全体の新規性 | **高**（未知アンプリコン + マルチクラスタ + 信頼度 + 編集解析） |
| 実用優位性 | **高**（既存ツールの組み合わせでは再現しにくいワークフロー） |
| 学術的新規性（論文化） | demux 単体では弱い。パイプライン設計 + ベンチマークで十分狙える |

#### demux 実装への示唆

MiniBar（[calacademy-research/minibar](https://github.com/calacademy-research/minibar)）の Index + Primer + edit distance は実績がある。ゼロから作るより **改善・拡張** が現実的:

- MiniBar 相当の探索をベースに
- **best/second-best + 4 段階分類** を追加
- その後の **clustering → multi-consensus** パイプラインと統合

---

## 7. 現行コードの改修方針

### 7.1 削除・置き換え対象

| 関数 / 機能 | 方針 |
|------------|------|
| `doDemultiplex()` | `doDemultiplex2()` に置き換え |
| `prepAmpliconDB()` | Edit-calling Module 内の `prepRefTargets()` に統合・再設計 |
| `doAlign()` | Edit-calling Module 内の `alignConsensusToRef()` に置き換え（全リード BLAST を廃止） |
| `doEditcall()` | Edit-calling Module 内の `callEditsFromConsensus()` に再設計 |
| `miaoEditcall()` | 新パイプライン `miaoPipeline()` に再構成（コア + オプションモジュール） |

### 7.2 維持・改修対象

| 関数 / 機能 | 方針 |
|------------|------|
| `doEditcall()` の編集型分類ロジック | `callEditsFromConsensus()` に移植・改修 |
| `editViewer()` | 維持・拡張（クラスタ別編集結果の可視化） |
| `doBasecall()` | 維持（Dorado ラッパー） |
| `evalMiao()` | 改修（新メトリクス・クラスタ情報に対応） |

### 7.3 新規実装

#### コアパイプライン

| 関数（仮名） | 役割 |
|-------------|------|
| `doDemultiplex2()` | Index+Primer、edit distance、4 段階分類 |
| `doCluster()` | サンプル内 read clustering（3 パイプライン切替: meshclust / mmseqs2 / umap_meshclust） |
| `doConsensus()` | クラスタごとのコンセンサス + pileup |
| `evalConfidence()` | 信頼度メトリクス算出 |
| `miaoPipeline()` | コアパイプラインのオーケストレーション |

#### Edit-calling Module（`prepAmpliconDB` + `doAlign` + `doEditcall` の統合後継）

| 関数（仮名） | 現行との対応 | 役割 |
|-------------|------------|------|
| `prepRefTargets()` | `prepAmpliconDB` | 参照ゲノム上の標的座標・野生型配列の定義 |
| `alignConsensusToRef()` | `doAlign` | コンセンサス → 参照ゲノムのアラインメント |
| `callEditsFromConsensus()` | `doEditcall` | コンセンサス単位の編集パターン判定 |
| `selectEditCluster()` | （新規） | 複数コンセンサスから編集報告の代表を選択 |
| `doEditCalling()` | `miaoEditcall` の編集部分 | 上記 4 関数のオーケストレーション |

#### その他

| 関数（仮名） | 役割 |
|-------------|------|
| `doTaxonomy()` | オプション taxonomy モジュール |
| `evaluateIndexSet()` | index 設計評価（最小 edit distance 等） |

---

## 8. 開発ロードマップ

### Phase 1: Demultiplex 刷新（最優先）

- [ ] `doDemultiplex2()` プロトタイプ
  - Index + Primer クエリ構築
  - 末端ウィンドウ semi-global alignment（`pwalign`）
  - edit distance ベーススコアリング
  - 4 段階分類 + TSV 出力
- [ ] 現行 BLAST 方式とのベンチマーク
- [ ] パラメータ最適化（T_high, T_low, δ, 末端ウィンドウ長）

### Phase 2: Clustering + Consensus

- [ ] `doCluster()` — 3 パイプライン切替実装
  - [ ] `"meshclust"`: meshclust ラッパー（本番デフォルト）
  - [ ] `"mmseqs2"`: MMseqs2 `easy-linclust` ラッパー（高速プロトタイプ用）
  - [ ] `"umap_meshclust"`: UMAPclust → meshclust 統合パイプライン
  - [ ] 共通 Stage 2: minimap2 リファインメント（`refine` パラメータで ON/OFF）
- [ ] 3 パイプライン × Stage 2 有無のベンチマーク（模擬データ）
- [ ] `doConsensus()` — pileup consensus
- [ ] `evalConfidence()` — 信頼度メトリクス
- [ ] 統合パイプライン `miaoPipeline()`（仮名）

### Phase 3: Edit-calling Module（統合再設計）

- [ ] `prepRefTargets()` — 現行 `prepAmpliconDB` の再設計
  - プライマー → 参照ゲノムマッピング（1 回限り）
  - gRNA 切断サイト座標の定義（`check_window` は廃止）
  - 野生型参照配列の抽出
- [ ] `alignConsensusToRef()` — 現行 `doAlign` の置き換え
  - コンセンサス配列のみを `pwalign` / minimap2 でアラインメント
  - 全リード BLAST の廃止
  - 速度・精度のベンチマーク（現行 vs 新方式）
- [ ] `callEditsFromConsensus()` — 現行 `doEditcall` の再設計
  - コンセンサス単位の編集パターン判定
  - クラスタ信頼度メトリクスの統合
- [ ] `selectEditCluster()` — 新規
  - 複数コンセンサスからの代表選択スコアリング
- [ ] `doEditCalling()` — モジュール統合オーケストレーション
- [ ] `editViewer()` 拡張 — クラスタ別編集結果の可視化

### Phase 4: Taxonomy + QC・ユーティリティ

- [ ] `doTaxonomy()` — オプション分類
- [ ] `evaluateIndexSet()` — index 設計評価
- [ ] 品質レポート自動生成（demux 率、クラスタ数、未分類率、コンタミ候補）
- [ ] chimera 検出（オプション）

---

## 9. 技術的決定事項（要検討）

| 項目 | 候補 | 現時点の推奨 |
|------|------|-------------|
| 末端 alignment（demux） | pwalign, Edlib (Rcpp), vmatchPattern | pwalign（既存依存）でプロトタイプ |
| Clustering | meshclust, MMseqs2, UMAPclust+meshclust（切替可能） | 本番: `"meshclust"`、高精度: `"umap_meshclust"`、開発: `"mmseqs2"` + 共通 minimap2 リファインメント |
| Consensus polish | なし, Racon, Medaka | まず pileup のみ、必要に応じて Medaka |
| コンセンサス → 参照アラインメント | pwalign (global), minimap2 | pwalign（コンセンサス同士の比較に最適） |
| プライマー → 参照ゲノムマッピング | minimap2, blastn-short | minimap2（`makeMMI` インフラ活用） |
| 編集代表クラスタ選択 | 新規スコアリング関数 | Phase 3 で設計・検証 |
| パッケージ名 | miaoseq 継続 or 新名称 | 要議論 |

---

## 10. 参考リンク

### Demultiplex 一般

- [Dorado barcoding algorithm](https://nanoporetech.com/document/data-analysis) — 改変 Needleman-Wunsch
- [Porechop](https://github.com/rrwick/Porechop) — barcode_threshold + barcode_diff
- [cutadapt algorithms](https://cutadapt.readthedocs.io/en/stable/algorithms.html) — edit distance
- [Barbell](https://github.com/rickbeeloo/barbell) — pattern-aware demux, custom flanks
- [TagGen](https://f1000research.com/articles/15-642) — Levenshtein barcode design + demux
- [MysteryMaster](https://pmc.ncbi.nlm.nih.gov/articles/PMC12487470/) — Cola aligner, position windows
- [MiniBar](https://github.com/calacademy-research/minibar) — Index + Primer + edit distance demux

### 関連ツール（§6.4 参照）

- [OpenPortablePipeline](https://github.com/c2997108/OpenPortablePipeline) — NGS 解析 GUI プラットフォーム
- [ONTbarcoder](https://github.com/asrivathsan/ONTbarcoder) — [論文 (BMC Biology 2021)](https://doi.org/10.1186/s12915-021-01141-x)、[ONTbarcoder 2.0 (bioRxiv 2023)](https://doi.org/10.1101/2023.06.26.546538)
- [DuBA.flow](https://github.com/RGSchindler/DuBA.flow) — [論文 (ACS Synth Biol 2023)](https://doi.org/10.1021/acssynbio.3c00522)

### クラスタリング

- [LACA](https://github.com/yanhui09/laca) — [論文 (Gut Microbes 2025)](https://doi.org/10.1080/19490976.2025.2516703)
---

## 11. 用語整理

| 用語 | 定義（本ツールにおいて） |
|------|------------------------|
| **インデックス** | サンプル識別用の既知配列（PCR アダプター等） |
| **プライマー** | PCR 増幅に使用する既知配列 |
| **増幅対象領域** | プライマー内側の **未知** 配列 |
| **コンセンサス** | クラスタ内リードから推定した代表配列 |
| **コアパイプライン** | 参照ゲノム・分類 DB 不要の処理（demux → cluster → consensus → confidence） |
| **Edit-calling Module** | 参照ゲノム + gRNA を用いた編集パターン同定（`prepRefTargets` → `alignConsensusToRef` → `callEditsFromConsensus` → `selectEditCluster`） |
| **下流モジュール** | Taxonomy または Edit-calling（ユーザー選択） |
