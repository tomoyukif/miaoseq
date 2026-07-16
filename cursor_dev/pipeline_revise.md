# パイプライン再設計計画

作成日: 2026-07-08  
改訂:

- 2026-07-09 — Phase C0–D 完了を反映。`miaoEditcall` は新フローでは不要（廃止予定）
- 2026-07-09 — **Phase B 完了**: basecall / `miaoEditcall` / `doAlign` / BLAST 撤去。`prepAmpliconDB` を `matchPattern` 実装に置換
- 2026-07-09 — §15: 類似度クラスタリングによる editcall 偽陽性（シーケンシングエラー集約）の問題と対策案を追記
- 2026-07-09 — §15.3.G + Phase G: `assembly_backend = "overlap_graph"`。当時の `allele_source = "reads"` は Phase I で経路 B 本線化（切替廃止）
- 2026-07-14 — §16 / Phase H: 最終クラスタ未寄与リードの再検証（コンセンサス同士・未割当リードの類似度）を開発項目として追加
- 2026-07-14 — §15.3.A を確定・実装: パイプラインを **Assemble 経路** と **Editcall 経路** に分岐。`doEditcall(gene_assign=…)` は **reads 専用**（`allele_source` 切替なし）。クラスタ純度フィルタ（§15.3.E）は Assemble に実装（Phase I 完了）
- 2026-07-14 — §15.7: 大 indel 向け **適応編集窓** を実装（`anchor_bp` / `max_expand`、`ref_seq` 列、アライン失敗リード棄却）
- 2026-07-14 — §16 Phase H: バックエンド比較拡張（`edlib` / `vsearch` / `blastn` / `mmseqs` → `summary_compare.tsv`）
- 2026-07-14 — 監査対応: `consensus_backend` を **spoa 固定**（majority/quality 廃止）; amplicon-ref fallback 削除; 期待長はトリム後 insert; `prepAmpliconDB` に multi-hit/長区間停止 + `amplicon_fasta`; demux の `check_prefix`/`high_confidence`/`require_prefix` 削除; indel ラベルを gap bp 厳密カウントに修正
- 2026-07-14 — §17: 監査「中程度」所見への対応方針を確定（`require_unique_pair` 削除、組合せ demux は利用者責任、vsearch↔edit は要検討、`fraction_sample` は gene バケツとして文書化、PAM 染色体は genome 一致・複数ガイドは案 A）
- 2026-07-14 — §7.4 / §17: 同一 amplicon 複数ガイドの **案 A′（ジョイント両 cut 検証）** を採用（ガイド別窓＋介在欠失イベント層）
- 2026-07-14 — §18: Phase J 実装前の **未決事項の説明＋回答欄** を追加。関連 md（demux_revise / editcall_adapter / ampliconresolve_plan）と README を現行方針に同期
- 2026-07-14 — §18.3: J1b–J11b 記入レビュー。追加回答欄 **J1c**（窓中心 vs cut・列5）のみ残す
- 2026-07-14 — **Phase J 実装完了**: `require_unique_pair` 削除、PAM seqname/`matchPattern` 写像・列5 strand、案 A+A′（`min_span_bp`/`excision_tol_bp`・joint CSV）、extdata・README・単体テスト。vsearch↔edit は文書一行のみ（J11b 案2）
- 2026-07-15 — §19 Phase K: クラスタリング/コンセンサス基盤の刷新。`cluster_backend = "internal"`（edlib greedy）を廃止し `"mmseqs"`（connected component）を追加。`spoa` を `abpoa` に置換。Racon 1パスポリッシュを追加
- 2026-07-16 — §19 Phase K 確定: **正式採用は** `vsearch` **+** `abpoa`**（**`-Q` **なし・racon なし）**。`mmseqs` / linclust / racon / abPOA `-Q` は評価のうえ廃止
- 2026-07-16 — §19.8 / §6: 精査 `20260716_1056` 反映。`min_cluster_identity`（既定 0.95）へ置換、`max_clusters=Inf`、`fraction_sample`→`fraction_bucket`、`clusters.fasta`=メンバー配列、run_stats 拡充、文書同期
関連:
- [cursor_dev/demux_revise.md](demux_revise.md) — demux 詳細
- [cursor_dev/ampliconresolve_plan.md](ampliconresolve_plan.md) — Step 2 / Step 3 詳細
- [cursor_dev/editcall_adapter.md](editcall_adapter.md) — editcall 接続メモ
- [cursor_dev/ontbarcoder_demux.md](ontbarcoder_demux.md) — ONTbarcoder 参照
- [cursor_dev/code_reorganize_plan.md](code_reorganize_plan.md) — 過去のファイル分割（参考）

---



## 0. 実装状況サマリ（2026-07-16）

- **最終確認**: 2026-07-16
- **対象コード**: `R/demultiplex.R`, `R/editcall.R`, `R/amplicon_resolve.R`, `man/doDemultiplex.Rd`, `man/doEditcall.Rd`, `inst/extdata/agr8_pam_list.csv`, `tests/test-editcall-phase-j.R`, `README.md`

**結論:** 本線は **2 経路に分岐**する。アンプリコン代表配列の復元（Assemble）と編集頻度推定（Editcall）は、§15.2 のとおり **クラスタリング設定を共有しない**。レガシー（`miaoEditcall` / `doAlign` / `doBasecall` / BLAST）は Phase B で削除済み。**Phase J（§17–§18）も実装済み。** **Phase K（§19）は確定:** Assemble は **vsearch + abpoa**（他バックエンド廃止）。


| Phase  | 内容                                                                  | 状態                                       |
| ------ | ------------------------------------------------------------------- | ---------------------------------------- |
| A      | Demux 刷新（edlib）                                                     | **完了**                                   |
| **B**  | basecall / レガシー / BLAST 撤去                                          | **完了**                                   |
| C      | Amplicon resolve（C0–C3）                                             | **完了**                                   |
| D      | Editcall 初期（`amplicon_out` 消費）                                      | **完了**（API は Phase I で再設計）               |
| G      | `overlap_graph` + per-read editcall 原型                              | **完了**                                   |
| **I**  | **2 経路分離 +** `gene_assign`**→**`doEditcall` **+ 純度フィルタ（§15.3.A/E）** | **完了**                                   |
| **I′** | **§15.7 適応編集窓**（`anchor_bp` / `max_expand`）                         | **完了**                                   |
| **H**  | 最終クラスタ未寄与リードの再検証（§16）MVP                                            | **完了**（edlib・Q1/Q2・`amplicon_reassess/`） |
| **J**  | 監査中程度対応（§17–§18）: demux API・PAM / 案 A+A′・文書注意（vsearch）              | **完了**                                   |
| **K**  | クラスタリング/コンセンサス刷新（§19）: **vsearch + abpoa**（§19.8 精査反映済み）              | **完了**                                   |




### 計画からの主な同期内容（Phase J）

- `require_unique_pair` 削除済み（`R/demultiplex.R` / man）
- PAM: seqname 完全一致、`matchPattern` 写像、列5 strand、cut 中心窓（`R/editcall.R`）
- 案 A（ガイド別窓）+ 案 A′（`editcall_joint*.csv`、`min_span_bp=30` / `excision_tol_bp=20`）
- `inst/extdata/agr8_pam_list.csv` を `chr01` 等形式へ更新
- 単体テスト `tests/test-editcall-phase-j.R`
- vsearch↔絶対 edit: README「Clustering identity note」1行のみ（技術揃えは別 Phase・todo）



### 差分・注意（定義は未変更）

- 案 A の「guide 欠落時 `gene_1`/`gene_2` 自動採番」は**未実装**。複数行 gene は列4 guide **必須**で停止（§18 回答どおりの厳格側）
- `editViewer` / `evalMiao` の joint 集計・plate 表示は意図どおり Phase J 外（todo）
- §3.1 の `doReassessAssemblies` 表記を「未着手」から Phase H MVP **完了**に修正（計画ずれ）



### 0.1 推奨フロー（2 経路）

```
[ユーザー] Dorado 等 → FASTQ
    ↓
doDemultiplex → demultiplex/
    ↓
doAssignGenes → amplicon_assign/gene_assignments.tsv
    ↓
    ├────【経路 A: Assemble】─────────────────────────────
    │    doAssembleAmplicons（+ クラスタ純度フィルタ §15.3.E）
    │      → amplicon/{sample}/consensus.fasta, cluster_counts.tsv
    │    （任意）doReassessAssemblies（§16）
    │
    └────【経路 B: Editcall】─────────────────────────────
         doEditcall(gene_assign = …)   # §15.7 + Phase J（案 A/A′）
           → editcall/editcall_*.csv（+ joint CSV・複数ガイド時）
         evalMiao / editViewer
```


| 経路             | 目的                      | 主入力                               | 主出力         |
| -------------- | ----------------------- | --------------------------------- | ----------- |
| **A Assemble** | フル長アンプリコン代表・QC          | `gene_assign`                     | `amplicon/` |
| **B Editcall** | cut site 周辺の変異頻度（リード単位） | `gene_assign` + FASTQ / by_sample | `editcall/` |


- 経路 B は `doAssembleAmplicons` **を必須としない**。
- Editcall は常にリード単位 edit-window 集計（Assemble consensus は用いない）。
- 経路 A: `cluster_backend = "vsearch"`（§19 Phase K 確定）+ **abpoa コンセンサス**（`-Q` / racon / mmseqs / spoa / internal は廃止）
- 経路 B: `doEditcall`（per-read; §15.3.A）。大 indel 向けの適応編集窓は **§15.7（実装済み）**。複数ガイドは **Phase J（案 A+A′）**
- Gene 割当は **プライマー edlib のみ**（amplicon NW fallback は削除）。`no_ref` / 非 `assigned` は Assemble/Editcall では使わない
- `amplicon_fn` があるときの期待長フィルタは **トリム後 insert 長**（生リード長ではない）
- `prepAmpliconDB`: ゲノム導出時は multi-hit または span > `2 * expected_length` で停止。任意で `amplicon_fasta`（配列名 = gene ID）を直接入力可
- Demux: prefix 検証（`check_prefix` / `require_prefix` / `high_confidence`）は廃止。組合せプレートで `allow_single_end=TRUE` は起動時停止。**組合せ dual の誤割当リスクは利用者責任**（§5.1 / §17）。`require_unique_pair` は **削除済み**
- Assemble 出力の `fraction_bucket` は **サンプル全体ではなく sample×gene バケツ内比率**（旧名 `fraction_sample`。§6.3）。列名は意図的に直観非依存とし README で定義
- 経路の役割: **Assemble** = フル長代表の復元（クラスタ内は均される）。**Editcall** = 既知 cut 周辺の局所パターン。16S 等の多様領域の全体像には Assemble；ゲノム編集の局所頻度には Editcall。Assemble の不確実性は ONT アンプリコン解析に共通で完全解決はない（§19.8）
- Editcall PAM: 染色体名は genome FASTA の seqname と **一致必須**（`chr%02d` 付与は **廃止済み**）。同一 amplicon 上の複数 PAM は **案 A**（ガイド別窓）+ **案 A′**（ジョイント両 cut／介在欠失）（§7.4 / §17）**実装済み**

---



## 1. 目的

miaoseq を「R から basecall まで全部やる巨大ラッパー」から、**ユーザーが用意した FASTQ を起点に、demux → gene 割当 →（Assemble または Editcall）** するパッケージへ再設計する。

**2026-07-14:** §15 を受け、Assemble と Editcall を直列必須にしない（類似度クラスタ経由の偽陽性集約を回避）。

```
[ユーザー] Dorado などで basecall → FASTQ 準備
    ↓
[miaoseq]  doDemultiplex
    ↓ 任意
[miaoseq]  splitDemultiplexReads
    ↓
[miaoseq]  doAssignGenes          … 両経路の共通前提
    ↓
    ├── 経路 A: doAssembleAmplicons → アンプリコン配列
    └── 経路 B: doEditcall → 編集頻度 → evalMiao / editViewer
```

---



## 2. 旧パイプラインの問題（解消済み）

旧 `miaoEditcall()` — **新規解析では使わない**:

```
doBasecall → doDemultiplex(BLAST) → doAlign(BLAST) → doEditcall
```


| 問題                          | 旧              | 新                                      |
| --------------------------- | -------------- | -------------------------------------- |
| basecall を R が抱え込む          | `doBasecall`   | ユーザーが Dorado 等で FASTQ 準備               |
| demux と amplicon assign が独立 | BLAST 巨大表      | `assignments.tsv` + サンプル単位 resolve     |
| サンプル別配列がない                  | 全リード一括 BLAST   | `amplicon/{sample_id}/consensus.fasta` |
| editcall が read 単位 join     | demult × align | `doEditcall(amplicon_out)`             |


---



## 3. 新パイプライン概要



### 3.1 ステップ構成（2 経路）


| Step     | 関数                        | 入力                                                                               | 主出力                                     | 状態                           |
| -------- | ------------------------- | -------------------------------------------------------------------------------- | --------------------------------------- | ---------------------------- |
| 0（R外）    | Dorado 等                  | pod5                                                                             | FASTQ                                   | パッケージ対象外                     |
| 1        | `doDemultiplex`           | FASTQ, index_list                                                                | `assignments.tsv` 等                     | **完了**                       |
| 1b（任意）   | `splitDemultiplexReads`   | FASTQ + assignments                                                              | `by_sample/*.fq.gz`                     | **完了**                       |
| 2a       | `doAssignGenes`           | assignments + FASTQ / by_sample                                                  | `amplicon_assign/gene_assignments.tsv`  | **完了**（両経路の共通）               |
| **2b-A** | `doAssembleAmplicons`     | `gene_assign`                                                                    | `amplicon/{sample}/` consensus・clusters | **完了**（+ 純度フィルタ Phase I）     |
| **2b-B** | `doEditcall`（per-read）    | `gene_assign` + PAM / genome / amplicon ref + FASTQ                              | `editcall/`（+ joint CSV）                | **完了**（Assemble 非依存・Phase J） |
| 2c（任意）   | `doReassessAssemblies`    | Assemble 出力 + 未寄与リード                                                             | 類似度 / 再割当候補                             | **完了**（§16 MVP）              |
| —        | `evalMiao` / `editViewer` | `out_dir`（`demultiplex/` + `amplicon_assign/` および/または `amplicon/` + `editcall/`） | レポート / Viewer                           | **完了**（`amplicon_assign` 優先） |
| —        | `prepAmpliconDB`          | genome + primers                                                                 | `ref/amplicon.fa`                       | editcall 用・継続                |


**分岐規則:** Step 2a の後、経路 A（2b-A）と経路 B（2b-B）は独立に実行できる。Editcall 本線は 2b-B であり、2b-A を前提としない。

### 3.2 廃止済み（Phase B で削除）


| 対象                           | 扱い                            |
| ---------------------------- | ----------------------------- |
| `miaoEditcall()`             | **削除**（`R/pipeline.R`）        |
| `doAlign()`                  | **削除**（`R/amplicon_assign.R`） |
| `doBasecall()` / `makeMMI()` | **削除**（`R/basecall.R`）        |
| `R/blast_utils.R`            | **削除**（BLAST 依存なし）            |
| demux の `blast_path`         | 以前より削除済み                      |
| 旧 `demultiplex_list.csv` 列互換 | 維持しない                         |


`prepAmpliconDB()` は **Biostrings::matchPattern**（両鎖完全一致）で primer 位置を特定。BLAST 不要。

---



## 4. ディレクトリ構造（新）

```
{out_dir}/
  demultiplex/
    assignments.tsv
    summary_by_sample.tsv
    unassigned.tsv
    index_layout.tsv
    barcode_conflicts.tsv
    design_check.tsv
    by_sample/                    # 任意
      {sample_id}.fq.gz
  amplicon/                       # Step 2
    summary_by_sample.tsv
    gene_assignments.tsv
    {sample_id}/
      consensus.fasta
      clusters.fasta
      cluster_counts.tsv
      stats.tsv
  editcall/                       # Step 3（任意）
    editcall_all.csv
    editcall_filtered.csv
    editcall_summary.csv
    intact_seq.fa
    editcall_joint.csv            # 複数ガイド時（§7.4 案 A′・実装済み）
    editcall_joint_summary.csv    # 同上
  ref/                            # prepAmpliconDB 出力（任意）
    amplicon.fa
  report/                         # 任意（evalMiao / editViewer）
    ...
```

入力 FASTQ は `{out_dir}` 外でもよい（ユーザーがパスを渡す）。

---



## 5. Step 1: Demultiplex

詳細は [demux_revise.md](demux_revise.md)。監査中程度所見への対応は **§17**。

### 契約

```r
# 入力
fastq        # character: 1つ以上の FASTQ パス
demult_dir
index_list
sample_list = NULL
split_reads = FALSE
# ※ require_unique_pair は削除済み（C++ に渡らず no-op だった）

# 出力（常時）
assignments.tsv
summary_by_sample.tsv
unassigned.tsv

# 出力（split_reads=TRUE または後から splitDemultiplexReads）
by_sample/{sample_id}.fq.gz
```



### 5.1 インデックスレイアウト方針（決定）


| 方針      | 内容                                                                  |
| ------- | ------------------------------------------------------------------- |
| 推奨      | **1:1**（各 F/R ID が一意）またはバーコード間 Hamming ≥ `2 * max_barcode_edit + 1` |
| 組合せプレート | 利用可。ただし「誤 dual が **別の合法ウェル**に入る」リスクは **利用者が負う**                     |
| ツール側    | 配列デザインの代わりに topology 検査で誤割当を潰すことはしない（suffix/prefix 厳格化は任意チューニング）    |
| 単端      | 組合せで `allow_single_end=TRUE` は起動時停止（実装済み）                           |


README にも「組合せレイアウトでは誤 dual が合法ウェルに入りうる」と明記する（反映済み）。

### `splitDemultiplexReads()` の位置づけ

```r
# demux と同時
doDemultiplex(..., split_reads = TRUE)

# あとから
doDemultiplex(..., split_reads = FALSE)
splitDemultiplexReads(
  fastq       = fastq,
  assignments = file.path(demult_dir, "assignments.tsv"),
  out_dir     = file.path(demult_dir, "by_sample")
)
```

下流 Step 2 は **どちらでも動く**こと:

1. assignments + 元 FASTQ からサンプル単位でオンデマンド抽出
2. 既存 `by_sample/*.fq.gz` を直接読む（HPC / 外部ツール向け）

---



## 6. Step 2: gene割当 + サンプル別アンプリコン再構成 ✅ 完了

旧 `doAlign`（全リード × amplicon DB BLAST）を **サンプル単位の配列再構成**に置換済み。詳細は [ampliconresolve_plan.md](ampliconresolve_plan.md)。

### 6.1 採用方針

**C. ref-aware** — primer edlib で gene 割当後にクラスタ + コンセンサス。

**正式採用（Phase K 確定; §19）:** `cluster_backend = "vsearch"` + **abpoa コンセンサス**（`-Q` なし・racon なし）。`mmseqs` / linclust / racon / abPOA `-Q` / `internal` / `spoa` は廃止。

**精査 20260716 反映（§19.8）:**
- 閾値: `min_cluster_identity`（既定 **0.95**）。旧 `max_cluster_edit`（edit→%id 換算）は廃止
- `max_clusters` 既定 **`Inf`**（全クラスタ候補を保持）。マイナー除外は `min_cluster_reads`（既定 5）のみ
- 列名: `fraction_sample` → **`fraction_bucket`**（sample×gene バケツ内比率）
- `clusters.fasta`: コンセンサスの複製ではなく **クラスタメンバー配列**（ヘッダ `>{read_id} {cluster_id}`）
- `run_stats.txt`: R / miaoseq 版 + 主要パラメータを記録

### 6.2 API（実装済み）

```r
doAssignGenes(
  assignments,
  out_dir,
  fastq = NULL,
  sample_fastq_dir = NULL,
  primer_list = NULL,
  amplicon_fn = NULL,          # 期待長フィルタ用（幅のみ; fallback 割当は無し）
  max_primer_edit = 5L,
  end_window = 150L,
  length_tolerance = 0.25,    # トリム後 insert vs amplicon.fa 幅
  samples = NULL,
  n_core = 1L,
  overwrite = FALSE,
  stats_unassign = FALSE
)

doAssembleAmplicons(
  gene_assign,
  out_dir,
  primer_list = NULL,
  method = c("both", "cluster", "consensus"),
  min_reads = 5L,
  min_cluster_reads = 5L,       # これ未満のクラスタは除外（マイナー除外の主レバー）
  max_clusters = Inf,           # サイズ上位打ち切りなし（Inf）
  min_cluster_identity = 0.95,  # vsearch --id（旧 max_cluster_edit 廃止）
  max_primer_edit = 5L,
  n_core = 1L,
  cluster_backend = "vsearch",  # 固定
  assembly_backend = c("cluster", "overlap_graph"),
  overlap_min_identity = 0.90,
  strict_end_trim = FALSE,
  min_cluster_purity = NULL,  # デフォルト無効; 出力列 cluster_purity のみ
  overwrite = FALSE
)
```

**Phase K + 精査反映:** クラスタリングは **vsearch のみ**（`--id = min_cluster_identity`）。コンセンサスは **abpoa（FASTA・`-Q` なし）のみ**。racon なし。qual はコンセンサスに渡さない。

### 6.3 出力契約

```
amplicon/{sample_id}/
  consensus.fasta         # >{gene_id}|cluster_{k};size={n};sample={id}
  clusters.fasta          # メンバー配列; ヘッダ >{read_id} {cluster_id}
  cluster_counts.tsv    # sample_id, gene_id, cluster_id, n_reads,
                        # fraction, fraction_bucket, n_reads_gene,
                        # method, cluster_purity, ...
  stats.tsv
  unassigned_to_cluster.tsv
```



#### `fraction` / `fraction_bucket`（誤解防止）


| 列                 | 意味（実装）                                               |
| ----------------- | --------------------------------------------------------- |
| `fraction`        | 同 sample×gene でクラスタ通過後のアレル比率（分母は残クラスタ合計 ≈ `n_reads_gene`） |
| `fraction_bucket` | **サンプル全体ではなく**、同 sample×gene バケツへの入力リード数に対する比率（旧名 `fraction_sample`） |


列名 `fraction_bucket` は意図的に直観非依存（README 確認を促す）。旧 `fraction_sample` は廃止（後方互換なし）。

#### クラスタ閾値

正式バックエンドは **vsearch のみ**（`--id = min_cluster_identity`、`--iddef 2`）。既定 0.95。旧 `max_cluster_edit`（中央値長からの換算）は **廃止**。

`max_clusters = Inf` が既定。クラスタ数の上限打ち切りは行わず、`min_cluster_reads`（既定 5）未満のみ除外する。網羅的 OTU 定量ではなく、十分な支持のある代表復元が目的。

旧 `internal` / `mmseqs` クラスタは **Phase K で廃止**（§19.3）。

---



## 7. Editcall（経路 B）— Phase I 再設計 + Phase J（複数ガイド）



### 7.1 契約（確定）


| 項目    | 内容                                                                                                                                                                                                                                 |
| ----- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 本線    | `doAssignGenes` → `doEditcall` → `evalMiao` / `editViewer`                                                                                                                                                                         |
| 主入力   | `gene_assign`（`GeneAssignResult` / `amplicon_assign/` / `gene_assignments.tsv`）                                                                                                                                                    |
| 必須参照  | `pam_list`, `genome_fn`, `amplicon_fn`, `primer_list`                                                                                                                                                                              |
| 配列ソース | `fastq` または `sample_fastq_dir`（`by_sample/`）                                                                                                                                                                                       |
| 窓抽出   | §15.7 適応編集窓（初期 `± check_window`、両端 `anchor_bp` 一致、失敗時 `max_expand` 棄却）。窓中心は **cut**（Phase J）                                                                                                                                       |
| 出力    | `editcall/editcall_all.csv`（`read_seq`, `ref_seq`, `intact`, `allele_source`, …）, `editcall_filtered.csv`, `editcall_summary.csv`, `intact_seq.fa`。複数ガイド時は `editcall_joint.csv` / `editcall_joint_summary.csv`（§7.4 案 A′・**実装済み**） |


```
gene_assignments（assigned）+ FASTQ
  → primer-trim → 参照 insert へ global align
  → 適応編集窓で read_seq / ref_seq を抽出（§15.7）
     ※ gene に紐づく PAM が複数ならガイドごとに繰り返し（§7.4 案 A）
     ※ さらに隣接ガイド対でジョイント両 cut／介在欠失を判定（§7.4 案 A′）
  → intact = identical(read_seq, ref_seq)
  → sample × target_gene × read_seq の完全一致集計
  → （複数ガイド時）editcall_joint.* に event_class 集計
  → min_count / filtered ルール → summary
```



### 7.2 API

```r
doEditcall(
  gene_assign = NULL,           # 主入力（必須; または amplicon_out から解決）
  amplicon_out = NULL,          # 移行互換: gene_assignments のフォールバック
  pam_list,
  genome_fn,
  amplicon_fn,
  primer_list,
  editcall_dir,
  sample_list = NULL,
  check_window = 10L,
  anchor_bp = 5L,
  max_expand = 50L,
  min_span_bp = 30L,            # Plan A′: excision 判定の最小期待 cut–cut span
  excision_tol_bp = 20L,        # Plan A′: |del_span - expected| の許容（bp 固定）
  min_count = 5L,
  fastq = NULL,
  sample_fastq_dir = NULL,
  max_primer_edit = 5L,
  n_core = 1L
)
```

`gene_assign` が必須（または `amplicon_out` / その `gene_assignments` から解決できること）。  
Assemble consensus からの editcall は提供しない（パラメータ切り替えなし）。

### 7.3 `pam_list` 契約（決定・**実装済み**）

CSV（ヘッダーなし）:


| 列   | 内容                                                                                                |
| --- | ------------------------------------------------------------------------------------------------- |
| 1   | gene ID（`primer_list` / `amplicon.fa` の gene と一致）                                                 |
| 2   | **chromosome / seqname** — `genome_fn` の FASTA header と **文字列一致**（`chr` や `%02d` をツール側で付けない）      |
| 3   | PAM 開始位置（1-based）                                                                                 |
| 4   | （任意だが複数行時は**必須**）guide ID                                                                         |
| 5   | （任意）strand `+` / `-`。`+` → cut = pam−3、`-` → cut = pam+3（**genome 上**）。欠落時は cut = pam 開始（オフセット 0） |


- 列2が genome に無ければ **エラーで解析停止**（`chr%02d` 強制は **廃止済み**）。
- 破壊的変更: 既存 `agr8_pam_list.csv` は `chr01` 等へ更新済み（`inst/extdata/`）。
- 編集窓・案 A′ の `c_i` はいずれも **cut 中心**（J1c）。amplicon↔genome は `matchPattern`（prep と同様）。
- 写像実装: `.parse_pam_list` / `.map_amplicon_to_genome` / `.build_editcall_metadata`（`R/editcall.R`）。



### 7.4 複数ガイド（決定: 案 A + 案 A′・**実装済み**）


| ケース                                 | 運用                                                               | ツール                    |
| ----------------------------------- | ---------------------------------------------------------------- | ---------------------- |
| **1** ガイドごとに別 primer / 別 amplicon   | gene ID を `SD1_1` / `SD1_2` のように分ける                              | 従来どおり gene 単位          |
| **2** 同一 amplicon・同一 primer に複数 PAM | gene ID は共有（例: `SD1`）。`pam_list` に同 gene・異 position・異 guide の複数行 | **案 A** + **案 A′**（下記） |


**案 A（ケース2・ガイド別窓・採用・実装済み）:**

1. `.build_editcall_metadata` は同 gene の PAM 行をすべて保持（`pos_idx[1]` 先頭固定を廃止済み）。
2. 各 PAM について genome→amplicon→insert 写像 → `aln_pos` 1行。集計キーは `target_gene = paste(gene, guide, sep="_")`（単一行で guide 欠落時は `target_gene = gene`）。
3. primers / `amplicon.fa` は **gene 名で共有**。
4. per-read: primer-trim は1回 → insert↔ref_insert の **global align は1回** → 紐づく PAM すべてで §15.7 窓抽出 → レコードをガイド数ぶん出力。キャッシュキーは `target_gene` **を含む**。
5. 検証: 同 gene で guide 重複 → 停止。guide 無しで同 gene が2行以上 → 停止。pam の gene が primer/amplicon に無い → 停止。

**案 A′（ケース2・ジョイント両 cut 検証・採用・実装済み）:**

案 A だけでは、両 cut による **介在大欠失（excision）** で片窓だけ discard／別アレル化することがある。同一アラインメント上にジョイント層を載せ、両 cut の同時発生を検証する。案 A（ガイド効率）と **相補**（置き換えではない）。

アルゴリズム（gene に PAM ≥ 2 のとき）:

1. cut 座標を insert 上の ungapped 座標 `c1 < c2 < …` に載せる（ソート: insert cut、同点は `guide_id`）。
2. **ガイド別窓**（案 A）はそのまま実行。
3. **隣接ガイド対**（ソート順の隣接ペアのみ。3ガイド以上の全組合せは後続可）ごとに:
  - 外側錨: `c_i − check_window` … `c_j + check_window`（§15.7 と同型の外側 expand）。
  - 介在区間の参照欠失長 `del_span` と期待長 `expected = c_j - c_i` を比較（許容: `|del_span - expected| ≤ excision_tol_bp`、既定 20）。
  - `event_class`:

    | class                   | 条件                                          |
    | ----------------------- | ------------------------------------------- |
    | `wt`                    | 両局所窓が intact                                |
    | `g_i_only` / `g_j_only` | 片方のみ非 WT                                    |
    | `both_local`            | 両局所非 WTだが介在はおおむね残存（または外側錨のみの部分ケース）          |
    | `both_cut_excision`     | 外側錨 OK かつ `del_span ≈ expected`（両 cut 介在欠失） |

  - cut 間が `min_span_bp`（既定 **30**）未満なら excision 判定は出さず `both_local` / `wt` 等へ倒す。
4. 外側錨も取れないリードは joint も discard（無理に呼ばない）。小窓が片側失敗しても、外側錨＋介在欠失が通れば `both_cut_excision` は数えられる。
5. excision 時: 案 A 行にも接合アレルを載せる（option B: gapped junction、`intact=FALSE`、`allele_source="excision"`）。ガイド間の二重計上は許容。

出力（実装済み）:


| ファイル                         | 内容                                                                                    |
| ---------------------------- | ------------------------------------------------------------------------------------- |
| `editcall_all.csv` 等（既存）     | 案 A: ガイド別 `target_gene` 集計（主）                                                         |
| `editcall_joint.csv`         | リード単位: `gene_id`, `guide_i`, `guide_j`, `event_class`, `del_span`, `expected_span`, … |
| `editcall_joint_summary.csv` | sample×gene×ペアの `excision_rate` 等（分母＝当該 gene で案 A 行が出たリード数）                           |


既知の限界（許容）: 接合アレルの細部分類、非隣接ペアの同時欠失、極めて近い二重 cut の局所／excision 境界はヒューリスティック依存。ドキュメント化する。

### 7.5 evalMiao / editViewer との整合


| 関数                    | 期待するディレクトリ                                                                                                                                                   |
| --------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `evalMiao(out_dir)`   | `demultiplex/` 必須。gene 集計は `amplicon_assign/gene_assignments.tsv` を優先し、無ければ `amplicon/gene_assignments.tsv`。editcall 集計は `editcall/`。案 A′ 後は任意で joint 集計を参照可 |
| `editViewer(out_dir)` | `editcall/editcall_summary.csv`（経路 B のみでも動作）。案 A 後は `target_gene`（ガイド単位）行が増える。ジョイント excision の plate 表示は後続でも可                                                |


Editcall 専用ランでは `amplicon/` が無くてもレポートが成立すること。

---



## 8. エントリポイント



### 8.1 推奨（2 経路）

`miaoEditcall()` **は使わない。** gene 割当後に目的に応じて分岐する。

```r
library(miaoseq)

fastq <- "/path/to/reads.fq"
out_dir <- "/path/to/run"
demult_dir <- file.path(out_dir, "demultiplex")
assign_dir <- file.path(out_dir, "amplicon_assign")
amplicon_dir <- file.path(out_dir, "amplicon")
editcall_dir <- file.path(out_dir, "editcall")

index_list <- system.file("extdata", "index_list.csv", package = "miaoseq")
primer_list <- system.file("extdata", "amplicon_primers.csv", package = "miaoseq")
pam_list <- system.file("extdata", "agr8_pam_list.csv", package = "miaoseq")
genome_fn <- "/path/to/genome.fa"
amplicon_fn <- prepAmpliconDB(primer_list, genome_fn, out_dir)

dem <- doDemultiplex(
  fastq = fastq,
  demult_dir = demult_dir,
  index_list = index_list,
  split_reads = TRUE          # editcall / assign では by_sample 推奨
)

ga <- doAssignGenes(
  assignments = file.path(demult_dir, "assignments.tsv"),
  out_dir = assign_dir,
  sample_fastq_dir = file.path(demult_dir, "by_sample"),
  primer_list = primer_list,
  amplicon_fn = amplicon_fn,
  overwrite = TRUE
)

# ---- 経路 B: Editcall（本線・Assemble 不要）----
edit <- doEditcall(
  gene_assign = ga,
  pam_list = pam_list,
  genome_fn = genome_fn,
  amplicon_fn = amplicon_fn,
  primer_list = primer_list,
  editcall_dir = editcall_dir,
  sample_fastq_dir = file.path(demult_dir, "by_sample"),
  check_window = 10L,
  anchor_bp = 5L,
  max_expand = 50L,
  min_count = 5L
)
evalMiao(out_dir)
editViewer(out_dir, sample_list = "/path/to/sample_list.csv")

# ---- 経路 A: Assemble（代表配列・QC、任意）----
amp <- doAssembleAmplicons(
  gene_assign = ga,
  out_dir = amplicon_dir,
  primer_list = primer_list,
  cluster_backend = "vsearch",
  min_cluster_purity = 0.8,
  overwrite = TRUE
)
```



### 8.2 `prepAmpliconDB`

```r
prepAmpliconDB(
  primer_list,
  genome_fn = NULL,
  out_dir,
  amplicon_fasta = NULL,
  expected_length = NULL   # ゲノム導出時は必須; span > 2*expected_length で停止
)
```

- **ゲノム導出:** F/R primer を `matchPattern`（完全一致・両鎖）で探索し、gene ごとに一意の F/R ペア区間を切り出す。同一 gene の複数ヒット、または span > `2 * expected_length` があれば **停止**。
- **ユーザー FASTA:** `amplicon_fasta` を渡すとゲノム探索を省略。配列名は primer の gene ID と一致必須。任意で `expected_length` による `2×` 上限チェック。
- BLAST 不要。



### 8.3 Editcall genotype ラベル（indel）と editViewer

アライン済み `read_seq` / `ref_seq` について:

- 欠失 bp = `read_seq` 内の `-` 個数
- 挿入 bp = `ref_seq` 内の `-` 個数
- ラベル: `delN` / `insN` / `indelD-I`（例: `indel5-3`）

`editViewer` はラベルを構造化パースする（`ref` / `sub` / `delN` / `insN` / `indelD-I`）。

- `sub`（SNP）は **frameshift 扱いにしない**
- in-frame 色分けは **net indel 長 mod 3**（CDS 位相モデルではない）
- 複合 indel は `|n_ins - n_del| %% 3`

---



## 9. パッケージ境界の再定義



### miaoseq がやる / やらない


| やる                    | やらない                  |
| --------------------- | --------------------- |
| demux（edlib）          | Dorado basecall       |
| サンプル別 FASTQ 分割        | pod5 管理               |
| サンプル別クラスタ / コンセンサス    | minimap2 index 構築の必須化 |
| editcall / report（任意） | 実験室 LIMS              |




### README に明記する前提

```bash
# 例: ユーザーが先に basecall
dorado duplex sup pod5/ --min-qscore 10 > basecall.bam
samtools fastq basecall.bam > reads.fq
# 任意: 長さフィルタ等

# その後 R
doDemultiplex(fastq = "reads.fq", ...)
```

あわせて README に以下を書く（一部反映済み）:

- 組合せ index では誤 dual が合法ウェルに入りうること（§5.1）
- `fraction_sample` は sample×gene バケツ比率であること（§6.3）
- `pam_list` 染色体は genome seqname 一致・複数ガイドは §7.3–7.4（案 A + 案 A′）

---



## 10. ファイル構成（現状）

```
R/
  demultiplex.R          # doDemultiplex, splitDemultiplexReads
  amplicon_resolve.R     # doAssignGenes, doAssembleAmplicons
  prep_amplicon_db.R     # prepAmpliconDB（vmatchPattern）
  editcall.R             # doEditcall
  report.R               # evalMiao, editViewer

src/
  amplicon_core.cpp      # primer 割当 / trim / cluster / consensus
```

NAMESPACE export:

- `doDemultiplex`, `splitDemultiplexReads`
- `doAssignGenes`, `doAssembleAmplicons`, `doEditcall`
- `prepAmpliconDB`, `evalMiao`, `editViewer`

Phase B で削除: `pipeline.R`, `amplicon_assign.R`, `basecall.R`, `blast_utils.R`

---



## 11. データ契約サマリ

```r
# Step 1 → Step 2
DemultiplexResult <- list(
  assignments = data.frame(read_id, sample_id, index_pair_id, ...),
  summary     = data.frame(sample_id, n_reads, ...),
  unassigned  = data.frame(read_id, reason, ...),
  demult_dir  = character(1)
)

# Step 1b（任意）
# demultiplex/by_sample/{sample_id}.fq.gz

# Step 2 → Step 3
AmpliconResolveResult <- list(
  samples   = character(),
  out_dir   = character(1),
  table     = data.frame(),          # cluster 縦結合（seq 列含む）
  gene_assignments = data.frame()
)

# Step 3 出力
# editcall/editcall_summary.csv  （editViewer 互換形式を維持）
```

---



## 12. 実装フェーズ（パイプライン全体）



### Phase A — Demux 刷新 ✅ 完了

- [x] [demux_revise.md](demux_revise.md) に沿って `doDemultiplex` 実装
- [x] `splitDemultiplexReads` 実装 + demux からの呼び出し
- [x] 旧 BLAST demux 除去（edlib 本線）
- [x] README 新フロー反映



### Phase B — basecall / レガシー / BLAST 撤去 ✅ 完了

- [x] `R/basecall.R` 削除（`doBasecall`, `makeMMI`）
- [x] `R/pipeline.R` 削除（`miaoEditcall`）
- [x] `R/amplicon_assign.R` 削除（`doAlign`）
- [x] `R/blast_utils.R` 削除
- [x] `prepAmpliconDB` を `Biostrings::matchPattern` 実装に置換（`blast_path` 引数廃止）
- [x] `doEditcall` からレガシー `demult_out` + `align_out` 経路を削除
- [x] NAMESPACE / man / README / `sample_script.R` 更新



### Phase C — Amplicon resolve ✅ 完了

- [x] 入出力契約（[ampliconresolve_plan.md](ampliconresolve_plan.md)）
- [x] C0: 配線・完全一致クラスタ
- [x] C1: primer edlib + greedy cluster + majority consensus
- [x] C2: `doAlign` deprecated、フルラン確認
- [x] C3: quality / spoa / vsearch / strict_end_trim（ablation 済み）



### Phase D — Editcall ✅ 完了（初期）

- [x] `doEditcall(amplicon_out, ...)` 実装
- [x] vsearch+spoa 経路でフル検証（386 samples）
- [x] `evalMiao` / `editViewer` 更新



### Phase I — 2 経路分離 + 純度フィルタ ✅ 完了

詳細は §0.1・§7・§15.3.A/E・§15.4。

- [x] `doEditcall(gene_assign = …)` を経路 B 本線 API にする（Assemble 非依存・reads 専用）
- [x] `evalMiao` が `amplicon_assign/gene_assignments.tsv` を優先読込
- [x] `doAssembleAmplicons` に `min_cluster_purity` / `cluster_purity` 列（§15.3.E）
- [x] `allele_source` / consensus editcall を廃止（切り替えなし）
- [x] README / `sample_script.R` / `editcall_adapter.md` を 2 経路に更新
- [ ] 回帰: gene_assign → editcall → evalMiao / editViewer のスモーク（データ依存・手動）



### Phase I′ — 適応編集窓（§15.7）✅ 完了

- [x] `doEditcall(anchor_bp, max_expand)` API
- [x] `.extract_adaptive_edit_window`（両端アンカー／外側拡張／棄却）
- [x] `read_seq` + `ref_seq` 同一窓切り出し、`intact = identical(read_seq, ref_seq)`
- [x] genotype ins サイズを同窓 `ref_seq` 長で計算
- [x] 長さ安全網 `2*check_window + 2*max_expand + 20`
- [x] 単体テスト `cursor_dev/tool_test/test_adaptive_edit_window.R`
- [ ] 実データ full 回帰（大 indel／アライン失敗件の分布確認）



### Phase H — 最終クラスタ未寄与リードの再検証 ✅ MVP 完了

詳細は **§16**。

- [x] U1 定義と assemble 時 `unassigned_to_cluster.tsv` 出力
- [x] `doReassessAssemblies()`: Q1 + Q2、`amplicon_reassess/`（本線非破壊）
- [x] バックエンド比較: `edlib` / `vsearch` / `blastn` / `mmseqs` + `summary_compare.tsv`
- [ ] Q3（未寄与サブクラスタ）/ Q4（閾値スイープ）
- [ ] 代表実データランでの回帰



### Phase J — 監査中程度対応（§17–§18）✅ 完了

詳細は **§17**（方針）と **§7.3–7.4**（実装契約）。回答記録は **§18**。

- [x] 方針文書化（§5.1 / §6.3 / §7.3–7.4 / §17）
- [x] README: 組合せ demux / `fraction_sample` / PAM 注意の明記
- [x] 案 A′（ジョイント両 cut）を §7.4 / §17 に採用記載
- [x] 関連 md 同期（demux_revise / editcall_adapter / ampliconresolve_plan）
- [x] **§18 未決事項への回答記入**（J1–J11 / J1b–J11b / J1c）
- [x] `require_unique_pair` 削除（`R/demultiplex.R` / man）
- [x] PAM: genome seqname 一致 + `chr%02d` 廃止 + extdata 更新
- [x] 複数ガイド案 A 実装（ガイド別窓・`target_gene` キャッシュ）
- [x] 複数ガイド案 A′ 実装（`editcall_joint*.csv`、隣接ペア excision、`min_span_bp`/`excision_tol_bp`）
- [x] vsearch ↔ 絶対 edit: README 文書一行のみ（J11b 案2）。技術揃えは **別 Phase・todo**
- [x] 単体テスト `tests/test-editcall-phase-j.R`
- [ ] 実データフルラン回帰（任意・データ依存）
- [ ] evalMiao / editViewer の joint 集計・plate（意図的に Phase J 外・todo）

---



## 13. 成功基準（全体）


| #   | 基準                                                         | 状態                                         |
| --- | ---------------------------------------------------------- | ------------------------------------------ |
| 1   | R から Dorado を呼ばず demux → resolve まで通る                      | ✓                                          |
| 2   | demux 主出力が assignments、split は任意                           | ✓                                          |
| 3   | サンプル単位で consensus / cluster を出せる                           | ✓                                          |
| 4   | 旧「全リード BLAST → join」依存を解消                                  | ✓                                          |
| 5   | README 上の前提（ユーザー用意 FASTQ）が一貫                               | ✓                                          |
| 6   | editcall が `gene_assign` を消費できる（Assemble 非依存）              | ✓（Phase I）                                 |
| 7   | vsearch + spoa で Assemble E2E が通る                          | ✓（full 検証）                                 |
| 8   | BLAST 依存なし（prepAmpliconDB 含む）                              | ✓（Phase B）                                 |
| 9   | `gene_assign` → `doEditcall` → evalMiao / editViewer       | ✓（Phase I）                                 |
| 10  | Assemble にクラスタ純度フィルタ                                       | ✓（Phase I / §15.3.E）                       |
| 11  | Step 2c で未寄与リード再検証レポートを出せる                                 | ✓（§16 MVP: `reassignable_frac` / merge 成分） |
| 12  | 適応編集窓（§15.7）が経路 B 既定で動作する                                  | ✓（単体テスト済・実データ回帰は任意）                        |
| 13  | §17: PAM seqname 一致・複数ガイド案 A + A′・`require_unique_pair` 削除 | ✓（Phase J）                                 |


---



## 14. 決定済み / 残タスク



### 決定済み


| 項目                                        | 決定                                                                      |
| ----------------------------------------- | ----------------------------------------------------------------------- |
| クラスタリング                                   | C1 既定: greedy edlib。**推奨:** vsearch                                     |
| コンセンサス                                    | **（歴史・Phase C/D 当時）spoa 固定**。現行は **abpoa**（§19） |
| マルチ gene 分離                               | primer edlib（`primer_list`）                                             |
| 処理方針                                      | ref-aware（C）                                                            |
| `prepAmpliconDB`                          | **matchPattern**（両鎖完全一致）。multi-hit / 長区間で停止。`amplicon_fasta` 可。BLAST 廃止 |
| `miaoEditcall` / `doAlign` / `doBasecall` | **Phase B で削除済み**                                                       |
| Assemble / Editcall                       | **2 経路分岐**（Phase I）。Editcall は Assemble 非依存                             |
| Editcall `read_seq`                       | **適応編集窓**（§15.7）。`anchor_bp=5`, `max_expand=50`                         |
| 組合せ demux のリスク（§5.1 / §17）                | **利用者責任**。1:1 または高 Hamming を推奨。ツール側で topology 救済しない                     |
| `require_unique_pair`                     | **削除**（no-op だった）                                                       |
| `fraction_sample`                         | **実装据え置き**（gene バケツ比率）。ドキュメントで明記                                        |
| PAM 染色体（§7.3）                             | genome seqname と一致必須。不一致は停止。`chr%02d` 付与**廃止済み**（Phase J）               |
| 複数ガイド（§7.4）                               | ケース1=gene 分割。ケース2=**案 A** + **案 A′・実装済み**（Phase J）                      |
| vsearch identity ↔ 絶対 edit（§17）           | Phase J では **README 注意のみ**。バックエンド揃え／換算式は **要検討（別 Phase）**               |




### 残タスク


| 項目                                                                                       | 優先度                             |
| ---------------------------------------------------------------------------------------- | ------------------------------- |
| `evalMiao` / `editViewer` 更新（旧 out_dir 契約）                                               | 完了                              |
| ~~**editcall 偽陽性対策（実装）**~~ per-read `doEditcall`                                         | **完了**                          |
| ~~**Step 2 graph assembly**~~ `assembly_backend = "overlap_graph"`                       | **完了**                          |
| ~~**Phase I: 2 経路分離 +**~~ `gene_assign`~~**→**~~`doEditcall` ~~**+ 純度フィルタ（§15.3.A/E）**~~ | **完了**                          |
| resolve サンプル並列                                                                           | 低                               |
| vsearch `id` と内部絶対 edit の揃え方（§17 #3）                                                     | **要検討**（別 Phase；文書一行は済）         |
| ~~異常~~ `read_seq` ~~/ アライン失敗棄却~~（§15.3.F + §15.7）                                        | **完了**（`max_expand` 棄却 + 長さ安全網） |
| ~~**適応編集窓**（§15.7）~~                                                                     | **完了**                          |
| ~~**最終クラスタ未寄与リードの再検証**（§16 / Phase H）MVP~~                                               | **完了**（edlib Q1/Q2）             |
| Phase H 拡張（Q3/Q4・実データ回帰）                                                                 | 任意                              |
| 実データ回帰（適応窓・gene_assign → editcall → evalMiao・Phase J）                                    | 手動・データ依存                        |
| `require_unique_pair` ~~削除~~                                                             | **完了**（Phase J）                 |
| ~~**PAM: seqname 一致 +**~~ `chr%02d` ~~**廃止** +~~ `agr8_pam_list.csv`                     | **完了**（Phase J）                 |
| ~~**複数ガイド案 A / A′**（§7.4）~~                                                              | **完了**（Phase J）                 |
| ~~**§18 未決事項への回答**~~                                                                     | **完了**（回答→実装済）                  |
| evalMiao / editViewer の joint 集計・plate 表示                                                | 後続（J6: Phase J は joint CSV のみ）  |
| `fraction_sample` の誤解防止文言を README 外ドキュメントにも揃える                                           | **済**（2026-07-14）               |


---



## 15. 類似度クラスタリングと editcall 偽陽性



### 15.1 問題の概要

旧パイプライン（`doAlign` → `doEditcall`）と新パイプライン（`doAssignGenes` → `doAssembleAmplicons` → `doEditcall`）を同一データで比較すると、新フローの `editcall_all.csv` では **シーケンシングエラーに由来しそうな** `read_seq` **パターンが、旧より多くのリードで支持されている** ように見えることがある。

**比較例（同一 FASTQ 系統）:**


| 出力  | パス                                                   | 行数（allele 行） | 集計単位                                         |
| --- | ---------------------------------------------------- | ------------ | -------------------------------------------- |
| 旧   | `…/20250918_pipeline_test/editcall/editcall_all.csv` | ~112,700     | リード単位の `read_seq` 完全一致                       |
| 新   | `…/run_report_test/editcall/editcall_all.csv`        | ~4,300       | クラスタ consensus から抽出した `read_seq` × `n_reads` |


行数は新の方が少ない（クラスタ集約のため）が、**同一 sample × gene ×** `read_seq` **の** `count` **を比較すると、新フローで大幅に増加したパターンが多数ある**（例: `miaoBC0289` / SD1 / `GCGTGCGCCACGCACGGGTTC` が旧 3 → 新 2,281）。また旧に存在しなかった非 intact パターンが `count ≥ 20` で 400 件以上、新フローにのみ現れる。

### 15.2 原因の整理

旧 `doEditcall` の集計は、実質的に **PAM / cut site 周辺** `check_window`**（既定 10 bp）の配列をリードごとに抽出し、文字列完全一致で** `table(read_seq)` していた。頻度の低いエラーパターンは自然に落ちる。

新フローでは Step 2 が先に **フル長アンプリコン insert 上で類似度クラスタリング**（greedy edlib または `vsearch --cluster_fast`、`max_cluster_edit = 12` 相当）と **コンセンサス**（推奨: spoa）を行い、Step 3 が **クラスタ代表配列から** `read_seq` **を再抽出して** `n_reads` **を重み付け** する。

この二段構えにより、次の「集約バイアス」が起きうる:

1. **異なるエラーパターンの合流** — フル長では数 bp の差しかないリードが同一クラスタに入り、コンセンサスが「中間的な」配列を生成する。旧では別 allele だったエラーが、新では 1 パターンに `n_reads` が積み上がる。
2. **真の編集とエラーの合流** — 新規編集配列と、参照に 1–2 bp 近いエラーリードがクラスタを共有すると、編集らしい consensus が形成され、変種の一つに見える。
3. **spoa による丸め込み** — クラスタ内のヘテロジェニティを多数決より強く平滑化し、エラー由来の配列が「支持されている変種」に見えることがある。
4. **（副次）アライン失敗** — consensus／リードと参照 insert の global alignment が崩れた場合、`read_seq` が `check_window` 長（~21 bp）を大きく超える異常行が出る（新フロー full 検証で 68 行が 100 bp 超）。これも偽陽性・解釈不能 allele の一因。
5. **（補足）固定窓と大 indel** — 参照上の固定 `± check_window`（~21 bp）内に大きな欠失があると、クエリ側の切り出しに **参照と相同な塩基がほとんど残らない**（ほぼ gap、または端が崩れた短い列）。大きな挿入では変異本体は入るが左右の「錨」が薄いことがある。固定長 `read_seq` だけではgenotype／可視化が不安定になりうる（対策は §15.7）。

要するに、**Step 2 の目的（アンプリコン代表配列の復元）と Step 3 の目的（編集窓の変異頻度推定）は、同じクラスタリング設定を共有すべきではない**。

### 15.3 対策の方向性（優先度順）



#### A. リード単位 edit-window 集計の復帰（**確定・経路 B 本線**）

gene 割当（`doAssignGenes`）までは共通とし、**Editcall は Assemble を経由せず**リード単位で `read_seq` を抽出・完全一致集計する（2026-07-14 確定）。

```
doDemultiplex → doAssignGenes → doEditcall
                                  → evalMiao / editViewer
```

- **長所:** 旧パイプラインの「低頻度エラーは自然に落ちる」性質を再現。§15.1 の集約バイアスを構造的に回避。
- **短所:** リード数に比例する処理（BLAST 全リード align よりは軽い）。
- **I/O:** `doEditcall` の唯一の入力契約は `gene_assign`（`GeneAssignResult` / `amplicon_assign/` / `gene_assignments.tsv`）+ FASTQ。`allele_source` 切替は持たない。
- **レポート:** `evalMiao` は `amplicon_assign/gene_assignments.tsv` を gene 集計に用いる（`amplicon/` が無い Editcall 専用ランでも成立）。



#### B. 二段階集計（ハイブリッド）

Step 2 のクラスタリングはアンプリコン代表・QC 用に維持し、editcall では **クラスタ内リードに戻って edit-window だけ完全一致集計**する。

- クラスタ `k` の `n_reads` をそのまま consensus に載せるのではなく、メンバーリードごとの `read_seq` を数えてから合算。
- クラスタ内に複数の edit-window ハプロタイプがある場合は **サブクラスタとして分裂**するか、最大ハプロタイプのみ採用（`cluster_purity` 閾値）。



#### C. パラメータ緩和（即効だが根本解ではない）

editcall 精度が主目的のランでは、少なくとも以下を試す:


| パラメータ                             | 現行推奨        | editcall 重視時の案                                           |
| --------------------------------- | ----------- | -------------------------------------------------------- |
| `cluster_backend`                 | `vsearch`   | `internal`（greedy）または **完全一致のみ**（`max_cluster_edit = 0`） |
| `max_cluster_edit`                | `12`        | `0`（完全一致）〜 `4`                                           |
| Consensus                         | **spoa 固定** | —（majority/quality 廃止）                                   |
| `min_cluster_reads` / `min_count` | `5`         | 引き上げ（例: `10`–`20`）                                       |


完全一致クラスタ（Phase C0 相当）は過剰分裂するが、**editcall 入力としては旧** `table(read_seq)` **に最も近い**。

#### D. edit-window 限定クラスタリング

フル長ではなく **primer-trimmed insert の cut site 周辺のみ**をクラスタリング距離の対象にする。アンプリコン末端の ONT エラーがクラスタ合流を起こさない。

- 実装案: `trim_amplicon_insert` 後に `check_window` 区間だけを `greedy_cluster_cpp` / vsearch に渡す。
- フル長 consensus は別途必要なら残す。



#### E. クラスタ純度フィルタ（post-cluster QC）— **Phase I 実装済み**

各クラスタについて、メンバー配列の **最頻ハプロタイプ占有率**（purity）を計算し、閾値未満のクラスタを棄却または代表を最頻列に置き換える。


| 適用先                            | 列 / 挙動                                                                                                                                 |
| ------------------------------ | -------------------------------------------------------------------------------------------------------------------------------------- |
| **経路 A** `doAssembleAmplicons` | `cluster_counts.tsv` に `cluster_purity`（insert 完全一致の最頻占有率）。任意で PAM 指定時は `edit_window_purity` も算出。`min_cluster_purity`（既定 0.8）未満は出力から除外 |
| **経路 B** `doEditcall`          | クラスタを経由しないため **純度フィルタ対象外**（リード単位集計が純度問題を回避）                                                                                            |


- spoa consensus を editcall に渡す経路は廃止。低純度クラスタの影響は **Assemble 側の純度フィルタ**で抑える。



#### F. 下流フィルタの強化

既存の `editcall_filtered.csv` 生成（`count > top/2`）に加え:

- 参照配列からの編集距離が 1 bp の sub のみ、かつ全体 read 深度に対する比率が低い → エラー候補として除外（**未実装・任意**）
- homopolymer 近傍の単一置換を ONT エラーモデルで減点（**未実装・任意**）
- **実装済み:** §15.7 の `max_expand` アンカー失敗棄却が主防御。加えて長さ安全網 `nchar(read_seq) > 2*check_window + 2*max_expand + 20`（既定 140 bp）で異常行を落とす



#### G. graph assembly（Step 2 向け — hifiasm 的アプローチの評価）

**対象:** `doAssembleAmplicons` のクラスタ + コンセンサス代替（Step 3 の editcall とは別レイヤ）。

**結論:** **hifiasm を sample×gene ごとに直呼するのは非推奨。** 入力は数十〜数千リード × ~200–400 bp の PCR アンプリコンであり、ゲノムアセンブラの前提（高深度・長リード・単一ジョブ）と合わない。384 samples × 8 genes ≈ 3,000 回の起動コストも過大。

一方、hifiasm が行う **overlap graph → バブル分離 → 低カバレッジ枝刈り** の考え方は Step 2 に有効。miaoseq 規模では次が現実的:


| 案                     | 内容                                                                                         | 採用              |
| --------------------- | ------------------------------------------------------------------------------------------ | --------------- |
| **overlap_graph**（推奨） | edlib 距離でリード間グラフを構築。高支持の **完全一致配列を path** とし、類似低支持 path を枝刈り。代表配列は spoa 平滑化ではなく **最頻実リード** | **Phase G1 実装** |
| abpoa マルチコンセンサス       | POA graph から複数 path 列挙                                                                     | G3（必要時）         |
| spoa バブル列挙            | 既存 spoa 拡張                                                                                 | G3              |
| hifiasm 直呼び           | sample×gene FASTQ → contig                                                                 | **非採用**         |


```r
doAssembleAmplicons(
  ...,
  assembly_backend = c("cluster", "overlap_graph"),  # 既定 "cluster"
  overlap_min_identity = 0.90
)
```

`assembly_backend = "overlap_graph"` 時は cluster + consensus をスキップし、`cluster_counts.tsv` / `consensus.fasta` 契約は維持（各 path = 1 cluster、配列は modal exact read）。

**Step 2 と Step 3 の責務分離（重要）:**


| 層               | 役割               | 推奨                                      |
| --------------- | ---------------- | --------------------------------------- |
| Step 2 graph    | フル長アンプリコン代表配列・QC | `overlap_graph` で偽合流を抑制                 |
| Step 3 editcall | cut site 周辺の変異頻度 | `doEditcall`（per-read・§15.7 適応窓・完全一致集計） |


Step 2 graph 化だけでは editcall 偽陽性は完全には解消しない。§15A と併用する。

#### H. 適応編集窓（大 indel）— **§15.7 実装済み**

経路 B の `read_seq` 抽出を、左右に連続 `anchor_bp`（既定 5）の参照一致が付くまで必要に応じて外側拡張する。両端とも初期窓で満たせば従来どおり。`max_expand` 超過は棄却。  
**固定窓で先にクラスタしてから伸ばすのではなく、リードごとに窓を決めてから完全一致集計する。** 詳細・intact 同規則・誤解防止は **§15.7**。

### 15.4 確定実装方針（Phase I）


| 優先  | 施策                                                                  | 状態                |
| --- | ------------------------------------------------------------------- | ----------------- |
| 1   | **A:** 経路 B = `doAssignGenes` → `doEditcall`（Assemble 非依存・reads 専用） | **確定・実装**         |
| 2   | consensus / `allele_source` 切替は **廃止**                              | **確定・実装**         |
| 3   | **E:** Assemble に純度フィルタ                                             | **実装**            |
| 4   | **F:** 異常 `read_seq` 棄却（長さ安全網）                                      | **実装**（§15.7 と併用） |
| 5   | **§15.7 / H:** 適応編集窓（アンカー 5 bp / `max_expand`）                      | **実装**            |
| 6   | C / D はドキュメント / preset                                              | 任意                |


API（経路 B 本線）:

```r
ga <- doAssignGenes(...)

edit <- doEditcall(
  gene_assign = ga,                 # 主入力（Assemble 不要）
  pam_list = pam_list,
  genome_fn = genome_fn,
  amplicon_fn = amplicon_fn,
  primer_list = primer_list,
  editcall_dir = file.path(out_dir, "editcall"),
  fastq = fastq,                    # または sample_fastq_dir
  check_window = 10L,
  anchor_bp = 5L,
  max_expand = 50L,
  min_count = 5L
)

evalMiao(out_dir)
editViewer(out_dir, sample_list = sample_list)
```

API（経路 A + 純度）:

```r
amp <- doAssembleAmplicons(
  gene_assign = ga,
  out_dir = file.path(out_dir, "amplicon"),
  primer_list = primer_list,
  min_cluster_purity = 0.8,         # §15.3.E
  ...
)
```

互換: `doEditcall(amplicon_out = amp, ...)` は、`amp` 内の `gene_assignments`（または同パス上の TSV）があれば gene 入力のフォールバックとして動作する（移行期間）。

### 15.5 検証基準

旧 `20250918_pipeline_test` と経路 B を同一 FASTQ で再比較し、以下を満たすこと:

1. sample × gene × `read_seq` の `count` 相関（Pearson / Spearman）が旧パイプライン（リード単位）と整合する
2. 旧に存在しない高支持（`count ≥ 20`）非 intact パターンの件数が、consensus 経路時代より減少していること（経路 B 本線化で構造的に対処済み）
3. `max_expand` 棄却後、長さ安全網を超える `read_seq` 行が 0（または警告ログのみ）
4. 既知の編集サンプルで真の編集 allele の検出感度が低下しないこと
5. （§15.7）大欠失で適応窓 `read_seq` が左右に連続 `anchor_bp` の参照一致を含むこと — **単体テスト済**（`test_adaptive_edit_window.R`）
6. （§15.7）非整列クエリは `max_expand` で棄却されること — **単体テスト済**



### 15.6 関連コード・ドキュメント


| 箇所                                                      | 内容                                                                           |
| ------------------------------------------------------- | ---------------------------------------------------------------------------- |
| `R/amplicon_resolve.R`                                  | `doAssignGenes` + `doAssembleAmplicons`、`assembly_backend = "overlap_graph"` |
| `R/editcall.R`                                          | 経路 B: per-read 抽出・集計 + §15.7（`.extract_adaptive_edit_window` 等）              |
| `cursor_dev/tool_test/test_adaptive_edit_window.R`      | §15.7 単体テスト（完全一致／大欠失／棄却）                                                     |
| [ampliconresolve_plan.md §6.4](ampliconresolve_plan.md) | Phase C0 完全一致クラスタの位置づけ                                                       |
| [editcall_adapter.md](editcall_adapter.md)              | editcall 接続メモ（reads 専用・適応窓）                                                  |




### 15.7 適応編集窓（大 indel 向けアンカー拡張）— **実装済み**

作成: 2026-07-14  
状態: **実装済み**（`R/editcall.R`；単体テスト `cursor_dev/tool_test/test_adaptive_edit_window.R`）  
対象: **経路 B** `doEditcall` の `read_seq` 抽出（Assemble／フル長類似度クラスタとは独立）

#### 15.7.1 動機

固定 `check_window`（既定 PAM／cut ±10 bp ≈ 21 bp）は小変異には十分だが、大きな欠失ではクエリ側の切り出しに参照相同塩基が残らず、大きな挿入では左右の錨が薄い（§15.2.5）。  
一方、gene 誤割当などで global alignment が崩れると窓切り出しが膨張し、超長い解釈不能 `read_seq` になる（§15.2.4）。  
→ **必要なときだけ外側に伸ばして左右アンカーを確保し、伸ばしきれないリードは棄却する。**

#### 15.7.2 用語の注意（誤解防止）

ここで言う「クラスタリング」は Step 2 の vsearch／edlib **類似度クラスタではなく**、抽出した `read_seq` 文字列の **完全一致グループ化**（従来の `table(read_seq)`）を指す。

**処理順（実装どおり）:**

1. 各リードについて適応窓で `read_seq` / `ref_seq` を決定する（または棄却）
2. **その後**、sample × gene × `read_seq`（+ `ref_seq` + `intact`）で完全一致集計する

固定 ±10 bp だけで先にクラスタし、代表列だけ後から伸ばす方式は **採用しない**。

#### 15.7.3 アルゴリズム（実装）

パラメータ（`doEditcall` 引数）:


| 名前             | 既定    | 意味                              |
| -------------- | ----- | ------------------------------- |
| `check_window` | `10L` | 初期編集窓: cut の左右 bp（参照・insert 座標） |
| `anchor_bp`    | `5L`  | 左右端に要求する参照との連続一致カラム数            |
| `max_expand`   | `50L` | 初期窓端から外側へ追加してよい最大 ungapped bp   |


実装関数:


| 関数                              | 役割                                               |
| ------------------------------- | ------------------------------------------------ |
| `.edit_seq_from_insert`         | primer-trimmed insert ↔ 参照 insert を global align |
| `.extract_adaptive_edit_window` | アンカー判定・外側拡張・棄却                                   |
| `.gapped_window_bounds`         | ungapped 参照座標 → ギャップ入りアライメント座標                   |
| `.alignment_end_has_anchor`     | 左／右端の連続 `anchor_bp` 一致判定                         |
| `.filter_abnormal_edit_records` | 長さ安全網                                            |


各リードについて:

1. メタデータ上の初期区間 `[start, end]`（PAM 窓・primer 長補正済み）でギャップ補正切り出しを開始。
2. **左端 / 右端:** アライン列の端から連続 `anchor_bp` カラムが「両非 gap かつ同一塩基」か。
3. **両端とも一致** → その区間の `query` / `ref` スライスを採用。
4. **不一致の側だけ** ungapped 座標を 1 bp ずつ外側へ広げ、再マップ（両側不一致なら両側）。
5. いずれかが `max_expand` **到達**でもアンカー不可 → **そのリードは棄却**（メッセージで件数を報告）。
6. 得られた `read_seq` で完全一致集計。

```
           ← expand left if needed ← │ initial ±check_window │ → expand right if needed →
… [===== 5bp match =====][ … variable indel / sub region … ][===== 5bp match =====] …
                              ↑ stop when each side has anchor_bp consecutive identity
                              ↑ abort read if either side hits max_expand without anchor
```



#### 15.7.4 intact・genotype との整合（実装どおり）

- アライン座標 `[left, right]` を決めたあと、**同一区間**から:
  - `read_seq = substr(query_aln, …)`
  - `ref_seq = substr(ref_aln, …)`
- `intact = identical(read_seq, ref_seq)`（ゲノム固定窓との文字列比較は使わない）
- genotype の ins サイズは `nchar(gsub("-", "", read_seq))` vs `nchar(gsub("-", "", ref_seq))`
- `intact_seq.fa` には従来どおりゲノム由来の既定 WT 窓（`± check_window`）を書き出す（Viewer / 参照用）
- `editcall_all.csv` / `editcall_filtered.csv` に `ref_seq` 列を出力



#### 15.7.5 「左右端が連続一致」の定義（実装）

- **一致カラム:** クエリ・参照がともに非 gap かつ同一文字。  
- **左端:** 窓の左から連続 `anchor_bp` カラムがすべて一致であること（途中で gap／不一致があれば失敗）。  
- 右端も対称。  
- 「窓のどこかに 5 bp 一致」ではない。  
- 拡張は必要最小限（初めて両端アンカーが揃った座標で停止）。



#### 15.7.6 §15.3.F・アライン失敗との関係


|     | F（長さカット）                                     | §15.7（アンカー）                     |
| --- | -------------------------------------------- | ------------------------------- |
| 主目的 | 異常に長い `read_seq` を事後棄却                       | 大 indel の錨確保＋伸ばしきれないリードを棄却      |
| 実装  | `2*check_window + 2*max_expand + 20`（既定 140） | `max_expand` 到達で `NULL` → リード棄却 |


両方が有効。アンカー失敗が主、長さカットは安全網。

#### 15.7.7 API（経路 B・実装済み）

```r
doEditcall(
  gene_assign = ga,
  ...,
  check_window = 10L,
  anchor_bp = 5L,
  max_expand = 50L,
  min_count = 5L
)
```



#### 15.7.8 非対象

- Step 2（Assemble）のフル長距離・spoa には載せない。  
- §15.3.D（edit-window 限定類似度クラスタ）とは別。§15.7 は **抽出文字列の定義**。  
- `max_expand` 到達後に部分一致のまま allele 採用しない。

---



## 16. Phase H: 最終クラスタ未寄与リードの再検証

作成: 2026-07-14  
状態: **MVP 実装済み**（`doReassessAssemblies()` / 複数バックエンド / Q1+Q2・診断専用）  
主目的: **召回の診断レポート**（成功目安: サンプルごとに `reassignable_frac` と近縁コンセンサスグループ数）。**自動マージは行わない**（採用判断はユーザー）。

### 16.1 背景と動機

Step 2b（`doAssembleAmplicons`）では、サンプルへ割当されたリードのうち、最終クラスタ（コンセンサス）に寄与する本数は典型的に少数に留まることがある。要因は複合的であり、少なくとも次を分離して議論する必要がある:

1. **gene / プライマー未割当**（`no_primer_hit` 等）
2. **プライマー割当済みだが、クラスタ閾値で落ちた**（`max_cluster_edit` 相当の類似度で主要クラスタに合流できず、かつ `min_cluster_reads` 未満のため破棄）
3. **真に別配列の稀少ハプロタイプ・キメラ・非特異産物**

現行パイプラインは (2)(3) を事後に区別する診断ステップを持たない。Step 2b 直後に **類似度に基づく再検証**を置くことで、(a) 既存コンセンサス同士がどの程度近いか、(b) 未寄与リードが既存コンセンサスへ再割当可能か、を定量し、閾値妥当性の根拠と「どのコンセンサスを採用するか」の判断材料とする。

本ステップは **元の** `amplicon/` **を書き換えない**。閾値内の近縁ペア／グループは **診断 TSV として提示するのみ**で、代表配列への自動統合や `merged_*.fasta` は出力しない（偽合流を避けるため、採用はユーザー判断）。

### 16.2 パイプライン上の位置

```
【経路 A】doAssignGenes → doAssembleAmplicons → 【NEW】再検証（Step 2c）
【経路 B】doAssignGenes → doEditcall(reads) → evalMiao / editViewer
          （Assemble 非経由。経路 A の再検証は経路 B と独立）
```


| 項目       | 内容                                                                                        |
| -------- | ----------------------------------------------------------------------------------------- |
| 関数名      | `doReassessAssemblies()`                                                                  |
| 入力       | `amplicon/`（読み取り専用）+ `gene_assign`（U1 配列復元）。assemble 時 `unassigned_to_cluster.tsv` があれば優先 |
| 出力ディレクトリ | `amplicon_reassess/`（診断 TSV のみ。マージ配列は出さない）。元の `amplicon/` は変更しない                          |
| バックエンド   | **edlib / vsearch / blastn / mmseqs**（複数指定可。欠けたバイナリは既定でスキップ）                              |
| 比較出力     | `summary_by_sample.tsv`（長形式）+ `summary_compare.tsv`（横断）                                   |
| 結果ツリー    | `amplicon_reassess/{backend}/{sample}/...`                                                |
| 実行タイミング  | Step 2b 後・任意                                                                              |




### 16.3 検証したい問い


| #   | 問い                                                                                    | 出力イメージ                                                   |
| --- | ------------------------------------------------------------------------------------- | -------------------------------------------------------- |
| Q1  | 既存コンセンサス同士はどの程度近いか？ユーザーが **採用／棄却を判断**できる類似度情報があるか？                                    | コンセンサス×コンセンサスの距離／identity、`pass_threshold`、近縁グループ（レビュー用） |
| Q2  | 最終クラスタに寄与しなかったリードは、既存コンセンサスのいずれかに **再割当可能**か？                                         | 各未寄与リードの最近傍コンセンサス ID、距離／identity、閾値内ヒット有無、ヒストグラム         |
| Q3  | 未寄与リード同士だけで **新規クラスタ候補**（既存に遠いが高支持の塊）があるか？                                            | 未寄与サブセットの簡易クラスタ要約（任意・計算コスト注意）                            |
| Q4  | Step 2b の `max_cluster_edit` / vsearch `id` / `min_cluster_reads` を変えた場合、再割当率がどう変わるか？ | 閾値スイープ表（感度解析）                                            |




### 16.4 「未寄与リード」の定義


| 層   | 定義                                         | 備考                                                                                          |
| --- | ------------------------------------------ | ------------------------------------------------------------------------------------------- |
| U0  | demux 割当済みだが `assign_status != "assigned"` | MVP 対象外。プライマー／長さ QC が主戦場                                                                    |
| U1  | gene 割当済みだが最終採用クラスタに含まれない                  | **MVP 中心**。assemble が `amplicon/{sample}/unassigned_to_cluster.tsv` に `read_id, reason` を書く |
| U2  | （任意）採用後に品質除外                               | 未実装                                                                                         |


`reason` 例: `below_min_cluster_reads`, `max_clusters_overflow`, `low_cluster_purity`, `overlap_unmatched`, `low_gene_reads`, `trim_or_empty`, `insufficient_after_trim`。  
TSV が無い旧 `amplicon/` に対しては、リード→consensus 距離からの推定フォールバックあり。

### 16.5 類似度バックエンド候補

数値化の中核は「クエリ配列 × ターゲット配列のアライメント距離または identity」。バックエンドは抽象化し、同一レポート形式に正規化する。

#### 16.5.1 ユーザ指定の候補


| ツール                | 出力の取り方                                          | 長所                                              | 短所・注意                                                                         | 本プロジェクトとの関係                                            |
| ------------------ | ----------------------------------------------- | ----------------------------------------------- | ----------------------------------------------------------------------------- | ------------------------------------------------------ |
| **edlib**（NW / HW） | 編集距離、任意で CIGAR                                  | demux / resolve で既に使用。依存追加が少ない。閾値を bp 単位で議論しやすい | 全対全は O(nm) で規模に弱い。長配列×多数クエリは要チャンク／近似                                          | **第一候補（小〜中規模・診断）**。Step 2b の `max_cluster_edit` と単位が揃う |
| **BLASTN**         | bitscore / e-value / pident / length            | 広く通じる報告言語。局所類似・部分ヒットの可視化に強い                     | 外部依存。閾値（e-value）は長さ依存で解釈が難しい。レガシー撤去方針（Phase B）と衝突しうるが **診断用オプション**としては再導入余地あり | 共同研究者向けレポートや外部ツール比較用                                   |
| **MMseqs2**        | easy-search / easy-cluster の identity・e-value 等 | 大規模に強い。クラスタ＋検索が同一系で揃う。ONT／メタゲノム系統で実績多い          | 外部バイナリ。パラメータ（感度 `-s`、被覆）の文書化が必要                                               | **大規模フルラン向け第一候補**                                      |




#### 16.5.2 その他妥当なツール（根拠つき）


| ツール                                                  | 根拠・適性                                                   | 備考                                                         |
| ---------------------------------------------------- | ------------------------------------------------------- | ---------------------------------------------------------- |
| **VSEARCH** (`--usearch_global` / `--cluster_fast`)  | Step 2b のクラスタリング本体と同一系統。再検証閾値を **現行** `id` **と直接比較**できる | 既に PATH 前提の推奨バックエンド。**再割当監査との一貫性で優先度高**                    |
| **Minimap2** (`-x map-ont` 等）                        | ONT リード長・エラープロファイル向けの実用標準マッパ。コンセンサスへリードを乗せる「再割当」に直観的    | identity は CIGAR / `--cs` から派生。クラスタ統合より **リード→コンセンサス割当**向き |
| **Parasail / Biostrings::pairwiseAlignment**         | edlib と同系の pairwise。R 内完結の代替                            | 速度は edlib C++ に通常劣る                                        |
| **USEARCH**                                          | VSEARCH の元系統。文献・旧パイプライン比較で用語を揃えやすい                      | ライセンス制約あり。新規依存としては VSEARCH 優先                              |
| **CD-HIT / CD-HIT-EST**                              | 同一性閾値クラスタの古典的標準。コンセンサス少数への「落とす」操作に近い                    | 短い差のハプロタイプ分離は苦手になりがち                                       |
| **Dedupe / BBTools** `clumpify`                      | 高類似リードの近傍探索・圧縮に強い実務ツール群                                 | 遺伝的解釈より前処理・負荷低減向き                                          |
| **Kmer 系（sourmash / Mash / MMseqs k-mer prefilter）** | 粗い候補絞りに有効。全対全の一次フィルタ                                    | 最終判定はアラインメント距離にフォールバック推奨                                   |
| **HMMER / Infernal** `nhmmer`                        | プロファイル参照がある場合の遠隔類似                                      | 本 amplicon パネルでは必須ではない                                     |
| **WFA2 / BiWFA**                                     | 長い ONT 配列の pairwise を高速化する編集団ライブラリ                      | edlib 代替として将来検討可                                           |


**推奨の使い分け（初期実装案）:**

1. **小規模診断・閾値スイープ:** edlib（既存）＋ 必要なら VSEARCH `--usearch_global`（Step 2b と同一 identity 定義）
2. **フルラン・サンプル多数:** MMseqs2 `easy-search`（未寄与 → コンセンサス）および `easy-cluster`（コンセンサス同士・未寄与サブセット）
3. **共同研究者・外部比較レポート:** BLASTN の pident / length を併記（解釈用。本線依存にはしない）
4. **リードをコンセンサスに「乗せる」操作の検証:** Minimap2 map-ont → 被覆率・不一致塩基数を再割当スコアに変換



### 16.6 出力契約

```
amplicon_reassess/                    # out_dir（元 amplicon/ とは別）
  run_stats.txt
  summary_by_sample.tsv               # n_consensus, n_unassigned_U1, reassignable_frac, n_near_duplicate_groups, ...
  summary_compare.tsv
  {backend}/{sample_id}/
    consensus_pairwise.tsv            # gene_id, cluster_i, cluster_j, tool, metric, value, pident, edit_distance, pass_threshold
    consensus_near_duplicate_groups.tsv  # レビュー用近縁グループ（自動マージしない）
    unassigned_to_consensus.tsv       # sample_id, gene_id, read_id, reason, layer, best_cluster_id, tool, metric, value, pident, edit_distance, pass_threshold
    metric_histograms.tsv             # 任意
    leftover_novel_clusters.tsv       # 任意（Q3）
```

- **元の** `amplicon/{sample}/consensus.fasta` **は変更しない。**
- **マージ配列は出力しない。** `consensus_pairwise.tsv` と近縁グループを見て、ユーザーが採用コンセンサスを判断する。

主要列の意味（例）:


| 列                   | 意味                                        |
| ------------------- | ----------------------------------------- |
| `metric`            | `edit_distance` / `pident` / `bitscore` 等 |
| `pass_threshold`    | Q1/Q2 のハイライト閾値を満たすか（自動統合フラグではない）          |
| `reassignable_frac` | U1 のうちいずれかのコンセンサスに `pass_threshold` した割合  |




### 16.7 API（実装済み）

```r
doReassessAssemblies(
  amplicon_out,                 # amplicon/ dir または assemble 戻り値（読み取り専用）
  gene_assign = ga,             # U1 配列復元に推奨
  out_dir = file.path(run_dir, "amplicon_reassess"),
  backend = c("edlib", "vsearch", "blastn", "mmseqs"),
  consensus_merge_max_edit = 12L,  # Q1 ハイライト閾値のみ（自動マージしない）
  read_assign_max_edit = 12L,
  min_identity = NULL,          # NULL なら edit から identity を推定
  skip_missing_backend = TRUE,  # PATH に無いツールは警告してスキップ
  primer_list = primer_list,
  overwrite = FALSE
)
```

出力の見方:


| ファイル                                                     | 内容                                                                     |
| -------------------------------------------------------- | ---------------------------------------------------------------------- |
| `summary_by_sample.tsv`                                  | sample × **backend** の `reassignable_frac` / `n_near_duplicate_groups` |
| `summary_compare.tsv`                                    | 同じ指標をバックエンド列で横並び（比較用）                                                  |
| `{backend}/{sample}/consensus_pairwise.tsv`              | 対の距離/`pident` + `pass_threshold`（主結果）                                  |
| `{backend}/{sample}/consensus_near_duplicate_groups.tsv` | 閾値内近縁グループ（レビュー用。配列は書かない）                                               |
| `{backend}/{sample}/unassigned_to_consensus.tsv`         | U1 の最近傍コンセンサス                                                          |


本線 `amplicon/` は変更しない。マージ fasta は出さない。

### 16.8 成功基準

1. 代表ラン（例: `tani_amplicon_20260709`）で、サンプルごとに `reassignable_frac`（U1）と近縁グループ数を算出できる
2. Step 2b の `max_cluster_edit` / vsearch `id` について、**レポート上は両尺度を併記**できること（「同一スケールに換算して比較する」ことは §17 #3 / §18-J11 で要検討。強制しない）
3. edlib と VSEARCH（または MMseqs2）で、同一サブセットの順位相関がおおむね一致すること（バックエンド健全性）
4. 実行後も元の `amplicon/{sample}/consensus.fasta` が変更されないこと。`merged_*.fasta` が出力されないこと
5. ドキュメント上、BLASTN は診断オプションであり本線必須依存ではないことが明記される



### 16.9 §15 との関係


|       | §15（editcall 偽陽性）       | §16（本節）                              |
| ----- | ----------------------- | ------------------------------------ |
| 失敗モード | 似すぎて **誤合流**            | 厳しすぎて **過剰分裂・過少召回**                  |
| 主目的   | 偽 allele の抑制            | 未寄与の内訳診断・閾値根拠・再割当候補・コンセンサス類似度の可視化    |
| 既定操作  | per-read `doEditcall` 等 | 診断レポートのみ（本線 `amplicon/` 非破壊・自動マージなし） |


両者は対立せず、**精度（§15）と召回／閾値説明（§16）を両側から支える**。

### 16.10 実装タスク分解


| 優先  | タスク                                                          | 状態    |
| --- | ------------------------------------------------------------ | ----- |
| 1   | U1 定義と assemble 時 `unassigned_to_cluster.tsv`                | **済** |
| 2   | edlib Q1/Q2 + `amplicon_reassess/`（本線非破壊）                    | **済** |
| 3   | **vsearch / blastn / mmseqs** バックエンド + `summary_compare.tsv` | **済** |
| 4   | 閾値スイープ（Q4）                                                   | 後続    |
| 5   | Q3 未寄与サブクラスタ                                                 | 後続    |
| 6   | 代表実データランでの回帰                                                 | 後続    |


---



## 17. 監査「中程度」所見への対応方針（2026-07-14）

コード監査で挙げた中程度所見について、議論のうえ次のとおり対応する。致命・高重大度の修正（spoa 固定、組合せ単端停止、indel bp 等）とは別レイヤの決定事項。

### 17.1 所見と決定一覧


| #   | 所見                                            | 決定                                                                                                              | 実装状態                                     |
| --- | --------------------------------------------- | --------------------------------------------------------------------------------------------------------------- | ---------------------------------------- |
| 1   | `require_unique_pair` は引数のみで未使用               | **削除**（API / roxygen / man / 計画書）                                                                               | **完了**（Phase J）                          |
| 2   | Dual は suffix が緩く「合法だが別人」を拾いやすい               | **バーコード／プレートデザインの問題として利用者責任**。1:1 または高 Hamming を勧める。組合せでは誤 dual → 合法ウェルがありうる旨を README に明記。ツール側で topology 救済はしない | README 反映済み。demux コード変更なし                |
| 3   | vsearch identity と内部絶対 edit は同値でない            | **要検討**として残す。Phase J では README 注意一文のみ（J11b 案2）。揃える・換算式改定は別 Phase                                                | 文書注意 **済**。技術揃えは **todo**                |
| 4   | `fraction_sample` は sample 全体ではなく gene バケツ内比率 | **実装は据え置き**。列名は変更しない。README / 計画書で定義を明記                                                                         | 計画書 §6.3・README・ampliconresolve_plan 更新済 |
| 5a  | PAM 染色体の `chr%02d` 強制                         | **廃止**。`pam_list` 列2は genome seqname と一致した文字列をユーザーが書く。不一致は **エラー停止**                                            | **完了**（Phase J）                          |
| 5b  | 複数ガイド時の先頭 PAM 固定                              | **案 A**（ガイド別窓）+ **案 A′**（ジョイント両 cut／介在欠失）。ケース1は gene 分割運用                                                       | **完了**（Phase J）                          |




### 17.2 問題2の補足（なぜツール非責とするか）

組合せプレートでは、正しい F に誤 R が付いても **別ウェルとして合法**になりうる（配列 Hamming だけでは防げない位相の問題）。miaoseq は:

- Hamming / `design_check.tsv` で **配列衝突**を可視化する
- 組合せ + `allow_single_end` は拒否する

一方、dual 経路での「合法クロスウェル」までは保証しない。精度が必要なら 1:1 レイアウトまたは十分大きい barcode 距離を選ぶ。

### 17.3 問題3の検討メモ（未決）

現状の換算: `id ≈ 1 - max_cluster_edit / median_len`（下限 0.70）、vsearch `--iddef 2`。長さヘテロや短い chimera で、絶対 edit と fractional identity の合格域がずれる。reassess の `approx_edit = (1-pident)*alnlen` も global NW ではない。

検討オプション（将来）:

1. 文書のみ（バックエンド間比較を非推奨）
2. vsearch も edit 尺度に寄せる（またはその逆）
3. 長さビン別閾値

いずれも本節の決定対象外（**要検討**）。

### 17.4 問題5・案 A / 案 A′ の実装チェックリスト

**案 A（ガイド別窓）**

- [x] `.build_editcall_metadata`: seqname をそのまま使用；genome 不一致で `stop`
- [x] 同 gene 複数 PAM 行を保持；`target_gene = gene_guide`
- [x] primers / amplicon は gene で共有；欠損は `stop`
- [x] `.reads_to_edit_records`: `pos_idx[1]` 廃止 → PAM 行ループ；insert align **1回/リード**
- [x] guide 欠落で複数行 / guide 重複 → `stop`
- [x] `inst/extdata/agr8_pam_list.csv` を `chr01` 等に更新
- [x] README / `doEditcall` man / `editcall_adapter.md` を §7.3–7.4 に合わせて更新
- [x] 単体テスト: 同 gene 2 PAM → 2× `target_gene`；seqname 不一致 → 停止（`tests/test-editcall-phase-j.R`）

**案 A′（ジョイント両 cut・§7.4）**

- [x] 隣接ガイド対の insert 座標 `c_i < c_j` と `expected_span`
- [x] 外側錨 + `del_span ≈ expected` → `both_cut_excision`（`min_span_bp=30` / `excision_tol_bp=20`）
- [x] `event_class`: `wt` / `g_i_only` / `g_j_only` / `both_local` / `both_cut_excision`
- [x] 出力: `editcall_joint.csv`, `editcall_joint_summary.csv`
- [x] 小窓片側 discard でも excision が取れる回帰テスト（合成大欠失リード）
- [x] 近接 cut（`< min_span_bp`）では excision を出さないテスト



### 17.5 関連ドキュメント


| ファイル                                               | 役割                                                     |
| -------------------------------------------------- | ------------------------------------------------------ |
| [README.md](../README.md)                          | 組合せ demux・`fraction_sample`・PAM・vsearch 注意（Phase J 同期） |
| [demux_revise.md](demux_revise.md)                 | §0 現行方針で prefix / `require_unique_pair` **削除済** を同期    |
| [editcall_adapter.md](editcall_adapter.md)         | pam・案 A / A′・joint 出力（Phase J 実装済）                     |
| [ampliconresolve_plan.md](ampliconresolve_plan.md) | `fraction_sample` 定義訂正・spoa 固定注記                       |
| **§18（本ファイル）**                                     | Phase J **実装前**の未決記録・回答欄（歴史；現状は実装済み）                   |


---



## 18. Phase J 実装前の未決事項（説明と回答欄）

作成: 2026-07-14  
目的: 方針（§17）は決まっているが、**実装仕様として曖昧な点**を一つずつ明確化し、回答をここに記入できるようにする。  
使い方: 各項目の「回答欄」に決定内容を書く。回答が揃った項目から実装に着手する。

> **ステータス（2026-07-14）:** 本節の回答に基づき **Phase J は実装済み**。以下の質問・回答欄は決定記録として残す。実装の正は §0 / §7.3–7.4 / `R/editcall.R`・`R/demultiplex.R`。

---



### J1. PAM 座標が「何を指すか」と、insert 上の点への写像

**何の話か（わかりやすく）**

`pam_list` の第3列には数字（ゲノム座標）が入る。Editcall ではその位置の周辺を切り出して「編集されたか」を見る。ところが実装には **二段の座標変換**が要る。

1. **ゲノム座標**（染色体上の位置）
2. → **amplicon / insert 上の位置**（プライマーで挟んだ PCR 断片の中での塩基番号）
3. → さらに §15.7 の **編集窓** `[start, end]` や、案 A′ の cut 点 `c_i`

いまのコードは「PAM 周辺の短い intact 配列を amplicon 全体にアラインし、ギャップの切れ目から窓を推定する」ヒューリスティックで、**1 gene・1 PAM** 向け。同 gene に PAM が2つあると名前衝突や誤った窓になりやすい。

加えて第3列が **PAM モチーフの先頭**なのか **Cas9 切断点（通常 PAM の約 3 bp 上流）**なのかが文書に無い。1 bp ずれると `expected`（cut 間距離）や小窓の中心がずれる。

**計画に既にあること**

- 列2は genome seqname と一致（`chr%02d` しない）  
- 案 A′ は `c1 < c2` と `expected ≈ c_j - c_i` を使う

**決めてほしいこと**

1. 第3列の生物学的意味（PAM 開始 / cut / その他）
2. ゲノム → insert 座標の正規手順（推奨例: genome 座標を amplicon ゲノム区間に相対化し、プライマー長を引いて insert 0-based/1-based を明示）
3. 旧「ギャップ島ヒューリスティック」を捨てるか、補助にするか

回答欄 J1（記入用）

```
決定:
  第3列の意味:PAM 開始
  写像手順:genomeにamplicon配列をマップして、ゲノム座標をアンプリコン上の座標に写像する。insert配列をアンプリコンにアラインして、アンプリコン上の座標を起点に窓を指定する。
  旧ヒューリスティック:捨てる
理由・補足:
日付:
記入者:
```

---



### J2. メタデータ構築：同 gene・複数 PAM を名前衝突させない手順

**何の話か**

案 A では gene=`SD1` にガイド `g1`,`g2` の2行を許す。一方 `amplicon.fa` と primers は `SD1` **が1本**。現行 `.build_editcall_metadata` は intact 配列の名前を gene 名にしたあと amplicon 名と `intersect` するため、**同名が2つあると潰れる／後勝ち**になる。

必要なデータの種類は実は2つ:


| キー                     | 例        | 使う場所                                      |
| ---------------------- | -------- | ----------------------------------------- |
| **gene（共有）**           | `SD1`    | primer trim、amplicon 参照、gene 割当 `gene_id` |
| **target_gene（ガイド単位）** | `SD1_g1` | 編集窓・集計・`intact_seq.fa` 名                  |


「いつ一意名にするか」を手順書に落とさないと、実装者が別々のやり方で壊れうる。

**決めてほしいこと**

1. 内部での一意化タイミング（例: PAM 行を読んだ直後に `target_gene` を確定し、amplicon は常に `target`=gene で紐づけ）
2. guide 列が空で gene が1行だけ → `target_gene = gene` でよいか
3. pam の gene が primer/amplicon に無いとき **即 stop** で確定か（§7.4 に書いてあるが確認）

回答欄 J2（記入用）

```
決定:
  一意化タイミング:PAM 行を読んだ直後に target_gene を確定し、amplicon は常に target=gene で紐づけ
  guide 欠落・単一行:target_gene = gene でよい
  primer/amplicon 不一致時: 即 stop で確定
理由・補足:
日付:
記入者:
```

---



### J3. 両 cut 欠失（excision）のとき、ガイド別出力（案 A）に何を書くか

**何の話か**

案 A′ が「このリードは両カットで真ん中が飛んだ（`both_cut_excision`）」と判定しても、案 A の小窓は左右の錨が取れず **片側または両側 discard** しがち。そのとき `editcall_all.csv`（ガイド効率の主表）に:

- 何も書かない（効率分母から消える）  
- discard として数える  
- 「excision があったのでこのガイドは編集扱い」と残す

のどれにするかで、**ガイド別編集率が大きく変わる**。

例: 100 リード中 40 が excision。小窓が全部失敗すると、案 A だけ見ると「編集なし／データ不足」、案 A′ だけ見ると「40% excision」になる。両方を見て初めて全体像が分かるか、案 A 側にもフォールバックで載せるかを決める必要がある。

**決めてほしいこと**

1. excision リードを案 A の各ガイド行にどう反映するか（載せる／載せない／フラグ列）
2. `event_class` の優先順位（例: 外側錨 OK かつ span 一致なら常に `both_cut_excision` を最優先）
3. 片方の小窓だけ成功＋介在欠失、のような中間ケースの扱い

回答欄 J3（記入用）

```
決定:
  案 A（editcall_all）への反映:載せる。excisionの場合は、各ガイドで編集が起こったものとしてカウント
  event_class 優先順位:外側錨 OK かつ span 一致なら常に both_cut_excision を最優先
  中間ケース:both_cut_excisionとして扱う。
理由・補足:
日付:
記入者:
```

---



### J4. 案 A′ のパラメータ・出力粒度・excision_rate の分母

**何の話か**

§7.4 には「±20 bp または 20%」「`min_span_bp` 例: 50」と例示があるが、`doEditcall` の引数になっていない。また:

- `editcall_joint.csv` が **リード1行ずつ**なのか、**集計済み**なのか曖昧  
- `excision_rate = both_cut_excision / ???` の `???` が未定義

分母の候補:


| 分母                      | 意味                            |
| ----------------------- | ----------------------------- |
| gene 割当できた全リード          | 失敗・discard も分母に入る → rate は控えめ |
| ジョイント判定に成功したリードだけ       | 技術失敗を除外 → rate は高め            |
| 案 A の少なくとも1ガイドで窓が取れたリード | 中間                            |


**決めてほしいこと**

1. 公開パラメータにするか、固定定数でよいか（値の推奨）
2. `editcall_joint.csv` = リード単位、`editcall_joint_summary.csv` = 集計、でよいか
3. `excision_rate`（および他 rate）の分母定義

回答欄 J4（記入用）

```
決定:
  パラメータ（min_span_bp / 許容幅 / 引数化の有無）:公開パラメータ
  joint CSV の粒度:editcall_joint.csv = リード単位、editcall_joint_summary.csv = 集計、でよい
  excision_rate 分母:案 A の少なくとも1ガイドで窓が取れたリード
理由・補足:
日付:
記入者:
```

---



### J5. ガイドの「隣接」の定義（ソートキー）

**何の話か**

案 A′ は当面「隣接ペアのみ」を見る。3ガイド `g1,g2,g3` があるとき、隣接とは通常 **insert 座標順**の隣どうし（`c1–c2`, `c2–c3`）であり、**CSV の行の順番**ではない。行順だとゲノム上の並びと食い違う。

**決めてほしいこと**

1. ソートキー = insert（または amplicon）座標の昇順で確定か
2. 座標が同値・極近接のときの tie-break（guide_id 文字列など）
3. 初期実装で非隣接ペア（`g1–g3`）を出さないことを確認

回答欄 J5（記入用）

```
決定:
  ソートキー:insert（または amplicon）座標の昇順で確定
  tie-break:guide_id 文字列
  非隣接ペア:非隣接ペア（g1–g3）を出さない
理由・補足:
日付:
記入者:
```

---



### J6. evalMiao / editViewer の Phase J スコープ

**何の話か**

案 A′ の結果を plate PDF や `evalMiao` の編集率にすぐ乗せるか、**まずは CSV だけ**出して viewer は後回しか。後者なら実装が軽いが、「主成果が viewer に出ない」期間が出る。

**決めてほしいこと**

1. Phase J MVP: joint CSV のみ / evalMiao 集計まで / editViewer まで
2. ガイド別 summary（案 A）は既存 viewer のまま動かすか

回答欄 J6（記入用）

```
決定:
  Phase J で必須の下流:joint CSV のみ
  viewer / evalMiao:既存 viewer のまま
理由・補足:
日付:
記入者:
```

---



### J7. 実装上のキャッシュキー（マルチガイド）

**何の話か**

現行 `.edit_seq_from_insert` はキャッシュキーに `gene_id` を使う。同一 insert 配列でもガイドが違えば窓が違うため、**キーに** `target_gene`**（または PAM 座標）を入れないと**、先に計算した別ガイドの窓を誤って再利用する。

これはほぼ実装バグ予防の固定仕様で、方針選択の幅は狭い。

**決めてほしいこと**

キーに `target_gene` を含める方針でよいか（推奨: はい）。

回答欄 J7（記入用）

```
決定:はい
理由・補足:
日付:
記入者:
```

---



### J8. 対応ヌクレアーゼ・対象外の明示

**何の話か**

座標オフセットや「両 cut → 介在欠失」モデルは **SpCas9 的な二重 DSB** を想定している。塩基エディタやプライムエディタ、片面だけ nick する系では「excision」意味が違う。非サポートを書いておかないと誤用される。

**決めてほしいこと**

1. サポート対象（例: Cas9 DSB 想定の二重ガイド）
2. 明示的非サポート

回答欄 J8（記入用）

```
決定:
  サポート:Cas9 DSB 想定の二重ガイド
  非サポート:それ以外は措定していないことを明示
理由・補足:
日付:
記入者:
```

---



### J9. 破壊的変更の移行（extdata・ユーザー pam_list）

**何の話か**

`chr%02d` 廃止後、同梱 `agr8_pam_list.csv` の `1,4,5…` は **そのままでは動かない**（genome は `chr01` 等）。ユーザー既存 CSV も同様。エラーメッセージだけで足りるか、マイグレーションノート／変換例が要るか。

**決めてほしいこと**

1. extdata を `chr01` 形式に更新する確認
2. エラー文面に「genome の header と一致させよ」例を出すか
3. 数値だけ書いた旧形式の自動変換は **しない** でよいか（現行方針どおり）

回答欄 J9（記入用）

```
決定:
  extdata 更新:更新する
  エラーメッセージ:genome の header と一致させよと出す。
  旧形式の自動変換:しない
理由・補足:
日付:
記入者:
```

---



### J10. Phase J の作業順序（パッケージ分割）

**何の話か**

依存関係を間違えると二重作業になる。推奨順の確認用。


| 順   | 作業                       | 依存             |
| --- | ------------------------ | -------------- |
| 1   | §18 回答確定                 | —              |
| 2   | `require_unique_pair` 削除 | 独立             |
| 3   | PAM seqname + extdata    | J1,J9          |
| 4   | 案 A（メタ＋ガイド別窓）            | J1,J2,J7       |
| 5   | 案 A′（joint）              | 案 A + J3,J4,J5 |
| 6   | 下流レポート（任意）               | J6             |
| 7   | 関連ドキュメント最終同期             | 全体             |


**決めてほしいこと**

この順序でよいか。2 と 3 の前後入れ替え可否。

回答欄 J10（記入用）

```
決定:
  採用する順序:この順序でOK
  変更点:無し
理由・補足:
日付:
記入者:
```

---



### J11. vsearch identity と絶対 edit（§17 #3）の「完了」条件

**何の話か**

「要検討」のままだと、いつまでもタスクが残る。いま取りうる閉じ方:

1. **文書注意のみで完了**（バックエンド間の数値比較を非推奨と README / man に書く）
2. 将来、換算式や長さビンを検討する（別 Phase）
3. 今すぐ揃える実装に入る

§16 の成功基準は「両尺度を併記」に緩め済み（本改訂）。

**決めてほしいこと**

Phase J の範囲では 1 で閉じるか。

回答欄 J11（記入用）

```
決定:別 Phase
理由・補足:
日付:
記入者:
```

---



### J12. 既に直した文書同期（参考・回答不要でも確認可）

以下は 2026-07-14 時点で同期済み。追加要望があれば回答欄へ。


| ファイル                      | 内容                                                        |
| ------------------------- | --------------------------------------------------------- |
| `demux_revise.md`         | §0 現行方針、prefix 廃止、`require_unique_pair` **削除済**、組合せリスク    |
| `editcall_adapter.md`     | pam / 案 A・A′・joint（Phase J 実装済）                           |
| `ampliconresolve_plan.md` | `fraction_sample` 訂正、spoa 固定注記                            |
| `README.md`               | 組合せ demux、`fraction_sample`、PAM seqname/strand、vsearch 注意 |
| 本ファイル §4                  | `editcall_joint*.csv` をツリーに追加                             |
| 本ファイル §16.8               | edit/vsearch「同一尺度」強制を緩和                                   |


回答欄 J12（追記要望があれば）

```
追加で直したい文書・箇所:特になし
日付:
記入者:
```

---



### 18.1 記入レビュー結果（2026-07-14）

ユーザー記入を検証した。**大方針は一貫しており実装可能なレベルに近い。** ただし次は実装前に数値・文字列まで落とす必要がある（下記 J1b / J3b / J4b / J5b / J11b）。


| ID  | 判定    | メモ                                                                               |
| --- | ----- | -------------------------------------------------------------------------------- |
| J1  | △ 要追記 | 「PAM 開始」「写像の骨格」「旧ヒューリスティック捨てる」はよい。写像の具体手順・strand・窓中心が未定 → **J1b**                |
| J2  | ○     | 一意化タイミング・単一行・即 stop で十分。複数行かつ guide 欠落は既存 §7.4（停止）で読める                           |
| J3  | △ 要追記 | excision を案 A に「載せる／編集扱い」はよい。`read_seq` に何を書くか未定 → **J3b**。J4 分母と衝突しうる → **J4b** |
| J4  | △ 要追記 | 公開パラメータ・CSV 粒度はよい。**デフォルト値と許容幅の単位が無い** → **J4b**                                 |
| J5  | △ 軽微  | ソート・tie-break・非隣接オフはよい。「insert **または** amplicon」が二択のまま → **J5b**                 |
| J6  | ○     | joint CSV のみ／viewer は既存（案 A summary）で明確                                          |
| J7  | ○     | `target_gene` をキーに含める                                                            |
| J8  | ○     | Cas9 DSB 二重ガイド；それ以外は想定外と明示（文面「措定」→実装時は「想定」）                                      |
| J9  | ○     | extdata 更新・一致を促すエラー・自動変換なし                                                       |
| J10 | ○     | 推奨順のまま                                                                           |
| J11 | △ 要確認 | 「別 Phase」＝Phase J では触れない、でよいか。文書一行も不要か → **J11b**                                |
| J12 | ○     | 追加なし                                                                             |


**J3 ↔ J4 の潜在矛盾（重要）**  
J3: 小窓失敗でも excision なら案 A に「編集」として載せる。  
J4: `excision_rate` 分母 =「案 A の少なくとも1ガイドで**窓が取れた**リード」。  
→ 両窓失敗＋excision のリードは案 A に載るのに「窓が取れた」には入らず、分母から落ちうる。J4b で定義を揃えること。

---



### J1b. 写像手順の具体化（J1 追記）

**なぜまだ曖昧か**

J1 の「genome に amplicon をマップ → ゲノム座標を amplicon 座標へ → insert を amplicon にアラインして窓」は正しい骨格だが、実装が分岐しうる点が残る。

1. **amplicon ↔ genome の載せ方**
  `prepAmpliconDB` と同様の完全一致（`matchPattern`）か、既存の pairwise か。多重ヒット時は既に prep で止める前提か。
2. **鎖**
  amplicon が genome と逆向きのとき、PAM（genome +）を amplicon 座標へどう折り返すか。
3. **窓の中心**
  第3列 = PAM 開始で確定。編集窓は `PAM_start ± check_window` を **amplicon 座標**で切るのか。Cas9 cut（PAM の約 3 bp 上流）へのオフセットは **入れない**のか（J8 は Cas9 想定だが J1 は PAM 開始）。
4. **案 A′ の** `c_i`
  PAM 開始の amplicon/insert 座標そのものでよいか（cut にずらさないか）。
5. **「insert を amplicon にアライン」**
  ここでの insert は（a）プライマー付き amplicon 参照からプライマーを除いた ref_insert か、（b）リードの primer-trim 後配列か。通常は参照 ref_insert 上に窓を定義し、リード側は同じ gapped 座標で切る。

**決めてほしいこと**

上記 1–5 を一文ずつ指定。

回答欄 J1b（記入用）

```
決定:
  amplicon↔genome の載せ方:prepAmpliconDB と同様の完全一致。多重ヒット時は既に prep で止める。
  鎖（RC）の扱い:amplicon が genome と逆向きのとき、ampliconをゲノムと同じ向きにする。
  編集窓の中心・幅（check_window との関係）:編集窓は PAM_start ± check_window を amplicon 座標で切る
  Cas9 cut オフセット（0 bp でよいか）:pam_listのpam開始は＋ー鎖の区別が無い。列5をchain列として+なら-3 bp、-なら+ 3bpのオフセットを追加する。
  案 A′ の c_i の定義:cut にずらす。
  「insert」の指す配列:ref_insert 
理由・補足:
日付:
記入者:
```

---



### J3b. excision を案 A に載せるときの `read_seq` / `intact`（J3 追記）

**なぜまだ曖昧か**

「各ガイドで編集が起こったものとしてカウント」は集計方針としてよいが、現行パイプラインは `read_seq` 文字列の **完全一致**でアレルを束ねる。載せ方によって summary の見え方が変わる。

候補例:


| 案   | `read_seq` に書くもの                              | 利点            | 欠点                  |
| --- | --------------------------------------------- | ------------- | ------------------- |
| A   | 固定トークン（例: `EXCISION` または `both_cut_excision`） | 集計が単純、ガイド間で揃う | 実際の接合配列は残らない        |
| B   | ジョイント外側錨で切った接合アレル（gapped 可）                   | 配列が見える        | ガイド行に同じ長アレルが二重計上される |
| C   | 各ガイド小窓が取れたときだけ実配列、失敗時のみトークン                   | 混在を最小化        | 実装分岐が増える            |


また `intact` は常に `FALSE` か。`count` はガイドごとに +1 か（同一リードが g1 と g2 の両方に乗る → **リードは二重計上**。J3 の意図どおりか）。

**決めてほしいこと**

1. `read_seq` / `ref_seq` の中身（A/B/C または別案）
2. `intact`
3. 同一リードが複数 `target_gene` 行に乗る二重計上を許容するか

回答欄 J3b（記入用）

```
決定:
  read_seq / ref_seq:B
  intact:FALSE
  ガイド間の二重計上:許容する
理由・補足:
日付:
記入者:
```

---



### J4b. パラメータ既定値と excision_rate 分母の再定義（J4 追記）

**なぜまだ曖昧か**

1. 「公開パラメータ」とだけあり、**名前・デフォルト・単位**が無いと実装できない。
  例が必要: `min_span_bp = 50`、`excision_tol_bp = 20` または `excision_tol_frac = 0.2`（どちらを使うか／両方か）。
2. J3 との整合: 小窓失敗でも excision を案 A に載せるなら、分母「窓が取れたリード」は不適切。
  推奨の言い直し例:  
  - **分母** = 当該 gene でジョイント判定を試行し、`event_class` が付いたリード（discard 除く）  
  - または = 案 A に **いずれかの target_gene 行が1行でも出た**リード（excision トークン含む）

**決めてほしいこと**

1. 引数名とデフォルト（数値）
2. 許容幅は bp 固定 / 割合 / `max(bp, frac×expected)` のどれか
3. `excision_rate`（および `both_local` 等の rate があれば）の分母を J3 と矛盾しない文言で

回答欄 J4b（記入用）

```
決定:
  引数名とデフォルト:min_span_bp = 30、excision_tol_bp = 20
  許容幅の定義:bp 固定 
  excision_rate 分母（J3 整合版）: 案 A に いずれかの target_gene 行が1行でも出たリード
理由・補足:
日付:
記入者:
```

---



### J5b. ソートに使う座標系（J5 追記）

**なぜまだ曖昧か**

「insert（または amplicon）座標」だと実装者がどちらか好きに選べる。J1 の写像で PAM を一度 **amplicon 座標**に載せ、primer 長を引いて **insert 座標**にするなら、隣接ソートはその **どちらに統一した座標**で行うかを決める。

推奨: **insert 座標（primer trim 後の参照）で昇順**（窓抽出と同じ系）。

回答欄 J5b（記入用）

```
決定:
  隣接ソートに使う座標:insert 座標（primer trim 後の参照）で昇順
理由・補足:
日付:
記入者:
```

---



### J11b. 「別 Phase」の意味（J11 追記）

**なぜまだ曖昧か**

「別 Phase」は次のどれにも読める。

1. Phase J では **何もしない**（README 注意も書かない）。vsearch↔edit は将来チケット
2. Phase J では **文書一行だけ**書いて技術検討は別 Phase
3. 別 Phase で **実装までやる**予定（換算式など）

Phase J の完了条件が変わる。

回答欄 J11b（記入用）

```
決定:
  Phase J での扱い（1 / 2 / 3 または別）:2
理由・補足:
日付:
記入者:
```

---



### 18.3 第2回レビュー（J1b–J11b 記入後・2026-07-14）

**総評:** 追記回答で実装ブロッカーは解消。J3↔J4 の分母矛盾も J4b で解消。**J1c 記入済み → Phase J 着手可。**


| ID   | 判定         | メモ                                                   |
| ---- | ---------- | ---------------------------------------------------- |
| J1b  | ○（J1c で補完） | matchPattern・RC・`ref_insert`・`c_i`=cut               |
| J1c  | ○          | 小窓中心も **cut**。列5欠落はオフセットなし。オフセットは **genome 上**       |
| J3b  | ○          | 案 B・`intact=FALSE`・二重計上許容                            |
| J4b  | ○          | `min_span_bp=30`・`excision_tol_bp=20`。分母＝案 A 行が出たリード |
| J5b  | ○          | insert 座標昇順                                          |
| J11b | ○          | Phase J は文書一行                                        |


~~残課題（J1c）~~ → **解消済み**（下記 J1c 回答・レビュー参照）。

---



### J1c. 編集窓の中心と pam_list 列5（J1b 追記）

**決めてほしいこと**

1. 案 A の §15.7 編集窓の中心は **PAM 開始のまま**か、**cut（列5オフセット後）**か
2. 列5（鎖 / `+`|`-`）の位置と必須性（例: `gene, seqname, pam_start, guide, strand`）。欠落時の挙動
3. オフセットの適用座標（genome 上でずらしてから amplicon へ写す／amplicon 上でずらす）— J1b の「先に amplicon を genome 向きに揃える」前提ならどちらでも等価になりやすいが、一文で固定したい

回答欄 J1c（記入用）

```
決定:
  案 A 編集窓の中心（PAM / cut）:cut
  列5の必須性・欠落時:欠落時はpam_startをそのまま使う
  オフセットを掛ける座標（genome / amplicon）:genome 上でずらす
理由・補足:
日付:
記入者:
```

**J1c レビュー:** 問題なし。案 A 小窓も A′ の `c_i` も **cut 中心で一致**。列5欠落時はオフセットなし（従来の PAM 座標中心）。オフセットは genome 上で適用してから amplicon へ写像。列5が不正値（`+`/`-` 以外）のときは実装時 `stop` でよい。

---



### 18.2 回答状況チェックリスト（記入後に更新）


| ID   | 題目                     | 状態                 |
| ---- | ---------------------- | ------------------ |
| J1   | PAM 意味と座標写像            | **記入済・問題なし**       |
| J1b  | 写像の具体化                 | **記入済・問題なし**       |
| J1c  | 編集窓中心と strand 列        | **記入済・問題なし**       |
| J2   | メタデータ一意化手順             | **記入済・問題なし**       |
| J3   | excision 時の案 A 連携      | **記入済・問題なし**       |
| J3b  | excision 載せるときの配列表現    | **記入済・問題なし**       |
| J4   | A′ パラメータ・粒度・分母         | **記入済・問題なし**       |
| J4b  | 既定値と分母の J3 整合          | **記入済・問題なし**       |
| J5   | 隣接のソート定義               | **記入済・問題なし**       |
| J5b  | insert vs amplicon     | **記入済・問題なし**       |
| J6   | evalMiao / viewer スコープ | **記入済・問題なし**       |
| J7   | キャッシュキー                | **記入済・問題なし**       |
| J8   | ヌクレアーゼ対象範囲             | **記入済・問題なし**       |
| J9   | extdata / 移行           | **記入済・問題なし**       |
| J10  | 作業順序                   | **記入済・問題なし**       |
| J11  | vsearch↔edit の閉じ方      | **記入済・問題なし**       |
| J11b | 「別 Phase」の解釈           | **記入済・問題なし**（文書1行） |
| J12  | 文書同期の追加要望              | **記入済（なし）**        |


**結論（2026-07-14）:** §18 の未決回答は揃い、**Phase J 実装まで完了**した。追加の回答欄は不要。

実装反映済み: §7.3 列5（strand）、`agr8_pam_list.csv`（`chr01` 等）、README / man / `editcall_adapter.md`、`tests/test-editcall-phase-j.R`。残る任意作業は実データ回帰と joint UI（J6）。

---



## 19. Phase K: クラスタリング/コンセンサス基盤の刷新



### 19.1 背景と動機

tani 16S データ（~1,500 bp × 57 サンプル）の解析で以下が判明した:

1. `cluster_backend = "internal"`**（edlib greedy）は長鎖アンプリコンに不適**: 絶対編集距離 12 bp は 1,500 bp 配列に対して 99.2% 同一性に相当し、ONT basecall 揺らぎにより大半のリードがクラスタ閾値を超えて細分化される。結果、57 サンプル中 7 サンプルでコンセンサスが得られず、クラスタ総数も vsearch の 289 に対し 74 に留まった。
2. **spoa コンセンサスの問題**: 旧 spoa 経路ではキメラ様の異常長コンセンサス（>1,600 bp、一部 >2,900 bp）が一定割合で生じた。abPOA への切替で長さ分布が正常化し、コンセンサス品質が改善した。
3. **代替バックエンドの評価が必要**: MMseqs2（connected component / linclust）、abPOA `-Q`（品質重み）、Racon 1 パスを比較し、正式採用を決定する必要があった。



### 19.2 正式採用（2026-07-16 確定）


| 項目      | 正式採用                                    | 廃止                                                         |
| ------- | --------------------------------------- | ---------------------------------------------------------- |
| クラスタリング | `vsearch`（`--cluster_fast`、identity 閾値） | `internal`（edlib greedy）、`mmseqs`（easy-cluster / linclust） |
| コンセンサス  | `abpoa -r 0 -m 0`（FASTA 入力、`-Q` なし）     | `spoa`、abPOA `-Q`（FASTQ 品質重み）                              |
| ポリッシュ   | **なし**                                  | Racon（minimap2 + racon 1 パス）                               |


**呼び出し:**

```bash
# クラスタリング（実装どおり）
vsearch --cluster_fast input.fa --id <min_cluster_identity> \
  --strand both --iddef 2 --uc clusters.uc

# コンセンサス
abpoa -r 0 -m 0 input.fa
# -r 0: consensus FASTA
# -m 0: global alignment
```



### 19.3 評価結果（不採用の根拠）



#### A. MMseqs2 connected component（`easy-cluster --cluster-mode 1`）

- フルラン比較では vsearch と生物学的に近いコンセンサスが得られたが、クラスタ数が増えプロセス起動オーバーヘッドで大幅に遅かった（vsearch ~44 min vs mmseqs ~145 min、tani フル）。
- 速度・運用コストの面で vsearch を上回る利点がなかった。



#### B. MMseqs2 linclust（`easy-linclust`）

- 32 サンプル比較: 割り当て率 **vsearch 13.9% vs linclust 3.6%**（linclust は約 1/5）。
- 速度は約 2.2× だが、k-mer 近似により感度が不足。アンプリコン（数千リード規模）には不適。



#### C. abPOA `-Q`（品質重み）

- 42 クラスタ比較: `-Q` あり/なしのコンセンサス差は平均 0.07 edit（ほぼ同一）。
- 十分なリード数では多数決 POA が品質重みを打ち消す。追加の FASTQ 入出力コストに見合う効果がなかった。



#### D. Racon 1 パス

- 42 クラスタ比較: メンバーへの平均 edit 距離が **abpoa 単独 12.0 → +racon 14.8** と悪化。
- minimap2 の末端 soft-clip によりコンセンサスが数 bp 短縮され、既に良好な POA コンセンサスを劣化させた。
- Racon は noisy de novo assembly 向けであり、abPOA 後の amplicon ポリッシュには不適。



### 19.4 API 変更

```r
doAssembleAmplicons(
  gene_assign,
  out_dir,
  primer_list = NULL,
  method = c("both", "cluster", "consensus"),
  min_reads = 5L,
  min_cluster_reads = 5L,
  max_clusters = Inf,
  min_cluster_identity = 0.95,
  max_primer_edit = 5L,
  cluster_backend = "vsearch",   # 固定（他バックエンド廃止）
  assembly_backend = c("cluster", "overlap_graph"),
  overlap_min_identity = 0.90,
  strict_end_trim = FALSE,
  min_cluster_purity = NULL,
  overwrite = FALSE,
  n_core = 1L
)
```



### 19.5 廃止する機能・依存


| 対象                                                     | 処置                                |
| ------------------------------------------------------ | --------------------------------- |
| `greedy_cluster_cpp` / `cluster_backend = "internal"`  | 削除済み                              |
| `cluster_backend = "mmseqs"` / `.cluster_ids_mmseqs()` | **廃止**（コード・API から除去）              |
| `.consensus_via_spoa()` / spoa バイナリ                    | 削除済み（abpoa 置換）                    |
| abPOA `-Q`（FASTQ 品質重み）                                 | **不採用**（FASTA モードのみ）              |
| `.polish_via_racon()` / minimap2 / racon               | **廃止**（コンセンサス経路から除去）              |
| mmseqs2 / racon / minimap2 コンテナ依存                      | Assemble 本線から外す（任意検証用に残す場合は文書で明示） |


**残す依存:** `vsearch`、`abpoa`（Apptainer ラッパー経由）。

### 19.6 実装手順（確定後の残作業）

1. `pull_sifs.sh` ~~に abpoa を追加~~（済）
2. `.consensus_via_abpoa()` ~~を実装（spoa 置換）~~（済）
3. `greedy_cluster_cpp` ~~/ spoa 関連コードを削除~~（済）
4. `cluster_backend` ~~を~~ `"vsearch"` ~~固定にし、~~`mmseqs` ~~分岐・~~`.cluster_ids_mmseqs()` ~~を削除~~（済）
5. `.consensus_for_members()` ~~から racon 呼び出しを削除（abpoa のみ返す）~~（済）
6. ~~DESCRIPTION / README /~~ `.tool_versions` ~~から mmseqs / racon / minimap2 の本線依存を除去~~（済）
7. ~~tani データで **vsearch + abpoa** の回帰確認~~（済: sample01–03, 28 clusters, method=abpoa）



### 19.7 期待される効果（確定構成）


| 指標       | 旧 (vsearch + spoa) | 正式採用 (vsearch + abpoa) |
| -------- | ------------------ | ---------------------- |
| コンセンサス精度 | キメラ様長配列あり          | 改善（長さ分布正常化）            |
| コンセンサス速度 | 基準                 | 2–15 倍高速 (abpoa)       |
| クラスタ感度   | 良好                 | 同左（vsearch 維持）         |
| 運用依存     | spoa               | **vsearch + abpoa のみ** |

### 19.8 精査 20260716 反映（API・I/O 改訂）

`codereview/20260716_1056_codereview.md` のユーザー方針に基づく追加確定:

| 項目 | 変更 |
|------|------|
| 閾値パラメータ | `max_cluster_edit` 廃止 → **`min_cluster_identity`（既定 0.95）** を vsearch `--id` に直接渡す。完了後に identity スイープ比較テストを行う |
| `max_clusters` | 既定 **`Inf`**（上位打ち切りなし）。除外は `min_cluster_reads`（既定 5）のみ |
| 列名 | `fraction_sample` → **`fraction_bucket`**（後方互換なし。README・計画・旧レポートを一括更新） |
| `clusters.fasta` | メンバー配列。ヘッダ `>{read_id} {cluster_id}`（スペース区切り） |
| quals | Assemble コンセンサス経路から受け渡し削除（`-Q` 不採用の徹底） |
| `run_stats.txt` | `r_version` / `miaoseq_version` + `min_cluster_identity` / `max_clusters` / `min_cluster_reads` 等 |
| テスト | `tests/test-assemble-phase-k.R`（バイナリ無ければ skip） |
| 文書 | `ampliconresolve_plan.md` を Phase K 同期；旧 spoa 節は「歴史・Phase C/D 当時」と注記；壊れた `test_abpoa_quality.R` は削除 |

**経路の役割（R4）:** Assemble = フル長代表（クラスタ内均しあり）。Editcall = 既知領域の局所パターン。16S 全体像に Editcall 局所分類を使うのはミスリード。ONT アンプリコンの代表復元の不確実性は分野共通で完全解はない。

---

