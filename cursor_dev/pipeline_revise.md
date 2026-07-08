# パイプライン再設計計画

作成日: 2026-07-08  
関連:
- [cursor_dev/demux_revise.md](demux_revise.md) — demux 詳細
- [cursor_dev/ontbarcoder_demux.md](ontbarcoder_demux.md) — ONTbarcoder 参照
- [cursor_dev/code_reorganize_plan.md](code_reorganize_plan.md) — 過去のファイル分割（参考）

---

## 1. 目的

miaoseq を「R から basecall まで全部やる巨大ラッパー」から、**ユーザーが用意した FASTQ を起点に、demux → サンプル別アンプリコン再構成 →（任意）editcall** するパッケージへ再設計する。

```
[ユーザー] Dorado などで basecall → FASTQ 準備
    ↓
[miaoseq]  doDemultiplex          … リードをサンプルへ分類（主出力: assignments）
    ↓ 任意
[miaoseq]  splitDemultiplexReads  … サンプル別 FASTQ 生成
    ↓
[miaoseq]  サンプルごとにクラスタリング / コンセンサス → アンプリコン配列
    ↓ 任意
[miaoseq]  editcall / report
```

---

## 2. 現行パイプラインの問題

現行 `miaoEditcall()`:

```
doBasecall → doDemultiplex(BLAST) → doAlign(BLAST) → doEditcall
                 ↑ 並列に basecall FASTA を消費
```

| 問題 | 説明 |
|------|------|
| basecall を R が抱え込む | Dorado は環境依存が強く、再実行・再開・ログ管理も重い。R 外で行う方が自然 |
| demux と amplicon assign が独立 | Step 4 でやっと結合。サンプル単位のコンセンサス導線になっていない |
| demux 出力が BLAST 風巨大表 | 下流が必要とするのは「read → sample」が主 |
| FASTA 中心 | コンセンサス / クラスタリングには quality 付き FASTQ の方が適切 |
| サンプル別配列がない | 下流ツールや HPC ジョブに渡しにくい |

---

## 3. 新パイプライン概要

### 3.1 ステップ構成

| Step | 関数 | 入力 | 主出力 | 備考 |
|------|------|------|--------|------|
| 0（R外） | Dorado 等 | pod5 | FASTQ | **パッケージ対象外** |
| 1 | `doDemultiplex` | FASTQ, index_list | `assignments.tsv` 等 | edlib 2段階 |
| 1b（任意） | `splitDemultiplexReads` | FASTQ + assignments | `by_sample/*.fq.gz` | demux 内からも呼べる |
| 2 | `doAmpliconResolve`（仮称） | assignments (+ FASTQ or by_sample) | サンプル別 consensus / clusters | **新規中核** |
| 3（任意） | `doEditcall` 系 | consensus / clusters + PAM 情報 | genotype summary | 既存を再設計 |
| — | `evalMiao` / `editViewer` | 出力 dir | レポート / Viewer | 出力仕様に合わせて更新 |

### 3.2 廃止するもの

| 対象 | 扱い |
|------|------|
| `R/basecall.R` | **削除** |
| `doBasecall()` | **廃止** |
| `makeMMI()` | **廃止**（必要ならユーザーが minimap2 を直接実行） |
| `miaoEditcall()` 内の basecall 呼び出し | **削除**。エントリポイントは再設計 |
| demux の `blast_path` | **削除** |
| 旧 `demultiplex_list.csv` 列互換 | **維持しない**（パイプライン再構築前提） |

BLAST ユーティリティ（`blast_utils.R`）は、Step 2 以降でまだ使う場合のみ残す。demux からは外す。

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
    {sample_id}/
      clusters.fasta              # または代表配列
      consensus.fasta
      cluster_counts.tsv
      stats.tsv
  editcall/                       # 任意 Step 3
    editcall_summary.csv
    ...
  report/                         # 任意
    ...
```

入力 FASTQ は `{out_dir}` 外でもよい（ユーザーがパスを渡す）。

---

## 5. Step 1: Demultiplex

詳細は [demux_revise.md](demux_revise.md)。

### 契約

```r
# 入力
fastq        # character: 1つ以上の FASTQ パス
demult_dir
index_list
sample_list = NULL
split_reads = FALSE

# 出力（常時）
assignments.tsv
summary_by_sample.tsv
unassigned.tsv

# 出力（split_reads=TRUE または後から splitDemultiplexReads）
by_sample/{sample_id}.fq.gz
```

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

## 6. Step 2: サンプル別アンプリコン再構成（新規）

現行 `doAlign`（全リード × amplicon DB BLAST → 表結合）を、**サンプル単位の配列再構成**に置き換える（または後段に移す）。

### 6.1 目的

demux 済みリードから、サンプルごとに:

1. （任意）アンプリコン種の分類（マルチプレックス遺伝子がある場合）
2. クラスタリングおよび/またはコンセンサス配列の取得
3. 代表配列・クラスタ頻度の出力

### 6.2 API 案

```r
doAmpliconResolve <- function(
  assignments,          # path or data.frame
  fastq = NULL,         # 元 FASTQ（by_sample が無いとき必須）
  sample_fastq_dir = NULL,  # demultiplex/by_sample
  out_dir,
  method = c("cluster", "consensus", "both"),
  amplicon_fn = NULL,   # 参考用 ref（任意）
  primer_list = NULL,   # 遺伝子分離用（任意）
  min_reads = 5,
  n_core = 1
)
```

### 6.3 処理方針（初期案・要議論）

| 選択肢 | 内容 | 用途 |
|--------|------|------|
| A. クラスタリング | isONclust / vsearch / custom OTU 風 | 編集多様性の可視化 |
| B. コンセンサス | 多数派やエラー補正付き合意配列 | 代表 genotype 配列 |
| C. ref-aware | amplicon DB / primer で遺伝子割当後に A/B | 現行 multi-gene 実験との接続 |

初期実装候補:

1. サンプル単位でリード抽出
2. （primer があれば）遺伝子ごとに分割
3. クラスタリング → 上位クラスタのコンセンサス
4. `cluster_counts.tsv` と `consensus.fasta` を出力

アルゴリズム選定（medaka / spoa / DECIPHER / 自前）は別イテレーションで決める。本ドキュメントでは **入出力契約を先に固定**する。

### 6.4 出力契約

```
amplicon/{sample_id}/
  consensus.fasta         # >cluster_id or gene_id
  clusters.fasta          # 任意: 代表 or 全クラスタ代表
  cluster_counts.tsv      # sample_id, cluster_id, gene_id?, n_reads, fraction
  stats.tsv               # n_reads_in, n_used, n_clusters, ...
```

---

## 7. Step 3: Editcall（再設計方針）

現行 `doEditcall` は demult 表 × align 表の結合前提。新フローでは:

```
consensus / cluster 配列
  → PAM / cut site 周辺の変異型判定
  → sample × gene の genotype summary
```

| 旧依存 | 新 |
|--------|----|
| `demult_out` BLAST 列 | 不要（サンプル単位で既に分離済み） |
| `align_out` 全リード表 | `amplicon/` の consensus / cluster 表に置換 |
| `strict`（barcode 端と aln 端の隣接チェック） | demux 側の end-trim / amplicon resolve 側で再定義 |

実装は Step 2 の出力が固まった後に着手。

---

## 8. エントリポイント

### 8.1 旧 `miaoEditcall()`

- basecall 依存を削除したうえで、**薄いオーケストレータ**として再定義する案  
- または **非推奨化**し、ユーザーが各 `do*` を順に呼ぶスタイルへ移行

推奨（明確さ優先）:

```r
# ユーザーが明示的にステップを呼ぶ
 dem <- doDemultiplex(fastq, demult_dir, index_list, sample_list, split_reads = FALSE)
 # 必要なら
 splitDemultiplexReads(fastq, dem$assignments, file.path(demult_dir, "by_sample"))
 amp <- doAmpliconResolve(dem$assignments, fastq = fastq, out_dir = amplicon_dir, ...)
 # 任意
 edit <- doEditcall(amp, pam_list, ...)
```

高水準ラッパーが必要なら後から `runMiaoseq()` のような新名で追加（旧名の破壊的再利用を避ける）。

### 8.2 `prepAmpliconDB`

参考配列構築用として **残す**（Step 2 の ref-aware モードや editcall 用）。  
入力はゲノム + primer/PAM。basecall とは独立。

---

## 9. パッケージ境界の再定義

### miaoseq がやる / やらない

| やる | やらない |
|------|----------|
| demux（edlib） | Dorado basecall |
| サンプル別 FASTQ 分割 | pod5 管理 |
| サンプル別クラスタ / コンセンサス | minimap2 index 構築の必須化 |
| editcall / report（任意） | 実験室 LIMS |

### README に明記する前提

```bash
# 例: ユーザーが先に basecall
dorado duplex sup pod5/ --min-qscore 10 > basecall.bam
samtools fastq basecall.bam > reads.fq
# 任意: 長さフィルタ等

# その後 R
doDemultiplex(fastq = "reads.fq", ...)
```

---

## 10. ファイル構成（目標）

```
R/
  demultiplex.R          # doDemultiplex, splitDemultiplexReads, 内部 edlib 関数
  amplicon_resolve.R     # doAmpliconResolve（新規）
  prep_amplicon_db.R     # 残す
  editcall.R             # 再設計
  report.R               # 出力仕様に合わせて更新
  pipeline.R             # 薄いラッパー or 簡素化
  blast_utils.R          # 必要なら残す（demux からは外す）

削除予定:
  R/basecall.R
  R/amplicon_assign.R    # doAlign を resolve に置換したら削除 or 大幅縮小
```

NAMESPACE:

- `export(doDemultiplex)`
- `export(splitDemultiplexReads)`
- `export(doAmpliconResolve)`（実装時）
- `doBasecall`, `makeMMI` を除去
- `doAlign` / `miaoEditcall` は移行期間の扱いは別途決定

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
  samples   = character(),           # 処理した sample_id
  out_dir   = character(1),
  table     = data.frame()           # cluster_counts を縦結合したもの等
)
```

---

## 12. 実装フェーズ（パイプライン全体）

### Phase A — Demux 刷新（本イテレーションの主対象）

- [ ] [demux_revise.md](demux_revise.md) に沿って `doDemultiplex` 実装
- [ ] `splitDemultiplexReads` 実装 + demux からの呼び出し
- [ ] 旧 BLAST demux 除去
- [ ] ドキュメント更新

### Phase B — basecall 撤去

- [ ] `R/basecall.R` 削除
- [ ] NAMESPACE / man / README / DESCRIPTION から basecall 記述を除去
- [ ] `pipeline.R` から `doBasecall` 呼び出し削除
- [ ] サンプルスクリプトを「ユーザー提供 FASTQ」前提に変更

### Phase C — Amplicon resolve

- [ ] 入出力契約の確定（本ドキュメント Section 6）
- [ ] プロトタイプ（1 サンプル・単純多数決コンセンサスでも可）
- [ ] マルチ gene 対応
- [ ] `doAlign` 経路の置換 / 廃止

### Phase D — Editcall / Report

- [ ] 新出力に合わせて `doEditcall` 再設計
- [ ] `evalMiao` / `editViewer` 更新

---

## 13. 成功基準（全体）

1. R から Dorado を呼ばなくても demux →（任意）split → amplicon resolve まで通る
2. demux 主出力が assignments 表であり、split は任意・後追い可能
3. 下流がサンプル単位でコンセンサス / クラスタを出せる
4. 旧「全リード BLAST → 後で join」への依存を解消できる道筋がある
5. README 上の前提（ユーザー用意 FASTQ）が一貫している

---

## 14. 未決事項（次の合意ポイント）

| 項目 | 選択肢 | メモ |
|------|--------|------|
| クラスタリング実装 | vsearch / isONclust / R 内 | 依存の重さ vs 精度 |
| コンセンサス実装 | spoa / medaka / DECIPHER / 自前 | quality 利用の有無 |
| マルチ gene 分離 | primer edlib / BLAST / なし | 現行 primer_list の再利用 |
| 旧 `miaoEditcall` | 改修 / 非推奨 / 新名 | 破壊的変更の許容度 |
| `doAlign` の扱い | 即削除 / deprecated 残置 | Phase C まで一時共存も可 |

これらは Phase C 着手前に決める。Phase A/B は demux + basecall 撤去を先行してよい。
