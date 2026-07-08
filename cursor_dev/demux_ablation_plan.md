# Demux feature ablation 計画（ideas 1–5）

作成: 2026-07-09  
前提: 現行 C++ dict 本線を **baseline** とし、提案 1–5 を **独立実験フラグ** として実装・単独評価する。

---

## 1. 決定事項


| 項目         | 決定                                                                       |
| ---------- | ------------------------------------------------------------------------ |
| デフォルト      | **全て OFF**（現状どおり）                                                        |
| 正解（oracle） | 旧 BLAST `demultiplex_list.csv`（`sseqid` ↔ `index_pair_id`、NA 除外・read 一意） |
| 評価単位       | **単独 ON**（組合せは第2段）                                                       |
| 主評価データ     | 20k → 有望なものだけフル                                                          |


### 指標（BLAST 正例集合 P、手法割当 A、全リード N）


| 指標                        | 定義                                       | 見るもの               |
| ------------------------- | ---------------------------------------- | ------------------ |
| assign_rate               | |A|/N                                    | 粗 recall           |
| BLAST_both                | |A \cap P|                               | BLAST と重なる割当数      |
| pair_agree                | both のうち pair ID 一致率                     | precision-like（合意） |
| false_vs_BLAST            | both のうち不一致                              | 悪化                 |
| **new_vs_baseline_agree** | baseline 未割当 → 機能 ON で新規に入った分の BLAST 一致率 | **その機能だけの精度**      |
| elapsed / n_edlib         | 時間・コスト                                   | 実用性                |


Idea 5 は割当集合を変えないので、`A_all` と `A_high`（tier=high のみ）の2系統で同じ指標を出す。

---



## 2. 独立性の設計

```text
Baseline (dict 本線)
  ├─ +pair_constrained      (idea 1)  … 向き判定を置換（別ラン比較）
  ├─ +sandwich_rescue       (idea 2)  … dict miss 端のみ
  ├─ +pass_b_full_index     (idea 3)  … unassigned のみ第2パス
  ├─ +one_end_constrained   (idea 4)  … 片端確定後の反対端のみ
  └─ +confidence_gate       (idea 5)  … 出力ラベルのみ（割当不変）
```


| Idea | Flag                  | 介入点              | 既存割当を書き換え？              |
| ---- | --------------------- | ---------------- | ----------------------- |
| 1    | `pair_constrained`    | 合法ペア空間で向き判定      | Yes（ロジック置換 → **別ラン比較**） |
| 2    | `sandwich_rescue`     | dict miss 端      | No                      |
| 3    | `pass_b_full_index`   | unassigned のみ    | No                      |
| 4    | `one_end_constrained` | 片端成功・反対端失敗時      | そこだけ                    |
| 5    | `confidence_gate`     | 出力 `assign_tier` | No                      |


**Idea 1 の注意:** 向き判定自体を変えるので「後処理」ではない。ablation では baseline と完全に別ランで比較する（標準的）。

---



## 3. 各機能の仕様（実装単位）



### Idea 1 — `pair_constrained`

- 各端 Top-K（K=3）dict EndHit
- F×R は `sample_map` **合法組だけ** スコア
- 採択: 最良が一意、かつ `second - best ≥ pair_margin`（default: barcode-edit 1 相当）
- さもなくば `ambiguous_pair` / unassigned



### Idea 2 — `sandwich_rescue`

- **dict 失敗端のみ**
- suffix 位置 → 上流 prefix HW → 挟み区間を親 barcode HW
- 条件: prefix OK + barcode 一意 + edit≤`max_barcode_edit` + 反対端と合法ペア



### Idea 3 — `pass_b_full_index`

- Pass A = baseline（確定割当は触らない）
- Pass B: unassigned のみ、両端で **フル index（prefix+barcode+suffix）** を合法 384 ペアだけ HW
- 一意 + margin で追加割当；`assign_method=pass_b`



### Idea 4 — `one_end_constrained`

- 片端が `barcode_edit==0` で一意 → 反対端はペア可能な **≤12 barcode** だけ探索
- 両端未確定 / 両端確定済みは no-op



### Idea 5 — `confidence_gate`

- アルゴリズム不変
- `assign_tier`: `high`（例: complete + 両端 edit0）/ `review`（その他 assigned）
- 評価は `A_all` と `A_high` の両方

---



## 4. 実装スコープ（ファイル）


| ファイル                                                         | 内容                                                    |
| ------------------------------------------------------------ | ----------------------------------------------------- |
| `src/demux_internal.h` / `demux_dict.cpp` / `demux_core.cpp` | フラグ分岐・各 idea の C++ 本体                                 |
| `R/demultiplex.R`                                            | 同名引数（default FALSE）、`assign_method` / `assign_tier` 列 |
| `R/demux_eval.R`（新規・internal）                                | BLAST 突合メトリクス                                         |
| `cursor_dev/tool_test/miao20250812/run_demux_ablation.R`     | 6 ラン実行                                                |
| `…/ablation_metrics.tsv` + `ablation_results.md`             | 結果表                                                   |


API イメージ:

```r
doDemultiplex(
  ...,
  pair_constrained = FALSE,
  sandwich_rescue = FALSE,
  pass_b_full_index = FALSE,
  one_end_constrained = FALSE,
  confidence_gate = FALSE,
  pair_margin = 1L,
  max_end_candidates = 3L
)
```

---



## 5. 実装・試験順

1. **Scaffold** — flags、列、eval harness、baseline 1 行
2. **Idea 5**（最軽量）
3. **Idea 4** → **Idea 2** → **Idea 1** → **Idea 3**（コスト増順）
4. **20k ablation 表**（baseline + 各単独）
5. 有望 1–2 個だけフル再確認



### 合格の目安（20k）

- `pair_agree` ≥ baseline − 0.5 pt（望ましくは ≥ baseline）
- `assign_rate` または `BLAST_both` が増加
- `new_vs_baseline_agree` ≥ 90%（弱いなら idea 5 の `review` 行き推奨）

---



## 6. Non-goals（このパス）

- 機能組合せの本格チューニング
- 以前撤去した緩い NW rescue の無条件復活
- `require_prefix=TRUE` をデフォルト化

---



## 7. 期待アウトプット

`ablation_results.md` に例えば:


| run                     | assign_rate | BLAST_both | pair_agree | new_agree | elapsed |
| ----------------------- | ----------- | ---------- | ---------- | --------- | ------- |
| baseline                | …           | …          | …          | —         | …       |
| +pair_constrained       | …           | …          | …          | …         | …       |
| +sandwich_rescue        | …           | …          | …          | …         | …       |
| +pass_b_full_index      | …           | …          | …          | …         | …       |
| +one_end_constrained    | …           | …          | …          | …         | …       |
| +confidence_gate (high) | …           | …          | …          | …         | …       |


これで「どれが precision を維持して recall を上げるか」を定量比較する。