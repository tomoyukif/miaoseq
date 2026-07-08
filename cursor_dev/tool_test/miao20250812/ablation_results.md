# Demux ablation results (20k)

Date: 2026-07-09 00:11  
FASTQ: `cursor_dev/tool_test/miao20250812/basecall_filt_20k.fq`  
n_core=16, end_window=120  
Oracle: BLAST `demultiplex_list.csv` (NA-free, unique `sseqid`)

| run | assign_rate | BLAST_both | pair_agree | new_agree | n_new | elapsed | edlib/read |
|-----|-------------|------------|------------|-----------|-------|---------|------------|
| baseline | 0.5216 | 8277 | 0.9675 | — | — | 1.07 | 8.0 |
| pair_constrained | 0.5487 | 8197 | 0.9584 | 0.235 | 1046 | 1.05 | 8.0 |
| sandwich_rescue | 0.5754 | 8458 | 0.9557 | 0.409 | 1074 | 1.63 | 255.5 |
| pass_b_full_index | 0.5933 | 8917 | 0.9596 | 0.858 | 1432 | 7.61 | 1476.9 |
| one_end_constrained | 0.5973 | 8664 | 0.9504 | 0.584 | 1513 | 1.04 | 33.0 |
| confidence_gate | 0.5216 | 8277 | 0.9675 | — | 0 | 0.92 | 8.0 |

### confidence_gate high-tier only

- n_high: 6649 (assign_rate_high=0.3325)
- BLAST_both_high: 6354  pair_agree_high: **0.9997**

### 合格目安との対比（計画 §5）

| 条件 | 結果 |
|------|------|
| `pair_agree` ≥ baseline − 0.5 pt（≥0.9625） | いずれも僅差で下回る（最良 pass_b 0.9596） |
| assign_rate / BLAST_both 増 | sandwich / pass_b / one_end で増加；pair_constrained は both がやや減 |
| `new_vs_baseline_agree` ≥ 90% | **どれも未達**。最良は pass_b **85.8%** |

### 読み取り

- **confidence_gate**: 割当不変。`high`（両端 edit0）だけ見ると BLAST 合意ほぼ完璧 → 後段フィルタとして有用。
- **pass_b_full_index**: recall / BLAST_both 最大、新規割当の合意も最良だが **~7× 時間・edlib 激増**。精度だけ見れば第1候補だがコスト重い。
- **one_end_constrained / sandwich_rescue**: assign は上がるが新規合意が弱い（60% / 41%）。単独採用は危険。
- **pair_constrained**: assign 微増だが新規合意が特に悪い（23%）。向き置換としては精度劣化寄り。

フル再確認は、計画どおり有望候補が出た場合に限る。現状は **gate を本線ラベルに足す**のが安全；recall 増は pass_b の閾値強化（margin / max_edit）を要検討。

Raw metrics: `cursor_dev/tool_test/miao20250812/ablation_metrics.tsv`  
Runner: `cursor_dev/tool_test/miao20250812/run_demux_ablation.R`
