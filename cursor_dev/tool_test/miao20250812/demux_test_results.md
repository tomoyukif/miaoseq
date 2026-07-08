# Demux test results: miao20250812

**日付:** 2026-07-08  
**実装:** C++ edlib core（suffix HW → dict barcode；rescue 撤去）  
**入力:** `data_check.md` の FASTQ / index_list

---

## 実装メモ

ホットパスは **Rcpp + vendored edlib**（`edlibR` なし）:

- suffix: edlib HW LOC（front/rear × F/R、各 fwd+RC → 8 HW/read）
- barcode: mutant `unordered_map` lookup（offset ±0..4）
- orientation: F@front+R@rear と F@rear+R@front を評価し、有効ペア優先
- 既定 `end_window = 120`

**barcode_rescue / Top-K salvage は撤去**（合意を守ると recall がほぼ増えず、コストだけ増えたため）。

---

## 20k サブセット（C++ ・基準）

| 項目 | 値 |
|------|-----|
| 入力 | `basecall_filt_20k.fq`（20,000） |
| パラメータ | `end_window=120`, `n_core=16` |
| 出力 | `demultiplex_20k_cpp/` |
| 時間 | **~1.3 s** |
| edlib/read | **8.0** |
| assigned | **10,433**（52.2%） |
| vs BLAST agree | **96.8%** |

---

## フル FASTQ（C++ ・基準）

| 項目 | 値 |
|------|-----|
| 入力 | `basecall_filt.fq`（**1,224,546** reads） |
| 出力 | `demultiplex_full/` |
| パラメータ | `n_core=32`, `end_window=120`, `chunk_size=100000` |
| 時間 | **23.3 s** |
| assigned | **617,752**（50.5%） |
| vs BLAST agree | **96.3%** |

### 採用方針

ablation（ideas 1–5）は 20k で単独評価したが、精度・コストのトレードオフが本線採用に足りず **採用せず**。  
demux 本線は上記 **baseline（dict + suffix HW）固定**。詳細は `ablation_results.md`。
