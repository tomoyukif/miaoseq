# テスト入力: miao20250812

**用途:** `devel` ブランチの demultiplex →（任意）split → amplicon resolve / editcall のテスト  
**実験:** `miao20250812`

---

## 必須入力

| 用途 | パス | 備考 |
|------|------|------|
| multiplexed FASTQ | `/home/ftom/01_wd/invitro_domestication/genomeEdit/edit_check/nanopore/20250918_pipeline_test/basecall/basecall_filt.fq` | 約 927 MB。size-filtered 済み。`doDemultiplex` の主入力 |
| index list | `/home/ftom/01_wd/invitro_domestication/genomeEdit/edit_check/nanopore/index_list.csv` | 384 ペア。ヘッダーなし 5 列。`inst/extdata/index_list.csv` と同内容 |

スモーク用 index（2 サンプル）: `inst/extdata/index_list_smoke.csv`

---

## 下流で使う入力（任意）

| 用途 | パス | 備考 |
|------|------|------|
| PCR primers | `/home/ftom/01_wd/invitro_domestication/genomeEdit/edit_check/nanopore/amplicon_primers.csv` | 8 アンプリコン × F/R。遺伝子分離・resolve 用。`inst/extdata/` と同内容 |
| PAM / cut sites | `/home/ftom/01_wd/invitro_domestication/genomeEdit/edit_check/nanopore/agr8_pam_list.csv` | editcall 用。`inst/extdata/` と同内容 |
| 参照ゲノム | `/home/ftom/01_wd/genomeData/rice/cultivar_sativa/nb_combined/version_2023/fasta/nb_genome.fa` | editcall / ref-aware 用。約 362 MB |

---

## 出力先（テスト用）

```
/home/ftom/01_wd/softDevel/miaoseq/cursor_dev/tool_test/miao20250812/
```

---

## 最小呼び出し例

```r
fastq      <- "/home/ftom/01_wd/invitro_domestication/genomeEdit/edit_check/nanopore/20250918_pipeline_test/basecall/basecall_filt.fq"
index_list <- "/home/ftom/01_wd/invitro_domestication/genomeEdit/edit_check/nanopore/index_list.csv"
demult_dir <- "/home/ftom/01_wd/softDevel/miaoseq/cursor_dev/tool_test/miao20250812/demultiplex"

# doDemultiplex(fastq = fastq, demult_dir = demult_dir, index_list = index_list, ...)
# splitDemultiplexReads(fastq = fastq, assignments = ..., out_dir = file.path(demult_dir, "by_sample"))
```
