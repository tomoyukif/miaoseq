# miaoseq

**Version 0.4.0** — see [inst/NEWS.md](inst/NEWS.md) for the changelog.

**miaoseq** is an R package for multiplexed indexed amplicon ONT sequencing (MIAOseq), focused on analysing CRISPR-Cas9 editing outcomes from Oxford Nanopore data.

The workflow starts from **user-provided FASTQ** (basecall outside R), then:

1. **Demultiplex** — assign reads to index pairs (`doDemultiplex`, edlib)
2. **Gene assignment** — assign/orient/filter reads (`doAssignGenes`)
3. Then either:
   - **Editcall (primary)** — per-read PAM-window alleles (`doEditcall`)
   - **Assemble (optional)** — clustering + consensus (`doAssembleAmplicons`) for representative amplicons / QC

## Overview

```text
[user] Dorado / other basecaller → FASTQ
   ↓
doDemultiplex(...)          → demultiplex/assignments.tsv
   ↓  optional
splitDemultiplexReads(...)  → demultiplex/by_sample/*.fq.gz
   ↓
doAssignGenes(...)          → amplicon_assign/gene_assignments.tsv
   ↓                            (length_outlier labeled when amplicon_fn set; not dropped)
   ├── Path B (editcall): doEditcall(gene_assign, ...) → editcall/
   └── Path A (assemble): doAssembleAmplicons(...) → amplicon/{index_pair_id}/
                              optional: doReassessAssemblies(...) → amplicon_reassess/
```

## Prerequisites

- **R ≥ 4.0** with Bioconductor packages listed in `DESCRIPTION` (`Biostrings`, `dplyr`, `GenomicRanges`, …)
- **External tools on `PATH`:** `vsearch` (clustering) and `abpoa` (consensus) for Assemble; optional `blastn` / `mmseqs` for Reassess. See `cursor_dev/containers/` for Apptainer wrappers.

## Installation

```r
devtools::install_github("tomoyukif/miaoseq")
library(miaoseq)
```

## Recommended usage

### Inputs

Prepare a FASTQ yourself, e.g. after Dorado basecalling and optional filtering.

```r
library(miaoseq)

fastq <- "/path/to/basecall_filt.fq"
out_dir <- "/path/to/run"
genome_fn <- "/path/to/genome.fa"
demult_dir <- file.path(out_dir, "demultiplex")
amplicon_dir <- file.path(out_dir, "amplicon")

index_list <- system.file("extdata", "index_list.csv", package = "miaoseq")
primer_list <- system.file("extdata", "amplicon_primers.csv", package = "miaoseq")
pam_list <- system.file("extdata", "agr8_pam_list.csv", package = "miaoseq")
```

### Step 0 — Amplicon reference (optional but recommended for editcall)

```r
amplicon_fn <- prepAmpliconDB(
  primer_list = primer_list,
  genome_fn = genome_fn,
  out_dir = out_dir,
  expected_length = 500L
)
# Or supply curated amplicons (names must match primer gene IDs):
# amplicon_fn <- prepAmpliconDB(
#   primer_list = primer_list,
#   out_dir = out_dir,
#   amplicon_fasta = "/path/to/amplicons.fa"
# )
# → {out_dir}/ref/amplicon.fa
```

Uses exact primer matching (`Biostrings::matchPattern`) against the reference genome. Stops on multi-locus hits or spans longer than `2 * expected_length`. No BLAST required.

### Step 1 — Demultiplex

Prefer a **1:1** index layout (each F/R ID used once) or barcodes with Hamming distance ≥ `2 * max_barcode_edit + 1`. Combinatorial plates (shared F/R IDs across wells) are allowed, but risk of a false dual hit landing in another legal well is on the user — on such layouts, a wrong dual assignment can still map to a valid well. `sample_list` `sample_name` values may be duplicated; demux FASTQ / summaries are always one row or file per `index_pair_id`.

```r
dem <- doDemultiplex(
  fastq = fastq,
  demult_dir = demult_dir,
  index_list = index_list,
  n_core = 8,
  split_reads = FALSE,
  allow_single_end = FALSE,
  return_tables = FALSE   # keep RAM low; tables stay on disk
)
# Optional later split (also C++ streaming):
# splitDemultiplexReads(fastq, file.path(demult_dir, "assignments.tsv"),
#                       file.path(demult_dir, "by_sample"), compress = TRUE)
```

### Step 2 — Gene assignment

```r
ga <- doAssignGenes(
  assignments = file.path(demult_dir, "assignments.tsv"),
  out_dir = file.path(out_dir, "amplicon_assign"),
  fastq = fastq,
  primer_list = primer_list,
  amplicon_fn = amplicon_fn,   # optional; labels trimmed-insert length outliers
  length_tolerance = 0.25,     # ±25% of (amplicon.fa width − F/R primer)
  n_core = 8,
  # store_sequences = FALSE,   # smaller TSV for Editcall-only (Assemble needs TRUE)
  overwrite = TRUE,
  stats_unassign = TRUE
)
# → amplicon_assign/gene_assignments.tsv
# → amplicon_assign/stats_unassigned.tsv  (when stats_unassign = TRUE)
```

When `amplicon_fn` is set, primer-assigned reads whose **trimmed insert** length
falls outside `expected ± length_tolerance` (or fails primer trim) are **kept**
with `assign_status = "length_outlier"`. Expected length is
`width(amplicon.fa) − F primer − R primer`. Primer misses (`no_primer_hit`) stay
in the table. `doEditcall` / `doAssembleAmplicons` process both `assigned` and
`length_outlier` by default so large excision products are not dropped upstream.
Use `n_core` for primer assignment and length-trim. Set `store_sequences = FALSE`
to omit oriented sequences from the TSV when only Editcall will follow (Assemble
/ Reassess reading from disk still need the sequences).

### Step 3 — Edit-calling (Path B; primary)

```r
edit <- doEditcall(
  gene_assign = ga,
  pam_list = pam_list,
  genome_fn = genome_fn,
  amplicon_fn = amplicon_fn,
  primer_list = primer_list,
  editcall_dir = file.path(out_dir, "editcall"),
  fastq = fastq,
  check_window = 10,
  anchor_bp = 5,
  max_expand = 50,
  min_count = 5,
  n_core = 8
)
```

Outputs: `editcall_all.csv`, `editcall_filtered.csv`, `editcall_summary.csv`, `intact_seq.fa`,
`run_stats.txt` (plus `editcall_joint*.csv` when a gene has multiple PAMs). Alleles and
summaries are keyed by `index_pair_id`. `sample_id` is a reference label only
(from `sample_list` / demultiplex `assignments.tsv`). Editcall keeps gene primers
in the aligned span (only outer adapters are removed) and places PAM windows on
the full `amplicon.fa` coordinates. After `min_count`, alleles with
`count > top/2` in each `index_pair_id × target_gene` are kept. `intact_seq.fa`
is the fixed genome ±`check_window` WT slice; compare alleles to CSV `ref_seq`,
not that FA.

### Step 3′ — Amplicon assembly (Path A; optional)

```r
amp <- doAssembleAmplicons(
  gene_assign = ga,
  out_dir = amplicon_dir,
  primer_list = primer_list,
  cluster_backend = "vsearch",
  min_cluster_identity = 0.95,  # vsearch --id; accepts (0, 1] (no 0.5 floor)
  min_reads = 5,
  min_cluster_reads = 5,
  # strict_end_trim = TRUE,  # default: cluster trimmed inserts only
  # min_cluster_purity = 0.8,  # default NULL (disabled); set for strict filtering
  # assembly_backend = "overlap_graph",  # optional graph-based modal paths
  overwrite = TRUE
)
```

Per `index_pair_id`, assemble writes `consensus.fasta`, `clusters.fasta`,
`cluster_counts.tsv`, `stats.tsv`, and `unassigned_to_cluster.tsv`
(U1 membership for Reassess).

### Step 2c — Reassess assemblies (Path A diagnostic; optional)

```r
re <- doReassessAssemblies(
  amplicon_out = amp,
  gene_assign = ga,
  out_dir = file.path(out_dir, "amplicon_reassess"),
  primer_list = primer_list,
  # backend = "edlib",  # default; independent edit-distance QC of vsearch clustering
  consensus_merge_max_edit = 12L,  # Q1 highlight threshold only (no auto-merge)
  read_assign_max_edit = 12L,
  # min_identity = NULL,  # if unset: 1 - edit/median_len, clamped to [0.50, 0.99]
  overwrite = TRUE
)
# → summary_by_sample.tsv (long; one row per index_pair_id×backend)
# → summary_compare.tsv (wide; reassignable_frac / near-dup groups side-by-side)
# → {backend}/{index_pair_id}/consensus_pairwise.tsv
# → {backend}/{index_pair_id}/consensus_near_duplicate_groups.tsv  (review only)
# → {backend}/{index_pair_id}/unassigned_to_consensus.tsv
# Does not modify amplicon/ and does not write merged sequences.
```

### Reporting (optional)

```r
evalMiao(out_dir, output_reads = FALSE)
# → miao_summary/read_stats.tsv, indexed_reads_per_index.tsv, gene_reads_per_gene.tsv

editViewer(out_dir, sample_list = "/path/to/sample_list.csv")
# → editviewer/edit_viewer_plate*.pdf
# Dual-guide / dual-cut: also open editcall/editcall_joint_summary.csv
# Uniq/Dup on the PDF compares wells that share sample_name (not index_pair_id).
# In-frame colors use |n_ins - n_del| mod 3 only (not CDS phase).
```

`sample_list` is a CSV with **plate well coordinates** (header optional), not barcode/index sequences:

```text
index_pair_id,sample_name,plate_id,row_id,col_id
miaoBC0001,sample01,plate1,A,1
miaoBC0009,sample02,plate1,A,2
```

Duplicate `sample_name` values are allowed. Downstream processing never
aggregates by `sample_name` / `sample_id`; the well key is always
`index_pair_id`. `sample_id` appears on editcall allele tables only as a
reference label.

- `row_id`: `A`–`H`
- `col_id`: `1`–`12`

Do **not** pass `index_list.csv` (barcode file) as `sample_list`.

## Input file formats

#### Primer list (`primer_list`)
CSV, two columns: primer ID (`*_F` / `*_R`), sequence. Header row optional.
Required for `doAssignGenes()` and `doAssembleAmplicons()` (`NULL` / no-ref mode is not supported).

#### Index list (`index_list`)
CSV, five columns: index pair ID, F index ID, F sequence, R index ID, R sequence. Header row optional. At least two unique F and two unique R index sequences are required.
On combinatorial layouts, a spurious dual barcode call may still match another plate well; prefer 1:1 pairing or high inter-barcode distance when cross-talk must be minimized.

#### PAM list (`pam_list`)
CSV (header optional):

1. gene ID matching primer / `amplicon.fa` names
2. chromosome / seqname matching genome FASTA headers **exactly** (no automatic `chr` / zero-pad)
3. PAM start (1-based)
4. optional strand `+` / `-` (Cas9: `+` → cut = pam−3, `-` → cut = pam+3 on genome; missing → cut = pam start)
5. optional guide ID (**required when a gene has multiple rows**)

A locus with **one PAM row** is labelled by gene only (`target_gene = SD1`), even if column 5 is filled. Multiple PAM rows on the same gene use `target_gene = gene_guide` (e.g. `SD1_g1`). Bundled `agr8_pam_list.csv` includes Cas9 strand (column 4) and is single-guide per gene.
Multiple guides on one amplicon: same gene ID + distinct guide column — guide-level alleles (Plan A) plus dual-cut excision (`editcall_joint*.csv`, Plan A′). Plan A may count an excision read on **both** guides; that double counting is intentional. With three or more guides, if both adjacent pairs are excision the middle guide Plan A allele is the token `---`. For pair-level excision rates read `editcall_joint_summary.csv` (editViewer does not load joint files). See `cursor_dev/pipeline_revise.md` §7.3–7.4 / §20.

#### Assemble column note (`fraction_bucket`)
In `cluster_counts.tsv`, `fraction_bucket` is the fraction of reads in the **`index_pair_id` × gene bucket after primer trim (cluster input)**, not of all demultiplexed reads for the well (formerly `fraction_sample`). `fraction` is within that gene after clustering. Minor clusters below `min_cluster_reads` (default 5) are excluded; `max_clusters` defaults to `Inf` (no top-N truncation).

#### Clustering identity note
Assemble clustering uses **vsearch** only (`--cluster_fast`, `--iddef 2`) with `--id = min_cluster_identity` (default **0.95**, range **(0, 1]** — values below 0.5 are allowed for aggressive merging sweeps). Consensus is **abpoa** (FASTA, no quality weighting, no racon polish). `clusters.fasta` holds member reads (`>{read_id} {cluster_id}`).

Rough tuning guide for `min_cluster_identity` (ONT ~1–2 kb inserts; exact numbers vary by dataset):

| Setting | Typical effect |
|---------|----------------|
| **~0.4–0.6** | Strong merging — often ~1 retained cluster per index_pair×gene, with most reads assigned. Good when a single dominant template is expected (e.g. clonal / colony DNA). |
| **~0.8–0.9** | Often a practical balance: high assignment rate, a few secondary clusters allowed. |
| **0.95 (default)** | More separation; more reads may fall below `min_cluster_reads` as tiny clusters. |
| **≥0.98–0.99** | Very strict — many ONT-error shards fail the minor-cluster filter, so assignment rate can collapse even when true diversity is low. |

Lower identity raises assigned fraction by merging noise into larger clusters (exact-sequence purity within a cluster usually drops). Higher identity splits clusters; retained count may peak then fall once shards drop below `min_cluster_reads`. Prefer Editcall for local allele frequencies around known cuts; use Assemble identity to control how aggressively full-length representatives are merged.

#### Pathway roles
**Assemble** restores full-length amplicon representatives (within-cluster variation is collapsed). **Editcall** estimates local edit patterns around known cut sites. For diverse loci such as 16S, use Assemble for global composition; do not treat Editcall local windows as a substitute for full-insert clustering. Residual uncertainty in ONT amplicon consensus is shared across tools and has no complete solution.

#### Reassess backend selection
Default `backend = "edlib"` provides an independent exact edit-distance evaluation of clustering results — a different perspective from vsearch's %identity used during clustering. Additional backends (`"vsearch"`, `"blastn"`, `"mmseqs"`) can be added for cross-metric comparison (`summary_compare.tsv`). For identity-based backends, when `min_identity` is `NULL` it is derived as `1 - edit / median_consensus_length` and clamped to **[0.50, 0.99]** (separate from Assemble's pass-through identity). Use `fraction` (not `fraction_bucket`) when ranking cluster importance within a gene.

## Output layout

```text
{out_dir}/
  ref/
    amplicon.fa                 # prepAmpliconDB
  demultiplex/
    assignments.tsv
    unassigned.tsv
    summary_by_sample.tsv
    stats_unassigned.tsv        # optional (stats_unassign = TRUE)
    by_sample/                  # optional; one FASTQ per index_pair_id
  amplicon_assign/
    gene_assignments.tsv        # keyed by index_pair_id
    stats_unassigned.tsv        # optional
  amplicon/                     # doAssembleAmplicons
    summary_by_sample.tsv
    run_stats.txt
    {index_pair_id}/
      consensus.fasta
      clusters.fasta
      cluster_counts.tsv
      stats.tsv
      unassigned_to_cluster.tsv
  amplicon_reassess/            # doReassessAssemblies (optional)
    summary_by_sample.tsv
    summary_compare.tsv
    {backend}/{index_pair_id}/...
  editcall/
    editcall_all.csv
    editcall_filtered.csv
    editcall_summary.csv
    editcall_joint*.csv         # multi-guide / dual-cut (when applicable)
    intact_seq.fa
    run_stats.txt
  miao_summary/                 # evalMiao
  editviewer/                   # editViewer PDFs
```

## Troubleshooting

1. **Demux assigned ≫ processable gene-assign rows** — with `amplicon_fn`, length outliers stay in `gene_assignments.tsv` as `assign_status = "length_outlier"` (not dropped, not counted as `no_primer_hit`). Downstream Editcall / Assemble use `assigned` + `length_outlier` + `trim_fail` by default. Check `stats_unassigned.tsv` for primer-miss reasons only. `trim_fail` is primer-span trim failure; `length_outlier` is a successful trim with unexpected insert length.
2. **Empty amplicon clusters** — check `skip_reason` in `stats.tsv` (`low_gene_reads`, `no_gene_assigned`, `no_clusters`).
3. **Slow resolve without `by_sample`** — provide `sample_fastq_dir`, or accept a single FASTQ pass.
4. **Tool binaries not found** — Assemble requires `vsearch` and `abpoa` on `PATH`. Use `cursor_dev/containers/pull_sifs.sh` to pull Apptainer images and add `cursor_dev/containers/bin/` to `PATH`. Reassess optionally uses `blastn` / `mmseqs`.

## Further documentation

- `cursor_dev/pipeline_revise.md` — pipeline redesign (§20 review lock, §21 speed, §22 `index_pair_id` grouping)
- `cursor_dev/demux_revise.md` — demultiplex design
- `cursor_dev/ampliconresolve_plan.md` — amplicon resolve phases
- `cursor_dev/editcall_adapter.md` — Editcall connection / API
