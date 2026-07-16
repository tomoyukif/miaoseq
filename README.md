# miaoseq

**Version 0.2.0** — see [inst/NEWS.md](inst/NEWS.md) for the changelog.

**miaoseq** is an R package for multiplexed indexed amplicon ONT sequencing (MIAOseq), focused on analysing CRISPR-Cas9 editing outcomes from Oxford Nanopore data.

The workflow starts from **user-provided FASTQ** (basecall outside R), then:

1. **Demultiplex** — assign reads to samples (`doDemultiplex`, edlib)
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
   ↓
   ├── Path B (editcall): doEditcall(gene_assign, ...) → editcall/
   └── Path A (assemble): doAssembleAmplicons(...) → amplicon/{sample}/
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

Prefer a **1:1** index layout (each F/R ID used once) or barcodes with Hamming distance ≥ `2 * max_barcode_edit + 1`. Combinatorial plates (shared F/R IDs across wells) are allowed, but risk of a false dual hit landing in another legal well is on the user — on such layouts, a wrong dual assignment can still map to a valid well.

```r
dem <- doDemultiplex(
  fastq = fastq,
  demult_dir = demult_dir,
  index_list = index_list,
  n_core = 8,
  split_reads = FALSE,
  allow_single_end = FALSE
)
```

### Step 2 — Gene assignment

```r
ga <- doAssignGenes(
  assignments = file.path(demult_dir, "assignments.tsv"),
  out_dir = file.path(out_dir, "amplicon_assign"),
  fastq = fastq,
  primer_list = primer_list,
  amplicon_fn = amplicon_fn,
  overwrite = TRUE
)
```

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
  min_count = 5
)
```

Outputs: `editcall_all.csv`, `editcall_filtered.csv`, `editcall_summary.csv`, `intact_seq.fa`.

### Step 3′ — Amplicon assembly (Path A; optional)

```r
amp <- doAssembleAmplicons(
  gene_assign = ga,
  out_dir = amplicon_dir,
  primer_list = primer_list,
  cluster_backend = "vsearch",
  min_reads = 5,
  min_cluster_reads = 5,
  # min_cluster_purity = 0.8,  # default NULL (disabled); set for strict filtering
  # assembly_backend = "overlap_graph",  # optional graph-based modal paths
  overwrite = TRUE
)
```

Per sample, assemble also writes `unassigned_to_cluster.tsv` (U1 membership for Phase H).

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
  overwrite = TRUE
)
# → summary_by_sample.tsv (long; one row per sample×backend)
# → summary_compare.tsv (wide; reassignable_frac / near-dup groups side-by-side)
# → {backend}/{sample}/consensus_pairwise.tsv
# → {backend}/{sample}/consensus_near_duplicate_groups.tsv  (review only)
# → {backend}/{sample}/unassigned_to_consensus.tsv
# Does not modify amplicon/ and does not write merged sequences.
```

### Reporting (optional)

```r
evalMiao(out_dir, output_reads = FALSE)
# → miao_summary/read_stats.tsv, indexed_reads_per_index.tsv, gene_reads_per_gene.tsv

editViewer(out_dir, sample_list = "/path/to/sample_list.csv")
# → editviewer/edit_viewer_plate*.pdf
```

`sample_list` is a headerless CSV with **plate well coordinates**, not barcode/index sequences:

```text
index_pair_id,sample_name,plate_id,row_id,col_id
miaoBC0001,sample01,plate1,A,1
miaoBC0009,sample02,plate1,A,2
```

- `row_id`: `A`–`H`
- `col_id`: `1`–`12`

Do **not** pass `index_list.csv` (barcode file) as `sample_list`.

## Input file formats

#### Primer list (`primer_list`)
CSV, two columns: primer ID (`*_F` / `*_R`), sequence.

#### Index list (`index_list`)
CSV, five columns: index pair ID, F index ID, F sequence, R index ID, R sequence.
On combinatorial layouts, a spurious dual barcode call may still match another plate well; prefer 1:1 pairing or high inter-barcode distance when cross-talk must be minimized.

#### PAM list (`pam_list`)
CSV (no header): (1) gene ID matching primer / `amplicon.fa` names, (2) chromosome/seqname matching genome FASTA headers **exactly** (no automatic `chr` / zero-pad), (3) PAM start (1-based), (4) optional guide ID (required when a gene has multiple rows), (5) optional strand `+`/`-` (Cas9: `+` → cut = pam−3, `-` → cut = pam+3 on genome; missing → cut = pam start).
Multiple guides on one amplicon: same gene ID + distinct guide column — Guide-level alleles (Plan A) plus dual-cut excision (`editcall_joint*.csv`, Plan A′). See `cursor_dev/pipeline_revise.md` §7.3–7.4 / §18.

#### Assemble column note (`fraction_bucket`)
In `cluster_counts.tsv`, `fraction_bucket` is the fraction of reads in the **sample × gene bucket**, not of all demultiplexed reads for the sample (formerly `fraction_sample`). `fraction` is within that gene after clustering. Minor clusters below `min_cluster_reads` (default 5) are excluded; `max_clusters` defaults to `Inf` (no top-N truncation).

#### Clustering identity note
Assemble clustering uses **vsearch** only (`--cluster_fast`, `--iddef 2`) with `--id = min_cluster_identity` (default **0.95**). Consensus is **abpoa** (FASTA, no quality weighting, no racon polish). `clusters.fasta` holds member reads (`>{read_id} {cluster_id}`).

#### Pathway roles
**Assemble** restores full-length amplicon representatives (within-cluster variation is collapsed). **Editcall** estimates local edit patterns around known cut sites. For diverse loci such as 16S, use Assemble for global composition; do not treat Editcall local windows as a substitute for full-insert clustering. Residual uncertainty in ONT amplicon consensus is shared across tools and has no complete solution.

#### Reassess backend selection
Default `backend = "edlib"` provides an independent exact edit-distance evaluation of clustering results — a different perspective from vsearch's %identity used during clustering. Additional backends (`"vsearch"`, `"blastn"`, `"mmseqs"`) can be added for cross-metric comparison (`summary_compare.tsv`). Use `fraction` (not `fraction_bucket`) when ranking cluster importance within a gene.

## Output layout

```text
{out_dir}/
  ref/
    amplicon.fa
  demultiplex/
    assignments.tsv
    summary_by_sample.tsv
    by_sample/                 # optional
  amplicon/
    summary_by_sample.tsv
    gene_assignments.tsv
    {sample_id}/
      consensus.fasta
      cluster_counts.tsv
      stats.tsv
  editcall/
    editcall_summary.csv
    ...
  miao_summary/              # evalMiao
  editviewer/                # editViewer PDFs
```

## Troubleshooting

1. **Empty amplicon clusters** — check `skip_reason` in `stats.tsv` (`low_gene_reads`, `no_gene_assigned`, `no_clusters`).
2. **Slow resolve without `by_sample`** — provide `sample_fastq_dir`, or accept a single FASTQ pass.
3. **Tool binaries not found** — Assemble requires `vsearch` and `abpoa` on `PATH`. Use `cursor_dev/containers/pull_sifs.sh` to pull Apptainer images and add `cursor_dev/containers/bin/` to `PATH`. Reassess optionally uses `blastn` / `mmseqs`.

## Further documentation

- `cursor_dev/pipeline_revise.md` — pipeline redesign
- `cursor_dev/demux_revise.md` — demultiplex design
- `cursor_dev/ampliconresolve_plan.md` — amplicon resolve phases
