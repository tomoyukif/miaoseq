# miaoseq

**miaoseq** is an R package for analyzing multiplexed indexed amplicon Oxford Nanopore sequencing (MIAOseq) data and characterizing CRISPR-Cas genome editing outcomes.

The package provides **two workflows**:

| Workflow | Entry point | Best for |
|----------|-------------|----------|
| **Legacy (read-level)** | `miaoEditcall()` | Pod5 → basecall → BLAST demux/align → per-read edit calling |
| **Refless (cluster-consensus)** | `miaoPipeline()` | Pre-basecalled FASTQ → edit-distance demux → clustering → consensus edit calling |

Both workflows share plate visualization via `editViewer()`.

---

## Table of contents

1. [System requirements](#system-requirements)
2. [Installation](#installation)
3. [Environment setup](#environment-setup)
   - [R dependencies](#r-dependencies)
   - [Legacy pipeline tools](#legacy-pipeline-tools)
   - [Refless pipeline tools](#refless-pipeline-tools)
   - [Apptainer image (recommended for refless)](#apptainer-image-recommended-for-refless)
   - [Environment variables](#environment-variables)
4. [Input file formats](#input-file-formats)
5. [Usage: refless pipeline](#usage-refless-pipeline)
6. [Usage: legacy pipeline](#usage-legacy-pipeline)
7. [Edit calling and visualization](#edit-calling-and-visualization)
8. [Output directories](#output-directories)
9. [Troubleshooting](#troubleshooting)
10. [Performance notes](#performance-notes)

---

## System requirements

- **OS**: Linux (tested on Ubuntu). macOS may work for R-only steps; external tools and Apptainer are Linux-oriented.
- **R**: ≥ 4.2 recommended.
- **RAM**: 16–32 GB for typical 384-well Flongle runs (legacy pipeline); refless pipeline is lighter per sample but clustering benefits from available memory.
- **CPU**: Multi-core recommended (`n_core` in legacy pipeline; mmseqs2/vsearch use multiple threads internally).

---

## Installation

```r
# From GitHub (refless branch or main):
devtools::install_github("tomoyukif/miaoseq", ref = "refless")

library(miaoseq)
```

For local development from a cloned repository:

```r
devtools::load_all()   # or devtools::install()
```

Release notes are in [`inst/NEWS.md`](inst/NEWS.md).

---

## Environment setup

### R dependencies

Installed automatically with the package:

| Package | Purpose |
|---------|---------|
| Biostrings | Sequence I/O and manipulation |
| pwalign | End-window alignment in `doDemultiplex2()` |
| GenomicRanges, IRanges | Genomic interval operations |
| dplyr | Data frame joins |
| BiocGenerics, methods, stats, tools, utils | Core utilities |

**Suggested** packages (not hard dependencies):

| Package | Purpose |
|---------|---------|
| ggplot2, tidyr | Used by `editViewer()` for plate heatmaps |
| devtools | Package installation from GitHub |

```r
install.packages(c("ggplot2", "tidyr", "devtools"))
```

### Legacy pipeline tools

Required only when running `miaoEditcall()`:

| Tool | Role | Installation |
|------|------|--------------|
| [Dorado](https://github.com/nanoporetech/dorado) | Basecalling pod5 → FASTQ | ONT releases |
| [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) | Demultiplexing and read alignment | `makeblastdb`, `blastn` |
| [samtools](https://github.com/samtools/samtools) | BAM processing | conda or system package |

Provide paths when calling `miaoEditcall()`:

```r
dorado_path  <- "/path/to/dorado"
blast_path   <- "/path/to/blast/bin"   # directory containing blastn
samtools_path <- "/path/to/samtools"
```

### Refless pipeline tools

Required when running `miaoPipeline()` and related functions:

| Tool | Role | Required for |
|------|------|--------------|
| [MMseqs2](https://github.com/soedinglab/MMseqs2) | Fast linclust (default) | `cluster_method = "mmseqs2"` |
| [minimap2](https://github.com/lh3/minimap2) | Cluster refinement; ref mapping in edit calling | `refine = TRUE`; `prepRefTargets()` |
| [vsearch](https://github.com/torognes/vsearch) | meshclust clustering; chimera detection | `cluster_method = "meshclust"`; `detect_chimeras = TRUE` |
| Python 3 + umap-learn, hdbscan, scikit-learn | UMAP pre-clustering | `cluster_method = "umap_meshclust"` |

**Resolution order**: host `PATH` → environment variable override → Apptainer container (see below).

Host install example (conda):

```bash
mamba create -n miaoseq-tools -c bioconda -c conda-forge \
    mmseqs2 minimap2 vsearch python=3.11 umap-learn hdbscan scikit-learn
conda activate miaoseq-tools
```

### Apptainer image (recommended for refless)

The repository ships an Apptainer definition that bundles mmseqs2, minimap2, vsearch, samtools, and Python UMAP dependencies.

**Prerequisites**: [Apptainer](https://apptainer.org/) (or Singularity) ≥ 1.1, and either `mksquashfs` (for `.sif`) or use the sandbox directly.

```bash
cd cursor_dev/apptainer
bash build.sh
```

This creates:

- `cursor_dev/apptainer/images/miaoseq-refless.sandbox` — always built
- `cursor_dev/apptainer/images/miaoseq-refless.sif` — created if SIF conversion succeeds

Point miaoseq to the image:

```bash
# Prefer SIF when available; fall back to sandbox
export MIAOSEQ_REFLESS_SIF=/path/to/miaoseq-refless.sif
# or
export MIAOSEQ_REFLESS_SIF=/path/to/miaoseq-refless.sandbox
```

If SIF conversion fails (`mksquashfs` error), the sandbox works identically:

```bash
export MIAOSEQ_REFLESS_SIF=cursor_dev/apptainer/images/miaoseq-refless.sandbox
```

### Environment variables

| Variable | Description |
|----------|-------------|
| `MIAOSEQ_REFLESS_SIF` | Path to Apptainer `.sif` or `.sandbox` image |
| `MIAOSEQ_MMSEQS_PATH` | Override mmseqs2 binary |
| `MIAOSEQ_MINIMAP2_PATH` | Override minimap2 binary |
| `MIAOSEQ_VSEARCH_PATH` | Override vsearch binary |
| `MIAOSEQ_UMAP_SCRIPT` | Override path to `umap_kmer_cluster.py` |
| `REFLESS_FASTQ` | Convenience variable used in sample scripts |
| `REFLESS_OUT` | Output directory for sample scripts |

---

## Input file formats

All CSV files below have **no header row**.

### Index list (`index_list`)

Five columns: `sample_id`, `forward_index_id`, `forward_index_seq`, `reverse_index_id`, `reverse_index_seq`.

```
miaoBC0001,miao_I7_index_001,CATACGAGATCGCTCAGTTCGTGACTGGAGTTCAGACGTGTG,miao_I5_index_001,CATACGAGATTCGTGGAGCGACACTCTTTCCCTACACGACG
```

Bundled files: `inst/extdata/index_list.csv` (full plate), `inst/extdata/index_list_smoke.csv` (2-sample smoke test).

### Primer list (`primer_list`)

Two columns: `primer_id` (must end with `_F` or `_R`), `sequence`.

Bundled: `inst/extdata/amplicon_primers.csv`.

### PAM / gRNA list (`pam_list`)

Used by edit-calling modules. Columns: `gene`, `chromosome`, `cut_position` [, optional `guide_id`].

Bundled: `inst/extdata/agr8_pam_list.csv`.

### Sample layout (for `editViewer`)

Five columns: `index_pair_id`, `sample_name`, `plate_id`, `row` (A–H), `col` (1–12).

```
miaoBC0001,Sample_A,1,A,1
miaoBC0002,Sample_B,1,A,2
```

---

## Usage: refless pipeline

The refless workflow processes **pre-basecalled multiplexed FASTQ** through demultiplexing, per-sample clustering, consensus calling, and optional QC reporting.

### Step 1 — Core pipeline

```r
library(miaoseq)

fastq_fn    <- "multiplexed.fastq.gz"
out_dir     <- "refless_output"
index_list  <- system.file("extdata", "index_list.csv", package = "miaoseq")
primer_list <- system.file("extdata", "amplicon_primers.csv", package = "miaoseq")
sif_path    <- Sys.getenv("MIAOSEQ_REFLESS_SIF", "cursor_dev/apptainer/images/miaoseq-refless.sandbox")

results <- miaoPipeline(
    fastq_fn         = fastq_fn,
    index_list       = index_list,
    primer_list      = primer_list,
    out_dir          = out_dir,
    cluster_method   = "mmseqs2",       # or "meshclust", "umap_meshclust"
    sif_path         = if(file.exists(sif_path)) sif_path else NULL,
    detect_chimeras  = TRUE,            # vsearch uchime_denovo per sample
    write_report     = TRUE,            # HTML + CSV quality report
    demultiplex_params = list(
        end_window = 200,               # bp window at each read end
        t_high     = 0.90,
        t_low      = 0.75,
        delta      = 0.05,
        n_core     = 4
    ),
    cluster_params = list(
        refine             = TRUE,      # minimap2 refinement
        mmseqs_min_seq_id  = 0.85
    ),
    resume = FALSE
)

table(results$demultiplex$class)
```

A runnable template is in `cursor_dev/sample_script_refless.R`.

### Step 2 — Edit calling (optional)

Consensus-based edit calling replaces per-read BLAST alignment:

```r
genome_fn <- "/reference/genome.fa"
pam_list  <- system.file("extdata", "agr8_pam_list.csv", package = "miaoseq")

edit_results <- doEditCalling(
    pipeline_results = results,
    genome_fn        = genome_fn,
    primer_list      = primer_list,
    pam_list         = pam_list,
    out_dir          = file.path(out_dir, "edit_calling"),
    sif_path         = sif_path,
    prefer_edited    = TRUE
)
```

### Step 3 — Index design check (optional)

```r
ix <- evaluateIndexSet(index_list)
ix$summary   # per-index-pair minimum edit distances
ix$pairwise  # full distance matrix
```

### Key `miaoPipeline()` parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `cluster_method` | `"mmseqs2"` | `"mmseqs2"`, `"meshclust"`, or `"umap_meshclust"` |
| `detect_chimeras` | `FALSE` | Run vsearch uchime_denovo; exclude chimeras from clustering |
| `write_report` | `TRUE` | Generate `report/quality_report.html` and CSV summaries |
| `resume` | `FALSE` | Skip steps when output files already exist |
| `samples` | all High+Low | Subset of sample IDs to process |

---

## Usage: legacy pipeline

The legacy workflow starts from **pod5** raw data and uses BLAST throughout.

### Step 1 — Prepare amplicon database

```r
library(miaoseq)

out_dir    <- "legacy_output"
genome_fn  <- "/reference/genome.fa"
primer_list <- system.file("extdata", "amplicon_primers.csv", package = "miaoseq")
blast_path  <- "/path/to/blast/bin"
n_core      <- 30

amplicon_fn <- prepAmpliconDB(
    blast_path  = blast_path,
    primer_list = primer_list,
    genome_fn   = genome_fn,
    out_dir     = out_dir,
    n_core      = n_core
)
```

### Step 2 — Run full pipeline

```r
in_dir      <- "/data/minion_run/pod5"
pam_list    <- system.file("extdata", "agr8_pam_list.csv", package = "miaoseq")
index_list  <- system.file("extdata", "index_list.csv", package = "miaoseq")

miaoEditcall(
    in_dir       = in_dir,
    out_dir      = out_dir,
    dorado_path  = "/path/to/dorado",
    samtools_path = "/path/to/samtools",
    blast_path   = blast_path,
    primer_list  = primer_list,
    pam_list     = pam_list,
    index_list   = index_list,
    genome_fn    = genome_fn,
    amplicon_fn  = amplicon_fn,
    size_sel     = c(300, 450),
    check_window = 10,
    n_core       = n_core,
    resume       = FALSE
)
```

### Step 3 — Summary statistics

```r
evalMiao(out_dir = out_dir, output_reads = FALSE)
```

---

## Edit calling and visualization

### `editViewer()` — plate heatmaps

Generates per-plate PDF heatmaps from `editcall/editcall_summary.csv`.

```r
sample_layout <- "sample_layout.csv"   # 5-column plate map (no header)

# Legacy wide-format summary
editViewer(out_dir, sample_layout)

# Refless long-format summary (auto-detected when cluster_id column present)
editViewer(out_dir, sample_layout, refless = TRUE)

# Show all clusters per well (primary marked with *)
editViewer(out_dir, sample_layout,
           refless = TRUE, cluster_level = TRUE)

# Single PDF for all plates
editViewer(out_dir, sample_layout, onefile = TRUE)
```

Genotype color key: `ref` (yellow), out-of-frame edits (blue shades), in-frame edits (green shades).

### Quality report (refless)

When `write_report = TRUE`, `miaoPipeline()` writes:

- `report/quality_report.html` — human-readable run summary
- `report/run_summary.csv` — total/assigned/unclassified read counts
- `report/sample_quality.csv` — per-sample cluster and chimera metrics
- `report/contaminant_clusters.csv` — low-fraction cluster flags

Regenerate manually:

```r
writeReflessReport(
    pipeline_results = results,
    pipeline_out     = out_dir,
    edit_results     = edit_results   # optional
)
```

---

## Output directories

### Refless (`miaoPipeline`)

```
out_dir/
├── demultiplex/
│   ├── demultiplex_assignments.csv
│   └── per_sample/<sample_id>.fastq
├── samples/<sample_id>/
│   ├── cluster/          # cluster_assignments.csv, chimera/ (optional)
│   ├── consensus/        # per-cluster FASTA + consensus_summary.csv
│   └── confidence/       # confidence_metrics.csv
├── report/               # quality_report.html, CSV summaries
└── miaoPipeline_results.rds
```

After `doEditCalling()`:

```
edit_calling/
├── ref_targets/
├── editcall/editcall_summary.csv
└── samples/<sample_id>/
    ├── align/
    ├── editcall/
    ├── primary_cluster.csv
    └── cluster_scores.csv
```

### Legacy (`miaoEditcall`)

```
out_dir/
├── basecall/
├── demultiplex/
├── align/
├── editcall/
├── ref/
└── miao_summary/
```

---

## Troubleshooting

### External tool not found (refless)

1. Confirm binaries are on `PATH`, or set `MIAOSEQ_MMSEQS_PATH` / `MIAOSEQ_MINIMAP2_PATH`.
2. Build and export the Apptainer image (see [Apptainer image](#apptainer-image-recommended-for-refless)).
3. Run the smoke test with `index_list_smoke.csv` and a small synthetic FASTQ before full plates.

### SIF build fails

If `build.sh` reports `mksquashfs` failure, use the sandbox path:

```bash
export MIAOSEQ_REFLESS_SIF=cursor_dev/apptainer/images/miaoseq-refless.sandbox
```

### Demultiplexing assigns few reads (refless)

- Increase `end_window` (default 200) if index+primer spans are longer.
- Lower `t_high` / raise `t_low` thresholds cautiously.
- Verify index sequences match the actual library design (`evaluateIndexSet()`).

### `editViewer()` errors

- Ensure `editcall/editcall_summary.csv` exists (from `miaoEditcall` or `doEditCalling`).
- Check that `sample_layout` index IDs match `index_pair_id` in the summary.
- Install `ggplot2` and `tidyr`.

### BLAST / Dorado errors (legacy)

- Verify `blastn` and `makeblastdb` are in `blast_path`.
- Confirm pod5 directory path and Dorado model availability.

---

## Performance notes

**Legacy pipeline** (30 cores, single Flongle cell): approximately **1 hour**, **15–20 GB RAM**.

**Refless pipeline**: demultiplexing is parallelized (`n_core` in `demultiplex_params`). Per-sample clustering scales with sample count; `mmseqs2` is fastest for prototyping, `meshclust` for production parity. Chimera detection adds one vsearch pass per sample.

---

## Citation and contact

Author: Tomoyuki Furuta (f.tomoyuki@okayama-u.ac.jp)

For bug reports and feature requests, open an issue on [GitHub](https://github.com/tomoyukif/miaoseq).
