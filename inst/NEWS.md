# miaoseq News

## miaoseq 0.1.0 (2026-07-03)

### Refless cluster-consensus pipeline (new)

- Added `miaoPipeline()` to orchestrate demultiplexing, per-sample clustering,
  consensus generation, and confidence evaluation.
- Added `doDemultiplex2()` with index + primer end-window alignment,
  edit-distance scoring, and four-tier read classification (High / Low /
  Ambiguous / Unclassified).
- Added `doCluster()` with three clustering backends: `mmseqs2`,
  `meshclust` (vsearch), and `umap_meshclust`, plus optional minimap2
  cluster refinement.
- Added `doConsensus()` for pileup-based per-cluster consensus sequences.
- Added `evalConfidence()` for per-cluster support and identity metrics.
- Added edit-calling module: `prepRefTargets()`, `alignConsensusToRef()`,
  `callEditsFromConsensus()`, `selectEditCluster()`, and `doEditCalling()`.
- Added `evaluateIndexSet()` for index design QC (minimum edit distance).
- Added Apptainer definition (`cursor_dev/apptainer/`) bundling mmseqs2,
  minimap2, vsearch, samtools, and UMAP/HDBSCAN dependencies.
- Added smoke-test index list (`inst/extdata/index_list_smoke.csv`) and UMAP
  pre-cluster script (`inst/scripts/umap_kmer_cluster.py`).

### Quality control and chimera detection

- Added `detectChimeras()` using vsearch `uchime_denovo`; integrates with
  `doCluster()` via `detect_chimeras = TRUE` to exclude chimera reads from
  clustering.
- Added `writeReflessReport()` to generate per-run CSV summaries and an HTML
  quality report (demultiplex rates, cluster counts, contaminant flags,
  chimera statistics).

### Visualization

- Extended `editViewer()` for refless long-format `editcall_summary.csv`:
  - `refless = TRUE` — auto-detect or force refless summary format.
  - `cluster_level = TRUE` — show multiple clusters per well with fraction
    labels; primary clusters marked with `*`.
  - `primary_only = TRUE` — show only primary clusters (default).
- `doEditCalling()` now writes `editcall/editcall_summary.csv` compatible
  with `editViewer()`.

### Other

- Added `makeMMI()` utility for minimap2 index generation.
- Removed unused `doAlign` argument.

---

## miaoseq 0.0.5 (2025-12-29)

- Improved `editViewer()` plate layout and labeling.
- Added `makeMMI()` for building minimap2 `.mmi` indexes from a reference
  FASTA.

---

## miaoseq 0.0.4 (2025-10-17)

- Added total read count to edit-calling summary and `editViewer()` plate
  annotations.
- Fixed total read count calculation in summary tables.
- Fixed minor bug in `doEditcall()`.
- Updated `editViewer()` to omit unnamed samples when checking genotype label
  uniqueness.

---

## miaoseq 0.0.3 (2025-10-13)

- Added `editViewer()` for per-plate PDF heatmaps of edit-calling results.
- General code cleanup and minor bug fixes.

---

## miaoseq 0.0.2 (2025-10-07)

- Minor bug fixes in demultiplexing and edit-calling steps.

---

## miaoseq 0.0.1 (2025-09-30)

- Initial public release of the miaoseq R package.
- Legacy read-level pipeline: `doBasecall()`, `doDemultiplex()`, `doAlign()`,
  `doEditcall()`, orchestrated by `miaoEditcall()`.
- Reference preparation with `prepAmpliconDB()`.
- Run evaluation with `evalMiao()`.
- Bundled example data in `inst/extdata/` (index list, primer list, PAM list).
