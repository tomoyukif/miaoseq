# miaoseq 0.2.0

## Major changes

* Redesign the end-to-end workflow around user-provided FASTQ:
  `doDemultiplex` → `doAssignGenes` → `doEditcall` (primary) and optional
  `doAssembleAmplicons` / `doReassessAssemblies`.
* Demultiplex with C++ edlib barcode assignment (`doDemultiplex`); drop the
  BLAST-based demux path.
* **Demultiplex I/O:** stream FASTQ from disk in chunks (no Biostrings full
  load); write `assignments.tsv` / `unassigned.tsv` incrementally; when
  `split_reads = TRUE`, write per-sample FASTQ in the same pass (C++ buffered
  gzip/plain). `splitDemultiplexReads()` also uses C++ streaming.
  New `return_tables` (default `TRUE`) can be set `FALSE` to avoid reloading
  large tables into R.
* Assemble with **vsearch** clustering and **abpoa** consensus only. Remove
  internal / mmseqs assemble backends, spoa, abPOA quality weighting (`-Q`),
  and racon polish.
* Split gene assignment (`doAssignGenes`) from assembly; add
  `doReassessAssemblies` for independent clustering QC.

## API / I/O changes

* `doEditcall()` retains gene primers in the aligned amplicon span (outer
  adapters only are trimmed via `include_primers = TRUE`). PAM / cut windows
  are coordinates on the full `amplicon.fa` sequence, so PAM sites near a
  primer end (e.g. RAE3) keep a full `check_window` rather than being clipped
  to a short insert-only fragment.
* `pam_list`: column 4 is strand, column 5 is guide ID. A gene with one PAM
  row is labelled by gene only. Bundled `agr8_pam_list.csv` includes Cas9
  strand. Editcall tables emit `sample_id` and `index_pair_id` (filled from
  demultiplex `assignments.tsv`).
* `doAssembleAmplicons()`: `min_cluster_identity` (default `0.95`) replaces
  `max_cluster_edit`; `max_clusters` defaults to `Inf`;
  `fraction_sample` renamed to `fraction_bucket`.
* `min_cluster_identity` accepts **`(0, 1]`** and is passed through to vsearch
  `--id` without a 0.5 floor (low-identity merging sweeps are supported).
* `clusters.fasta` now stores member reads with headers
  `>{read_id} {cluster_id}`.
* Expand Assemble `run_stats.txt` with R / package / tool versions and key
  clustering parameters.
* `doAssignGenes()` with `amplicon_fn`: primer-assigned reads whose trimmed
  insert length is outside `expected ± length_tolerance` (default 0.25) are
  labeled `length_outlier` and **kept**. Expected length is amplicon width
  minus F/R primer lengths. `doEditcall` / `doAssembleAmplicons` include
  these reads by default. Length labeling batches primer trim per gene
  (`n_core`). Primer search skips reverse-orientation / worse genes after a
  perfect hit. `store_sequences = FALSE` omits oriented sequences from the
  TSV (Editcall-only; Assemble from disk still needs `TRUE`).
* `doEditcall()` writes `run_stats.txt`. With three or more guides, a middle
  guide involved in two adjacent `both_cut_excision` events uses Plan A
  allele `---`.
* `doEditcall()` per-read path is C++ / OpenMP (`edlib` global NW, adaptive
  window, Plan A/A′). `n_core` applies to primer trim and edit extraction.
* `doReassessAssemblies()`: when `min_identity` is `NULL`, derive identity as
  `1 - edit / median_consensus_length` via `.vsearch_id_from_edit()` and clamp
  to `[0.50, 0.99]` (Assemble pass-through identity is unchanged).
* Remove legacy entry points (`doBasecall`, BLAST demux helpers, monolithic
  pipeline wrappers) in favour of the modular API above.

## Documentation

* Refresh README / man / `cursor_dev/pipeline_revise.md` (§20–§21) for
  redesigned paths, `amplicon_assign/` layout, Assemble identity tuning
  (including sub-0.5 values), gene-assign `length_outlier` labeling,
  `store_sequences`, Editcall `run_stats.txt` / C++ path, and Reassess
  identity derivation.

# miaoseq 0.1.0

* Initial modular pipeline toward demultiplex / amplicon resolve / editcall
  (pre–Phase K Assemble formalization).
