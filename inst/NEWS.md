# miaoseq 0.4.0

## Breaking / API

* `primer_list` is required for `doAssignGenes()` and `doAssembleAmplicons()`.
  `NULL` / `no_ref` / unknown-gene assignment is an error.
* Editcall column `cut_insert` renamed to `cut_amplicon`.
* Assemble default is `strict_end_trim = TRUE` (cluster trimmed inserts only;
  trim failures go to U1). Reassess Q2 prep uses the same rule.

## Behavior

* `prepAmpliconDB()` orients genome-derived (and user) amplicons so the
  forward primer is at the 5′ end. `amplicon.fa` is primer-inclusive F→R.
* Gene-assign splits `trim_fail` vs `length_outlier`. Downstream processable
  statuses are `assigned`, `length_outlier`, and `trim_fail`.
* Editcall summary / Viewer keep `---` / `excision` (not `indel3-3`).
* `max_barcode_edit > 2` is implemented, capped at barcode length.
* Demux / Assign write `run_stats.txt`. Reassess records R / miaoseq versions.
* List CSVs (`index_list`, `primer_list`, `pam_list`, `sample_list`) accept an
  optional header row.
* Legacy `sample_id` → `index_pair_id` copy stops when labels are duplicated.
* Incomplete trailing FASTQ records are counted (warning + demux `run_stats`).
* Post-hoc split maps on `read_id` + `source_file`. `by_sample` write/read
  share one filename sanitizer.
* evalMiao QC uses processable statuses; editcall yield is unique reads when
  `read_id` is available.

## Internals

* Assign / Assemble / Reassess internals that held well IDs are named
  `index_pair_id` (public `sample_id` labels unchanged).
* Removed unused C++ majority consensus / amplicon-ref assignment exports.
* Gapped edit windows use an ungapped walk; length safety-net keeps large
  insertions and excision tokens.

# miaoseq 0.3.1

## API / I/O changes

* Pipeline grouping key is **`index_pair_id`**. `sample_list$sample_name` /
  `sample_id` may be duplicated across wells. `doDemultiplex()` `by_sample/`
  FASTQ and `summary_by_sample.tsv`, `splitDemultiplexReads()`, assign /
  assemble / reassess, and editcall I/O all key on `index_pair_id`.
  `sample_id` remains a per-read label in `assignments.tsv` and is attached
  only as a reference label on editcall tables.
* `doAssignGenes()` writes `index_pair_id` (not `sample_id`) in
  `gene_assignments.tsv`. Legacy tables with only `sample_id` are still
  accepted and copied to `index_pair_id`.
* Assemble / Reassess per-well directories and summary columns use
  `index_pair_id` (`amplicon/{index_pair_id}/`, consensus header
  `index_pair=`).

## Documentation

* README, man pages, and `cursor_dev` plans (`pipeline_revise.md` §22,
  `demux_revise.md`, `ampliconresolve_plan.md`, `editcall_adapter.md`)
  synced for the `index_pair_id` grouping contract.

# miaoseq 0.3.0

## Major changes

* **Demultiplex I/O:** stream FASTQ from disk in chunks (no Biostrings full
  load); write `assignments.tsv` / `unassigned.tsv` incrementally; when
  `split_reads = TRUE`, write per-well FASTQ in the same pass (C++ buffered
  gzip/plain). `splitDemultiplexReads()` also uses C++ streaming.
  New `return_tables` (default `TRUE`) can be set `FALSE` to avoid reloading
  large tables into R.
* **Editcall speed:** per-read path is C++ / OpenMP (`edlib` global NW,
  adaptive window, Plan A/A′ in `src/editcall_core.cpp`). `n_core` applies to
  primer-span trim and edit extraction.
* **AssignGenes speed:** primer search skips reverse-orientation / worse genes
  after a perfect hit; length labeling batches primer trim per gene (`n_core`).

## API / I/O changes

* `doEditcall()` retains gene primers in the aligned amplicon span (outer
  adapters only are trimmed via `include_primers = TRUE`). PAM / cut windows
  are coordinates on the full `amplicon.fa` sequence, so PAM sites near a
  primer end keep a full `check_window`.
* `pam_list`: column 4 is strand, column 5 is guide ID. A gene with one PAM
  row is labelled by gene only (`target_gene = gene`). Bundled
  `agr8_pam_list.csv` includes Cas9 strand. Editcall tables emit `sample_id`
  and `index_pair_id` (filled from demultiplex `assignments.tsv`).
* `doEditcall()` writes `run_stats.txt`. With three or more guides, a middle
  guide involved in two adjacent `both_cut_excision` events uses Plan A
  allele `---`.
* `doAssignGenes()` with `amplicon_fn`: primer-assigned reads whose trimmed
  insert length is outside `expected ± length_tolerance` (default 0.25) are
  labeled `length_outlier` and **kept**. Expected length is amplicon width
  minus F/R primer lengths. `doEditcall` / `doAssembleAmplicons` include
  these reads by default.
* `doAssignGenes(store_sequences = FALSE)` omits oriented sequences from the
  TSV (Editcall-only; Assemble / Reassess from disk still need `TRUE`).
* `min_cluster_identity` accepts **`(0, 1]`** and is passed through to vsearch
  `--id` without a 0.5 floor (low-identity merging sweeps are supported).
* `doReassessAssemblies()`: when `min_identity` is `NULL`, derive identity as
  `1 - edit / median_consensus_length` via `.vsearch_id_from_edit()` and clamp
  to `[0.50, 0.99]` (Assemble pass-through identity is unchanged).

## Documentation

* Refresh README / man for `length_outlier`, `store_sequences`, Editcall
  `run_stats.txt` / C++ path, PAM column order, and Reassess identity
  derivation.

# miaoseq 0.2.0

## Major changes

* Redesign the end-to-end workflow around user-provided FASTQ:
  `doDemultiplex` → `doAssignGenes` → `doEditcall` (primary) and optional
  `doAssembleAmplicons` / `doReassessAssemblies`.
* Demultiplex with C++ edlib barcode assignment (`doDemultiplex`); drop the
  BLAST-based demux path.
* Assemble with **vsearch** clustering and **abpoa** consensus only. Remove
  internal / mmseqs assemble backends, spoa, abPOA quality weighting (`-Q`),
  and racon polish.
* Split gene assignment (`doAssignGenes`) from assembly; add
  `doReassessAssemblies` for independent clustering QC.

## API / I/O changes

* `doAssembleAmplicons()`: `min_cluster_identity` (default `0.95`) replaces
  `max_cluster_edit`; `max_clusters` defaults to `Inf`;
  `fraction_sample` renamed to `fraction_bucket`.
* `clusters.fasta` now stores member reads with headers
  `>{read_id} {cluster_id}`.
* Expand Assemble `run_stats.txt` with R / package / tool versions and key
  clustering parameters.
* Remove legacy entry points (`doBasecall`, BLAST demux helpers, monolithic
  pipeline wrappers) in favour of the modular API above.

## Documentation

* Refresh README for the redesigned paths, Assemble defaults, and Reassess
  usage.

# miaoseq 0.1.0

* Initial modular pipeline toward demultiplex / amplicon resolve / editcall
  (pre–Phase K Assemble formalization).
