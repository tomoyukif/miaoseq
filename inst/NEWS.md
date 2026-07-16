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
