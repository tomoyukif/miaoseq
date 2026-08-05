% Generated manually
\name{doAssignGenes}
\alias{doAssignGenes}
\title{Assign reads to genes (Step 2a)}
\usage{
doAssignGenes(
  assignments,
  out_dir,
  fastq = NULL,
  sample_fastq_dir = NULL,
  primer_list,
  amplicon_fn = NULL,
  max_primer_edit = 5L,
  end_window = 150L,
  length_tolerance = 0.25,
  samples = NULL,
  n_core = 1L,
  store_sequences = TRUE,
  overwrite = FALSE,
  stats_unassign = FALSE
)
}
\description{
Reads demultiplexed FASTQ records, assigns each read to a target gene by primer matching, applies orientation normalization, and optionally labels trim failures (\code{trim_fail}) and length outliers (\code{length_outlier}) against \code{amplicon_fn} (amplicon width minus primer lengths). Neither status is dropped. \code{primer_list} is required (\code{NULL} / \code{no_ref} is not supported). Writes \code{run_stats.txt}. Processing and \code{gene_assignments.tsv} are keyed by \code{index_pair_id} (not \code{sample_id}). Amplicon-ref fallback assignment has been removed. \code{n_core} parallelizes primer assignment and length-trim. \code{store_sequences = FALSE} writes empty oriented sequence columns (Editcall can reload FASTQ; Assemble from disk needs sequences).
}
\value{A GeneAssignResult list with \code{samples}, \code{out_dir}, \code{table}, and optional \code{stats_unassigned}.}
