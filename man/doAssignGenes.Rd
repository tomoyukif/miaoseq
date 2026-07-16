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
  primer_list = NULL,
  amplicon_fn = NULL,
  max_primer_edit = 5L,
  end_window = 150L,
  length_tolerance = 0.25,
  samples = NULL,
  n_core = 1L,
  overwrite = FALSE,
  stats_unassign = FALSE
)
}
\description{
Reads demultiplexed FASTQ records, assigns each read to a target gene by primer matching, applies orientation normalization, and optionally filters by trimmed insert length against \code{amplicon_fn} widths. Amplicon-ref fallback assignment has been removed.
}
\value{A GeneAssignResult list with \code{samples}, \code{out_dir}, \code{table}, and optional \code{stats_unassigned}.}
