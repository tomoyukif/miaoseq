% Generated manually
\name{doAssembleAmplicons}
\alias{doAssembleAmplicons}
\alias{doAmpliconResolve}
\title{Assemble amplicons from gene-assigned reads (Step 2b)}
\usage{
doAssembleAmplicons(
  gene_assign,
  out_dir,
  primer_list = NULL,
  method = c("both", "cluster", "consensus"),
  min_reads = 5L,
  min_cluster_reads = 5L,
  max_clusters = Inf,
  min_cluster_identity = 0.95,
  max_primer_edit = 5L,
  cluster_backend = "vsearch",
  assembly_backend = c("cluster", "overlap_graph"),
  overlap_min_identity = 0.90,
  strict_end_trim = FALSE,
  min_cluster_purity = NULL,
  overwrite = FALSE,
  n_core = 1L
)

doAmpliconResolve(gene_assign, out_dir, ...)
}
\description{
Consumes GeneAssignResult (from \code{doAssignGenes}) or \code{gene_assignments.tsv}
and performs per sample x gene clustering (\code{vsearch}) and consensus
reconstruction (\code{abpoa}). Clusters whose exact-insert modal purity is below
\code{min_cluster_purity} are discarded (\code{NULL} disables the filter).
\code{doAmpliconResolve} is an alias of this assembly-only function.
}
\arguments{
\item{min_cluster_reads}{Minimum reads per cluster (default 5). Minor clusters
below this size are excluded.}

\item{max_clusters}{Maximum clusters to keep per sample x gene. Default
\code{Inf} (no top-N truncation).}

\item{min_cluster_identity}{Minimum vsearch identity (\code{--id}, default 0.95).}

\item{cluster_backend}{Fixed to \code{"vsearch"} (Phase K).}

\item{min_cluster_purity}{Minimum fraction of cluster members that share the
most common exact insert sequence. Default \code{NULL} (disabled). Set to
e.g. \code{0.8} for strict filtering.}
}
\details{
\code{clusters.fasta} contains member reads with headers
\code{>{read_id} {cluster_id}} (space-separated). \code{consensus.fasta}
contains abPOA consensus sequences. In \code{cluster_counts.tsv},
\code{fraction_bucket} is the fraction within the sample x gene bucket
(formerly \code{fraction_sample}).
}
\value{A list with \code{samples}, \code{out_dir}, assembled \code{table}, and \code{gene_assignments}.
\code{cluster_counts.tsv} includes a \code{cluster_purity} column when purity was computed.}
