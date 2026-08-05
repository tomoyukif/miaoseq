% Generated manually
\name{prepAmpliconDB}
\alias{prepAmpliconDB}
\title{Prepare amplicon database from primer locations or a user FASTA}
\usage{
prepAmpliconDB(
  primer_list,
  genome_fn = NULL,
  out_dir,
  amplicon_fasta = NULL,
  expected_length = NULL
)
}
\arguments{
\item{primer_list}{Path to a CSV with primer ID and sequence columns.
Forward primer IDs must end with \code{_F}, reverse with \code{_R}.}

\item{genome_fn}{Path to the reference genome FASTA (required unless
\code{amplicon_fasta} is supplied).}

\item{out_dir}{Output directory; results are written under \code{ref/}.}

\item{amplicon_fasta}{Optional user-supplied amplicon reference FASTA.
Sequence names must equal the gene IDs from \code{primer_list}. When set,
genome-based locus discovery is skipped.}

\item{expected_length}{Expected amplicon length (bp). Required when deriving
amplicons from \code{genome_fn}; intervals longer than
\code{2 * expected_length} cause a hard stop. Optional when
\code{amplicon_fasta} is supplied.}
}
\value{Path to \code{ref/amplicon.fa}.}
\description{
Either locates forward/reverse primers in a reference genome with exact
matching and extracts amplicon intervals, or accepts a user-supplied amplicon
FASTA. Genome derivation stops on multi-locus hits or over-long spans.
Sequences are oriented so the forward primer is at the 5' end
(primer-inclusive F to R). User FASTA must match that contract.
}
\examples{
# amplicon_fn <- prepAmpliconDB(
#   primer_list = "inst/extdata/amplicon_primers.csv",
#   genome_fn = "/path/to/genome.fa",
#   out_dir = "/path/to/run",
#   expected_length = 500L
# )
}
