################################################################################
# Prepare indexed amplicon database
################################################################################

#' Prepare amplicon database from primer locations or a user FASTA
#'
#' Either (1) locates forward/reverse primers in a reference genome with exact
#' matching (`Biostrings::matchPattern` on both strands) and extracts amplicon
#' intervals, or (2) accepts a user-supplied amplicon FASTA whose sequence names
#' match gene IDs inferred from the primer list (`*_F` / `*_R` stripped).
#'
#' When deriving intervals from the genome, the run **stops** if any gene has
#' more than one F/R locus pair, or if any interval length exceeds
#' `2 * expected_length`.
#'
#' @param primer_list Path to a CSV with primer ID and sequence columns.
#'   Forward primer IDs must end with \code{_F}, reverse with \code{_R}.
#' @param genome_fn Path to the reference genome FASTA (required unless
#'   \code{amplicon_fasta} is supplied).
#' @param out_dir Output directory; results are written under \code{ref/}.
#' @param amplicon_fasta Optional user-supplied amplicon reference FASTA.
#'   Sequence names must equal the gene IDs from \code{primer_list}. When set,
#'   genome-based locus discovery is skipped and this FASTA is copied to
#'   \code{ref/amplicon.fa}.
#' @param expected_length Expected amplicon length (bp). Required when deriving
#'   amplicons from \code{genome_fn}; intervals longer than
#'   \code{2 * expected_length} cause a hard stop. Optional when
#'   \code{amplicon_fasta} is supplied (then used only to validate sequence
#'   lengths against the same \code{2×} rule).
#' @return Path to \code{ref/amplicon.fa}.
#' @export
#' @import Biostrings
#' @import dplyr
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom BiocGenerics width start end
#' @importFrom Biostrings DNAString DNAStringSet readDNAStringSet writeXStringSet
#' @importFrom Biostrings reverseComplement matchPattern
#' @importFrom utils read.csv write.csv
#'
#' @examples
#' # amplicon_fn <- prepAmpliconDB(
#' #   primer_list = "inst/extdata/amplicon_primers.csv",
#' #   genome_fn = "/path/to/genome.fa",
#' #   out_dir = "/path/to/run",
#' #   expected_length = 500L
#' # )
#' # # Or supply curated amplicons directly:
#' # amplicon_fn <- prepAmpliconDB(
#' #   primer_list = "inst/extdata/amplicon_primers.csv",
#' #   out_dir = "/path/to/run",
#' #   amplicon_fasta = "/path/to/amplicons.fa"
#' # )
prepAmpliconDB <- function(primer_list,
                           genome_fn = NULL,
                           out_dir,
                           amplicon_fasta = NULL,
                           expected_length = NULL) {
  gene_ids <- .primer_gene_ids(primer_list)

  ref_dir <- file.path(out_dir, "ref")
  dir.create(ref_dir, recursive = TRUE, showWarnings = FALSE)
  amplicon_fn <- file.path(ref_dir, "amplicon.fa")

  if (!is.null(amplicon_fasta)) {
    if (!file.exists(amplicon_fasta)) {
      stop("amplicon_fasta does not exist: ", amplicon_fasta)
    }
    amplicon_seq <- readDNAStringSet(amplicon_fasta)
    .validate_user_amplicon_fasta(amplicon_seq, gene_ids, expected_length)
    writeXStringSet(amplicon_seq, amplicon_fn)
    if (!is.null(genome_fn)) {
      if (!file.exists(genome_fn)) {
        stop("genome_fn does not exist: ", genome_fn)
      }
      genome <- readDNAStringSet(genome_fn)
      writeXStringSet(genome, file.path(ref_dir, "ref_genome.fa"))
    }
    return(amplicon_fn)
  }

  if (is.null(genome_fn) || !nzchar(genome_fn)) {
    stop("genome_fn is required when amplicon_fasta is not supplied.")
  }
  if (is.null(expected_length) || length(expected_length) != 1L ||
      is.na(expected_length) || as.numeric(expected_length) <= 0) {
    stop(
      "expected_length (positive scalar, bp) is required when deriving ",
      "amplicons from genome_fn. Intervals longer than 2 * expected_length ",
      "are rejected."
    )
  }
  expected_length <- as.integer(expected_length)
  max_span <- 2L * expected_length

  primers <- read.csv(primer_list, header = FALSE, stringsAsFactors = FALSE)
  f_idx <- grep("_F$", primers$V1)
  r_idx <- grep("_R$", primers$V1)
  if (length(f_idx) < 1L || length(r_idx) < 1L) {
    stop("primer_list must contain IDs ending with _F and _R.")
  }

  f_primers <- DNAStringSet(primers$V2[f_idx])
  names(f_primers) <- primers$V1[f_idx]
  r_primers <- DNAStringSet(primers$V2[r_idx])
  names(r_primers) <- primers$V1[r_idx]

  genome <- readDNAStringSet(genome_fn)
  if (length(genome) < 1L) {
    stop("genome_fn contains no sequences: ", genome_fn)
  }

  f_out <- .collect_primer_genome_hits(f_primers, genome)
  r_out <- .collect_primer_genome_hits(r_primers, genome)
  if (is.null(f_out) || is.null(r_out) || nrow(f_out) < 1L || nrow(r_out) < 1L) {
    stop("No exact primer matches found in genome_fn.")
  }

  f_out$qseqid <- sub("_F$", "", f_out$primer_id)
  r_out$qseqid <- sub("_R$", "", r_out$primer_id)

  amplicon_df <- full_join(
    f_out[, c("qseqid", "sseqid", "sstart", "send")],
    r_out[, c("qseqid", "sseqid", "sstart", "send")],
    by = c("qseqid", "sseqid"),
    suffix = c(".x", ".y")
  )
  amplicon_df <- amplicon_df[
    !is.na(amplicon_df$sstart.x) & !is.na(amplicon_df$sstart.y),
    ,
    drop = FALSE
  ]
  if (nrow(amplicon_df) < 1L) {
    stop("No F/R primer pairs could be joined on the same chromosome.")
  }

  coord_mat <- as.matrix(amplicon_df[, c("sstart.x", "send.x", "sstart.y", "send.y")])
  amplicon_df$start <- apply(coord_mat, 1L, min, na.rm = TRUE)
  amplicon_df$end <- apply(coord_mat, 1L, max, na.rm = TRUE)
  amplicon_df$span <- as.integer(amplicon_df$end - amplicon_df$start + 1L)

  n_per_gene <- table(amplicon_df$qseqid)
  multi <- names(n_per_gene)[n_per_gene > 1L]
  if (length(multi) > 0L) {
    stop(
      "Multiple F/R amplicon loci found for gene(s): ",
      paste(multi, collapse = ", "),
      ". Resolve primer specificity or supply amplicon_fasta.",
      call. = FALSE
    )
  }

  too_long <- amplicon_df$span > max_span
  if (any(too_long)) {
    bad <- amplicon_df[too_long, , drop = FALSE]
    stop(
      "Amplicon interval(s) exceed 2 * expected_length (", max_span, " bp): ",
      paste0(bad$qseqid, "=", bad$span, "bp", collapse = "; "),
      ". Supply a curated amplicon_fasta or adjust expected_length.",
      call. = FALSE
    )
  }

  missing_genes <- setdiff(gene_ids, amplicon_df$qseqid)
  if (length(missing_genes) > 0L) {
    stop(
      "No unique amplicon locus for gene(s): ",
      paste(missing_genes, collapse = ", "),
      call. = FALSE
    )
  }

  amplicon_gr <- GRanges(
    seqnames = amplicon_df$sseqid,
    ranges = IRanges(start = amplicon_df$start, end = amplicon_df$end),
    id = amplicon_df$qseqid
  )

  writeXStringSet(genome, file.path(ref_dir, "ref_genome.fa"))

  amplicon_seq <- genome[amplicon_gr]
  names(amplicon_seq) <- amplicon_gr$id

  write.csv(
    amplicon_df[, c("qseqid", "sseqid", "sstart.x", "send.x", "sstart.y", "send.y", "start", "end", "span")],
    file.path(ref_dir, "amplicon.csv"),
    row.names = FALSE
  )
  writeXStringSet(amplicon_seq, amplicon_fn)
  amplicon_fn
}

#' @keywords internal
.primer_gene_ids <- function(primer_list) {
  primers <- read.csv(primer_list, header = FALSE, stringsAsFactors = FALSE)
  f_ids <- primers$V1[grep("_F$", primers$V1)]
  r_ids <- primers$V1[grep("_R$", primers$V1)]
  if (length(f_ids) < 1L || length(r_ids) < 1L) {
    stop("primer_list must contain IDs ending with _F and _R.")
  }
  f_genes <- sub("_F$", "", f_ids)
  r_genes <- sub("_R$", "", r_ids)
  miss_r <- setdiff(f_genes, r_genes)
  miss_f <- setdiff(r_genes, f_genes)
  if (length(miss_r) || length(miss_f)) {
    stop(
      "primer_list has unpaired genes: missing _R for ",
      paste(miss_r, collapse = ", "),
      "; missing _F for ",
      paste(miss_f, collapse = ", ")
    )
  }
  sort(unique(f_genes))
}

#' @keywords internal
.validate_user_amplicon_fasta <- function(amplicon_seq, gene_ids, expected_length) {
  ids <- names(amplicon_seq)
  if (is.null(ids) || any(!nzchar(ids))) {
    stop("amplicon_fasta sequences must have non-empty names matching gene IDs.")
  }
  if (any(duplicated(ids))) {
    stop("amplicon_fasta has duplicated sequence names: ",
         paste(unique(ids[duplicated(ids)]), collapse = ", "))
  }
  missing <- setdiff(gene_ids, ids)
  extra <- setdiff(ids, gene_ids)
  if (length(missing) > 0L) {
    stop(
      "amplicon_fasta is missing gene ID(s) from primer_list: ",
      paste(missing, collapse = ", ")
    )
  }
  if (length(extra) > 0L) {
    stop(
      "amplicon_fasta has sequence name(s) not in primer_list gene IDs: ",
      paste(extra, collapse = ", ")
    )
  }
  if (!is.null(expected_length) && length(expected_length) == 1L &&
      !is.na(expected_length) && as.numeric(expected_length) > 0) {
    max_span <- 2L * as.integer(expected_length)
    lens <- as.integer(width(amplicon_seq))
    too_long <- lens > max_span
    if (any(too_long)) {
      stop(
        "amplicon_fasta sequence(s) exceed 2 * expected_length (", max_span,
        " bp): ",
        paste0(ids[too_long], "=", lens[too_long], "bp", collapse = "; ")
      )
    }
  }
  invisible(NULL)
}

.collect_primer_genome_hits <- function(primers, genome) {
  parts <- lapply(seq_along(primers), function(i) {
    primer_id <- names(primers)[i]
    primer <- primers[[i]]
    chr_hits <- lapply(seq_along(genome), function(j) {
      chr_name <- names(genome)[j]
      if (!nzchar(chr_name)) {
        chr_name <- paste0("chr", j)
      }
      .find_exact_primer_hits(primer, genome[[j]], chr_name)
    })
    chr_hits <- chr_hits[!vapply(chr_hits, is.null, logical(1))]
    if (length(chr_hits) < 1L) {
      return(NULL)
    }
    out <- do.call(rbind, chr_hits)
    out$primer_id <- primer_id
    out
  })
  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (length(parts) < 1L) {
    return(NULL)
  }
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out
}

.find_exact_primer_hits <- function(primer, chr_seq, chr_name) {
  pat <- DNAString(primer)
  subj <- DNAString(chr_seq)
  if (length(pat) < 1L || length(subj) < 1L) {
    return(NULL)
  }

  rows <- list()
  m_fwd <- matchPattern(pat, subj, fixed = TRUE)
  if (length(m_fwd) > 0L) {
    rows[[length(rows) + 1L]] <- data.frame(
      sseqid = chr_name,
      sstart = start(m_fwd),
      send = end(m_fwd),
      stringsAsFactors = FALSE
    )
  }

  subj_rc <- reverseComplement(subj)
  m_rev <- matchPattern(pat, subj_rc, fixed = TRUE)
  if (length(m_rev) > 0L) {
    chr_len <- length(subj)
    rows[[length(rows) + 1L]] <- data.frame(
      sseqid = chr_name,
      sstart = chr_len - end(m_rev) + 1L,
      send = chr_len - start(m_rev) + 1L,
      stringsAsFactors = FALSE
    )
  }

  if (length(rows) < 1L) {
    return(NULL)
  }
  do.call(rbind, rows)
}
