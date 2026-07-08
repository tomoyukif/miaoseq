################################################################################
# Step 2: Amplicon resolve (Phase C0 prototype)
################################################################################

#' Resolve amplicons per sample (Phase C0)
#'
#' Prototype implementation of Step 2 in the revised pipeline.
#' Takes demultiplexed assignments and FASTQ, groups reads by sample,
#' and performs trivial clustering (exact-sequence groups) to produce
#' per-sample consensus/cluster outputs.
#'
#' This Phase C0 version:
#' - Treats all reads as belonging to a single gene `unknown`
#'   (no primer/amplicon-based gene assignment yet)
#' - Uses exact-sequence grouping as clusters
#' - Uses the cluster representative sequence itself as the consensus
#'
#' @param assignments Path to `assignments.tsv` or a data.frame with at least
#'   `read_id` and `sample_id` columns.
#' @param out_dir Output directory, typically `{run_dir}/amplicon`.
#' @param fastq Character vector of source FASTQ paths. Required when
#'   `sample_fastq_dir` is NULL.
#' @param sample_fastq_dir Optional directory containing `by_sample/{sample_id}.fq`
#'   or `.fq.gz`. When provided and files exist, those are used preferentially.
#' @param method Character; one of `"both"`, `"cluster"`, `"consensus"`.
#'   Phase C0 treats `"both"` and `"cluster"` equivalently.
#' @param min_reads Integer; minimum reads per `sample_id` to process.
#' @param min_cluster_reads Integer; minimum reads per exact-sequence cluster
#'   to retain in outputs.
#' @param samples Optional character vector of `sample_id` to restrict processing.
#'   By default, all samples present in `assignments` are processed.
#' @param n_core Reserved for future parallelisation (currently unused).
#' @param overwrite Logical; if `FALSE` (default), stop when `out_dir` already
#'   exists and is non-empty.
#'
#' @return A list with fields:
#'   \item{samples}{Processed sample IDs}
#'   \item{out_dir}{Output directory path}
#'   \item{table}{Combined `cluster_counts.tsv` as a data.frame}
#'
#' @export
#' @importFrom utils write.table
#' @importFrom Biostrings DNAStringSet writeXStringSet
doAmpliconResolve <- function(
  assignments,
  out_dir,
  fastq = NULL,
  sample_fastq_dir = NULL,
  method = c("both", "cluster", "consensus"),
  min_reads = 5L,
  min_cluster_reads = 5L,
  samples = NULL,
  n_core = 1L,
  overwrite = FALSE
) {
  method <- match.arg(method)

  if (!overwrite && dir.exists(out_dir) && length(list.files(out_dir)) > 0L) {
    stop("out_dir already exists and is not empty. Set overwrite = TRUE to overwrite.")
  }
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  assn <- .as_assignment_df(assignments)
  if (!all(c("read_id", "sample_id") %in% names(assn))) {
    stop("assignments must contain read_id and sample_id columns.")
  }

  if (!is.null(samples)) {
    assn <- assn[assn$sample_id %in% samples, , drop = FALSE]
  }
  if (nrow(assn) < 1L) {
    warning("No assignments to process in doAmpliconResolve().")
    summary_path <- file.path(out_dir, "summary_by_sample.tsv")
    write.table(
      data.frame(),
      file = summary_path,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
    return(list(samples = character(), out_dir = out_dir, table = data.frame()))
  }

  sample_ids <- sort(unique(assn$sample_id))

  # Load sequences once: by_sample files preferred; otherwise one FASTQ pass
  sample_seqs <- .resolve_all_sample_sequences(
    sample_ids = sample_ids,
    assignments = assn,
    fastq = fastq,
    sample_fastq_dir = sample_fastq_dir
  )

  # Phase C0: single gene bucket per sample
  gene_id <- "unknown"

  combined_counts <- list()
  summary_rows <- vector("list", length(sample_ids))

  for (i in seq_along(sample_ids)) {
    sid <- sample_ids[[i]]
    sid_dir <- file.path(out_dir, sid)
    if (!dir.exists(sid_dir)) {
      dir.create(sid_dir, recursive = TRUE, showWarnings = FALSE)
    }

    seqs <- sample_seqs[[sid]]
    if (is.null(seqs)) seqs <- character()

    n_in <- length(seqs)
    if (n_in < min_reads) {
      stats_df <- .make_amplicon_stats(
        sample_id = sid,
        n_reads_in = n_in,
        n_reads_assigned_gene = 0L,
        n_genes_detected = 0L,
        n_clusters_total = 0L,
        skip_reason = "low_sample_reads"
      )
      utils::write.table(
        stats_df,
        file = file.path(sid_dir, "stats.tsv"),
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
      )
      summary_rows[[i]] <- data.frame(
        sample_id = sid,
        n_reads = n_in,
        n_clusters = 0L,
        skip_reason = "low_sample_reads",
        stringsAsFactors = FALSE
      )
      next
    }

    # Phase C0 clustering: exact sequence groups
    tbl <- table(seqs)
    tbl <- sort(tbl, decreasing = TRUE)
    keep <- tbl >= as.integer(min_cluster_reads)
    tbl_keep <- tbl[keep]

    if (length(tbl_keep) < 1L) {
      stats_df <- .make_amplicon_stats(
        sample_id = sid,
        n_reads_in = n_in,
        n_reads_assigned_gene = 0L,
        n_genes_detected = 0L,
        n_clusters_total = 0L,
        skip_reason = "no_clusters"
      )
      utils::write.table(
        stats_df,
        file = file.path(sid_dir, "stats.tsv"),
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
      )
      summary_rows[[i]] <- data.frame(
        sample_id = sid,
        n_reads = n_in,
        n_clusters = 0L,
        skip_reason = "no_clusters",
        stringsAsFactors = FALSE
      )
      next
    }

    # Build cluster table
    n_gene_reads <- sum(tbl_keep)
    clusters <- data.frame(
      sample_id = sid,
      gene_id = gene_id,
      cluster_id = seq_along(tbl_keep),
      seq = names(tbl_keep),
      n_reads = as.integer(unname(tbl_keep)),
      stringsAsFactors = FALSE
    )
    clusters$n_reads_gene <- n_gene_reads
    clusters$fraction <- clusters$n_reads / n_gene_reads
    clusters$fraction_sample <- clusters$n_reads / n_in
    clusters$method <- "exact"

    # consensus.fasta / clusters.fasta
    dna <- Biostrings::DNAStringSet(clusters$seq)
    names(dna) <- sprintf(
      "%s|cluster_%d;size=%d;sample=%s",
      clusters$gene_id,
      clusters$cluster_id,
      clusters$n_reads,
      clusters$sample_id
    )

    # In Phase C0, consensus == cluster representative
    Biostrings::writeXStringSet(
      dna,
      filepath = file.path(sid_dir, "consensus.fasta")
    )
    if (method %in% c("both", "cluster")) {
      Biostrings::writeXStringSet(
        dna,
        filepath = file.path(sid_dir, "clusters.fasta")
      )
    }

    counts_path <- file.path(sid_dir, "cluster_counts.tsv")
    utils::write.table(
      clusters[, c(
        "sample_id", "gene_id", "cluster_id",
        "n_reads", "fraction", "fraction_sample",
        "n_reads_gene", "method"
      )],
      file = counts_path,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )

    stats_df <- .make_amplicon_stats(
      sample_id = sid,
      n_reads_in = n_in,
      n_reads_assigned_gene = n_gene_reads,
      n_genes_detected = 1L,
      n_clusters_total = nrow(clusters),
      skip_reason = ""
    )
    utils::write.table(
      stats_df,
      file = file.path(sid_dir, "stats.tsv"),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )

    summary_rows[[i]] <- data.frame(
      sample_id = sid,
      n_reads = n_in,
      n_clusters = nrow(clusters),
      skip_reason = "",
      stringsAsFactors = FALSE
    )

    combined_counts[[sid]] <- clusters
  }

  # Combined outputs under out_dir
  summary_df <- do.call("rbind", summary_rows[!vapply(summary_rows, is.null, logical(1))])
  if (is.null(summary_df)) {
    summary_df <- data.frame(
      sample_id = character(),
      n_reads = integer(),
      n_clusters = integer(),
      skip_reason = character(),
      stringsAsFactors = FALSE
    )
  }
  utils::write.table(
    summary_df,
    file = file.path(out_dir, "summary_by_sample.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  combined_table <- if (length(combined_counts)) {
    do.call("rbind", combined_counts)
  } else {
    data.frame(
      sample_id = character(),
      gene_id = character(),
      cluster_id = integer(),
      seq = character(),
      n_reads = integer(),
      n_reads_gene = integer(),
      fraction = numeric(),
      fraction_sample = numeric(),
      method = character(),
      stringsAsFactors = FALSE
    )
  }

  list(
    samples = sample_ids,
    out_dir = out_dir,
    table = combined_table
  )
}

################################################################################
# Internal helpers (Phase C0)
################################################################################

#' @keywords internal
.make_amplicon_stats <- function(sample_id,
                                 n_reads_in,
                                 n_reads_assigned_gene,
                                 n_genes_detected,
                                 n_clusters_total,
                                 skip_reason) {
  data.frame(
    sample_id = sample_id,
    n_reads_in = n_reads_in,
    n_reads_assigned_gene = n_reads_assigned_gene,
    n_reads_unassigned_gene = n_reads_in - n_reads_assigned_gene,
    n_genes_detected = n_genes_detected,
    n_clusters_total = n_clusters_total,
    n_skipped_low_count = as.integer(skip_reason == "low_sample_reads"),
    skip_reason = skip_reason,
    stringsAsFactors = FALSE
  )
}

# Resolve sequences for all samples.
# Uses per-sample FASTQ when available; otherwise buckets from one FASTQ pass.
#' @keywords internal
.resolve_all_sample_sequences <- function(sample_ids,
                                          assignments,
                                          fastq,
                                          sample_fastq_dir = NULL) {
  result <- stats::setNames(vector("list", length(sample_ids)), sample_ids)
  need_fastq <- character()

  for (sid in sample_ids) {
    seqs <- .read_sample_fastq_if_present(sid, sample_fastq_dir)
    if (!is.null(seqs)) {
      result[[sid]] <- seqs
    } else {
      need_fastq <- c(need_fastq, sid)
      result[[sid]] <- character()
    }
  }

  if (length(need_fastq) > 0L) {
    if (is.null(fastq) || length(fastq) < 1L) {
      stop("Either sample_fastq_dir (with per-sample FASTQ) or fastq must be provided.")
    }
    id2sample <- stats::setNames(assignments$sample_id, assignments$read_id)
    buckets <- .bucket_sequences_from_fastq(
      fastq = fastq,
      id2sample = id2sample,
      target_samples = need_fastq
    )
    for (sid in need_fastq) {
      result[[sid]] <- buckets[[sid]] %||% character()
    }
  }

  result
}

#' @keywords internal
.read_sample_fastq_if_present <- function(sample_id, sample_fastq_dir) {
  if (is.null(sample_fastq_dir)) {
    return(NULL)
  }
  cand <- file.path(sample_fastq_dir, paste0(sample_id, ".fq.gz"))
  if (!file.exists(cand)) {
    cand <- file.path(sample_fastq_dir, paste0(sample_id, ".fq"))
  }
  if (file.exists(cand)) {
    return(.read_all_sequences_from_fastq(cand))
  }
  NULL
}

# One pass over FASTQ(s), bucket sequences by sample_id.
#' @keywords internal
.bucket_sequences_from_fastq <- function(fastq, id2sample, target_samples) {
  buckets <- stats::setNames(
    vector("list", length(target_samples)),
    target_samples
  )
  for (sid in target_samples) {
    buckets[[sid]] <- character()
  }

  for (fq in fastq) {
    con <- if (grepl("\\.gz$", fq, ignore.case = TRUE)) {
      gzfile(fq, open = "rt")
    } else {
      file(fq, open = "rt")
    }

    repeat {
      h <- readLines(con, n = 1L)
      if (length(h) == 0L) break
      s <- readLines(con, n = 1L)
      p <- readLines(con, n = 1L)
      q <- readLines(con, n = 1L)
      if (length(s) == 0L || length(p) == 0L || length(q) == 0L) break

      rid <- sub("^@", "", strsplit(h, "\\s+")[[1]][1])
      sid <- unname(id2sample[rid])
      if (length(sid) == 1L && !is.na(sid) && sid %in% target_samples) {
        buckets[[sid]] <- c(buckets[[sid]], s)
      }
    }
    close(con)
  }

  buckets
}

# Null-coalescing helper for lists.
#' @keywords internal
`%||%` <- function(x, y) if (is.null(x)) y else x

# Read all sequences from a (possibly gzipped) FASTQ file.
.read_all_sequences_from_fastq <- function(fastq) {
  con <- if (grepl("\\.gz$", fastq, ignore.case = TRUE)) {
    gzfile(fastq, open = "rt")
  } else {
    file(fastq, open = "rt")
  }
  on.exit(close(con), add = TRUE)

  seqs <- character()
  repeat {
    h <- readLines(con, n = 1L)
    if (length(h) == 0L) break
    s <- readLines(con, n = 1L)
    p <- readLines(con, n = 1L)
    q <- readLines(con, n = 1L)
    if (length(s) == 0L || length(p) == 0L || length(q) == 0L) break
    seqs <- c(seqs, s)
  }
  seqs
}

