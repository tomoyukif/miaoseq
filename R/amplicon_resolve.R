################################################################################
# Step 2: Amplicon resolve
################################################################################

#' Assign reads to genes (Step 2a)
#'
#' Reads demultiplexed FASTQ records, assigns each read to a target gene by
#' primer matching, applies orientation normalization, and optionally labels
#' **trimmed-insert** length outliers against `amplicon_fn` (amplicon width
#' minus primer lengths). Trim failures are `assign_status = "trim_fail"`;
#' length outliers (successful trim, unexpected insert length) are
#' `assign_status = "length_outlier"`. Neither is dropped.
#'
#' @param assignments Path to `assignments.tsv` or a data.frame with at least
#'   `read_id` and `index_pair_id` (legacy tables with only `sample_id` are
#'   accepted and copied to `index_pair_id`).
#' @param out_dir Output directory (typically `{run_dir}/amplicon_assign`).
#' @param fastq Character vector of source FASTQ paths. Required when
#'   per-`index_pair_id` FASTQ files are unavailable.
#' @param sample_fastq_dir Optional directory with `by_sample/{index_pair_id}.fq`
#'   or `.fq.gz`.
#' @param primer_list Required CSV of primer sequences (`primer_id`, `seq`).
#'   `NULL` / missing primers (`no_ref`) are not supported.
#' @param amplicon_fn Optional amplicon reference FASTA used only for
#'   expected-length labeling (sequence widths minus primer lengths). Gene
#'   assignment by amplicon NW fallback has been removed.
#' @param max_primer_edit Maximum edit distance for primer matching / trim.
#' @param end_window Window size (bp) at read ends for primer search.
#' @param length_tolerance Relative length tolerance when `amplicon_fn` is set
#'   (trimmed insert vs amplicon width − F/R primer). Failed trims are
#'   `trim_fail`; length outliers stay `length_outlier`. Neither is removed.
#' @param samples Optional subset of `index_pair_id` to process.
#' @param n_core Number of OpenMP threads for C++ assignment/trimming.
#' @param store_sequences If `FALSE`, write empty `oriented_seq` /
#'   `oriented_qual` (smaller TSV; Editcall can reload FASTQ). Assemble /
#'   Reassess need `TRUE` (default) when reading assignments from disk.
#' @param overwrite If `FALSE`, stop when `out_dir` is non-empty.
#' @param stats_unassign If `TRUE`, write `stats_unassigned.tsv`.
#'
#' @return `GeneAssignResult` list with `samples`, `out_dir`, `table`,
#'   and `stats_unassigned`.
#' @export
#' @importFrom utils read.csv write.table
#' @importFrom BiocGenerics width
#' @importFrom Biostrings DNAString
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings reverseComplement
#' @importFrom Biostrings readDNAStringSet
#' @importFrom Biostrings writeXStringSet
doAssignGenes <- function(
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
) {
  out_dir <- .abs_path(out_dir[[1L]], mustWork = FALSE)
  if (is.character(assignments) && length(assignments) == 1L) {
    assignments <- .abs_path(assignments, mustWork = TRUE)
  }
  if (!is.null(sample_fastq_dir)) {
    sample_fastq_dir <- .abs_path(sample_fastq_dir[[1L]], mustWork = TRUE)
  }
  if (!is.null(fastq)) {
    fastq <- .abs_path(fastq, mustWork = TRUE)
  }
  primer_list <- .require_primer_list(primer_list)
  if (!is.null(amplicon_fn)) {
    amplicon_fn <- .abs_path(amplicon_fn[[1L]], mustWork = TRUE)
  }

  if (!overwrite && dir.exists(out_dir) && length(list.files(out_dir)) > 0L) {
    stop("out_dir already exists and is not empty. Set overwrite = TRUE to overwrite.")
  }
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  out_dir <- .abs_path(out_dir, mustWork = TRUE)

  assn <- .ensure_index_pair_col(.as_assignment_df(assignments))
  if (!all(c("read_id", "index_pair_id") %in% names(assn))) {
    stop("assignments must contain read_id and index_pair_id columns.")
  }
  if (!is.null(samples)) {
    assn <- assn[assn$index_pair_id %in% samples, , drop = FALSE]
  }
  if (nrow(assn) < 1L) {
    warning("No assignments to process in doAssignGenes().")
    empty <- .empty_gene_assign_df()
    utils::write.table(
      empty,
      file = file.path(out_dir, "gene_assignments.tsv"),
      sep = "\t", row.names = FALSE, quote = FALSE
    )
    return(list(
      samples = character(),
      out_dir = out_dir,
      table = empty,
      stats_unassigned = NULL
    ))
  }

  primers <- .parse_primer_pairs(primer_list)
  expected_lengths <- if (!is.null(amplicon_fn)) {
    amplicon_refs <- Biostrings::readDNAStringSet(amplicon_fn)
    amp_w <- stats::setNames(
      as.integer(BiocGenerics::width(amplicon_refs)),
      names(amplicon_refs)
    )
    .expected_insert_lengths(amp_w, primers)
  } else {
    NULL
  }

  index_pair_ids <- sort(unique(assn$index_pair_id))
  pair_reads <- .resolve_all_sample_reads(
    index_pair_ids = index_pair_ids,
    assignments = assn,
    fastq = fastq,
    sample_fastq_dir = sample_fastq_dir
  )
  pair_reads <- .assign_all_sample_reads(
    pair_reads = pair_reads,
    primers = primers,
    end_window = end_window,
    max_primer_edit = max_primer_edit,
    n_core = n_core
  )

  parts <- vector("list", length(index_pair_ids))
  for (i in seq_along(index_pair_ids)) {
    pair_id <- index_pair_ids[[i]]
    reads <- pair_reads[[pair_id]]
    if (is.null(reads) || nrow(reads) < 1L) {
      parts[[i]] <- .empty_gene_assign_df()
      next
    }
    if (!"qual" %in% names(reads)) {
      reads$qual <- rep("", nrow(reads))
    }
    reads$index_pair_id <- pair_id
    oriented <- .orient_reads_with_qual(reads$seq, reads$qual, reads$flipped)
    reads$oriented_seq <- oriented$seq
    reads$oriented_qual <- oriented$qual
    if (!is.null(expected_lengths)) {
      reads <- .label_length_outliers(
        reads = reads,
        primers = primers,
        expected_lengths = expected_lengths,
        length_tolerance = length_tolerance,
        max_primer_edit = max_primer_edit,
        n_core = n_core
      )
    }
    if (!isTRUE(store_sequences)) {
      reads$oriented_seq <- rep("", nrow(reads))
      reads$oriented_qual <- rep("", nrow(reads))
    }
    parts[[i]] <- reads[, c(
      "read_id", "index_pair_id", "gene_id",
      "assign_status", "flipped", "total_edit", "f_edit", "r_edit",
      "oriented_seq", "oriented_qual"
    )]
  }

  gene_assignments <- do.call("rbind", parts[!vapply(parts, is.null, logical(1))])
  if (is.null(gene_assignments)) {
    gene_assignments <- .empty_gene_assign_df()
  }
  rownames(gene_assignments) <- NULL
  utils::write.table(
    gene_assignments,
    file = file.path(out_dir, "gene_assignments.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )

  pkg_ver <- tryCatch(
    as.character(utils::packageVersion("miaoseq")),
    error = function(e) NA_character_
  )
  write(
    paste0(
      "r_version=", R.version.string, "\n",
      "miaoseq_version=", pkg_ver, "\n",
      "primer_list=", primer_list, "\n",
      "amplicon_fn=", if (is.null(amplicon_fn)) "" else amplicon_fn, "\n",
      "max_primer_edit=", max_primer_edit, "\n",
      "end_window=", end_window, "\n",
      "length_tolerance=", length_tolerance, "\n",
      "store_sequences=", store_sequences, "\n",
      "n_core=", n_core, "\n",
      "n_reads=", nrow(gene_assignments), "\n"
    ),
    file = file.path(out_dir, "run_stats.txt")
  )

  stats_unassigned_df <- NULL
  if (isTRUE(stats_unassign)) {
    stats_unassigned_df <- .stats_unassign_amplicon(gene_assignments)
    utils::write.table(
      stats_unassigned_df,
      file = file.path(out_dir, "stats_unassigned.tsv"),
      sep = "\t", row.names = FALSE, quote = FALSE
    )
  }

  list(
    samples = index_pair_ids,
    out_dir = out_dir,
    table = gene_assignments,
    stats_unassigned = stats_unassigned_df
  )
}

#' Assemble amplicons from gene-assigned reads (Step 2b)
#'
#' Consumes `GeneAssignResult` (from [doAssignGenes()]) or `gene_assignments.tsv`
#' and performs per `index_pair_id` x gene clustering / consensus.
#'
#' @param gene_assign `GeneAssignResult` list, `gene_assignments.tsv` path, or
#'   a data.frame.
#' @param out_dir Output directory for per-`index_pair_id` consensus and cluster tables.
#' @param primer_list Required primer CSV used for insert trimming before
#'   clustering. `NULL` is not supported.
#' @param method One of `"both"`, `"cluster"`, or `"consensus"`.
#' @param min_reads Minimum reads per `index_pair_id` x gene bucket to process.
#' @param min_cluster_reads Minimum reads per cluster to retain (minor-cluster filter).
#' @param max_clusters Maximum clusters per sample x gene. Default `Inf` (no
#'   top-N truncation). Use a finite integer to keep only the largest clusters.
#' @param min_cluster_identity Minimum sequence identity for vsearch clustering
#'   (`--id`, `--iddef 2`). Default `0.95`.
#' @param max_primer_edit Maximum edit distance for primer trim.
#' @param cluster_backend Clustering method. Fixed to `"vsearch"` (Phase K).
#' @param assembly_backend Reconstruction backend (`"cluster"` or `"overlap_graph"`).
#' @param overlap_min_identity Identity threshold for overlap-graph backend.
#' @param strict_end_trim If `TRUE` (default), cluster only successfully
#'   trimmed inserts; trim failures become U1. If `FALSE`, mix oriented
#'   full-length failed reads into the same cluster input (not recommended).
#' @param min_cluster_purity Minimum fraction of cluster members that share the
#'   most common exact insert sequence. Clusters below this purity are
#'   discarded. Default `NULL` (disabled; purity is still reported). Set to
#'   e.g. `0.8` for strict filtering.
#' @param overwrite If `FALSE`, stop when `out_dir` is non-empty.
#' @param n_core OpenMP threads for primer trim.
#'
#' @return A list with `samples`, `out_dir`, and assembled `table`
#'   (data.frame). Key columns in `table`: `fraction` = reads in cluster /
#'   reads in same gene bucket after clustering; `fraction_bucket` = reads in
#'   cluster / total reads entering the **index_pair_id x gene** bucket after primer
#'   trim (NOT total demultiplexed reads for the well).
#' @export
doAssembleAmplicons <- function(
  gene_assign,
  out_dir,
  primer_list,

  method = c("both", "cluster", "consensus"),
  min_reads = 5L,
  min_cluster_reads = 5L,
  max_clusters = Inf,
  min_cluster_identity = 0.95,
  max_primer_edit = 5L,
  cluster_backend = "vsearch",
  assembly_backend = c("cluster", "overlap_graph"),
  overlap_min_identity = 0.90,
  strict_end_trim = TRUE,
  min_cluster_purity = NULL,
  overwrite = FALSE,
  n_core = 1L
) {
  method <- match.arg(method)
  cluster_backend <- match.arg(cluster_backend, choices = "vsearch")
  if (!is.numeric(min_cluster_identity) || length(min_cluster_identity) != 1L ||
      is.na(min_cluster_identity) || min_cluster_identity <= 0 ||
      min_cluster_identity > 1) {
    stop("min_cluster_identity must be a single numeric in (0, 1].")
  }
  if (!(is.infinite(max_clusters) ||
          (is.numeric(max_clusters) && length(max_clusters) == 1L &&
             !is.na(max_clusters) && max_clusters >= 1))) {
    stop("max_clusters must be Inf or a number >= 1.")
  }
  if (!is.null(min_cluster_purity) &&
      (!is.numeric(min_cluster_purity) ||
         length(min_cluster_purity) != 1L ||
         is.na(min_cluster_purity) ||
         min_cluster_purity < 0 ||
         min_cluster_purity > 1)) {
    stop("min_cluster_purity must be NULL or a single numeric in [0, 1].")
  }
  assembly_backend <- match.arg(assembly_backend)
  primer_list <- .require_primer_list(primer_list)
  if (identical(assembly_backend, "cluster") &&
      method %in% c("both", "consensus") &&
      !nzchar(Sys.which("abpoa"))) {
    stop("Consensus assembly requires the abpoa binary on PATH.")
  }
  if (identical(assembly_backend, "cluster") &&
      !nzchar(Sys.which("vsearch"))) {
    stop("Clustering requires the vsearch binary on PATH.")
  }
  if (!overwrite && dir.exists(out_dir) && length(list.files(out_dir)) > 0L) {
    stop("out_dir already exists and is not empty. Set overwrite = TRUE to overwrite.")
  }
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  ga <- .as_gene_assign_df(gene_assign)
  if (nrow(ga) > 0L &&
      all(is.na(ga$oriented_seq) | !nzchar(as.character(ga$oriented_seq)))) {
    stop(
      "oriented_seq is empty in gene_assign. ",
      "doAssignGenes(store_sequences=FALSE) cannot be used with Assemble; ",
      "re-run assignment with store_sequences=TRUE.",
      call. = FALSE
    )
  }
  index_pair_ids <- sort(unique(ga$index_pair_id))
  primers <- .parse_primer_pairs(primer_list)

  combined_counts <- list()
  combined_members <- list()
  summary_rows <- vector("list", length(index_pair_ids))
  for (i in seq_along(index_pair_ids)) {
    pair_id <- index_pair_ids[[i]]
    pair_dir <- file.path(out_dir, pair_id)
    if (!dir.exists(pair_dir)) {
      dir.create(pair_dir, recursive = TRUE, showWarnings = FALSE)
    }
    reads <- ga[ga$index_pair_id == pair_id, , drop = FALSE]
    n_in <- nrow(reads)
    if (n_in < 1L) {
      .write_sample_skip(pair_dir, pair_id, 0L, 0L, "empty_sample")
      summary_rows[[i]] <- .make_summary_row(pair_id, 0L, 0L, 0L, "empty_sample")
      next
    }
    reads <- reads[.is_processable_gene_row(reads$gene_id, reads$assign_status),
                   , drop = FALSE]
    gene_ids <- unique(reads$gene_id)
    if (length(gene_ids) < 1L) {
      .write_sample_skip(pair_dir, pair_id, n_in, 0L, "no_gene_assigned")
      summary_rows[[i]] <- .make_summary_row(pair_id, n_in, 0L, 0L, "no_gene_assigned")
      next
    }

    pair_clusters <- list()
    pair_members <- list()
    pair_unassigned <- list()
    n_assigned <- 0L
    n_clusters_total <- 0L
    n_genes_low_reads <- 0L
    n_genes_no_clusters <- 0L
    for (gid in gene_ids) {
      gene_reads <- reads[reads$gene_id == gid, , drop = FALSE]
      if (nrow(gene_reads) < min_reads) {
        n_genes_low_reads <- n_genes_low_reads + 1L
        pair_unassigned[[length(pair_unassigned) + 1L]] <- .unassigned_rows(
          pair_id, gid, gene_reads$read_id, "low_gene_reads"
        )
        next
      }
      prepared <- .prepare_gene_seqs_for_cluster(
        gene_reads = gene_reads,
        primers = primers,
        gene_id = gid,
        max_primer_edit = max_primer_edit,
        n_core = n_core,
        strict_end_trim = strict_end_trim
      )
      # Reads dropped by trim / empty sequence.
      dropped_ids <- setdiff(
        as.character(gene_reads$read_id),
        as.character(prepared$read_id)
      )
      if (length(dropped_ids) > 0L) {
        pair_unassigned[[length(pair_unassigned) + 1L]] <- .unassigned_rows(
          pair_id, gid, dropped_ids, "trim_or_empty"
        )
      }
      if (length(prepared$seqs) < min_reads) {
        n_genes_no_clusters <- n_genes_no_clusters + 1L
        pair_unassigned[[length(pair_unassigned) + 1L]] <- .unassigned_rows(
          pair_id, gid, prepared$read_id, "insufficient_after_trim"
        )
        next
      }
      assembled <- if (identical(assembly_backend, "overlap_graph")) {
        .assemble_overlap_graph_reads(
          seqs = prepared$seqs,
          gene_id = gid,
          index_pair_id = pair_id,
          min_cluster_reads = min_cluster_reads,
          max_clusters = max_clusters,
          min_cluster_identity = min_cluster_identity,
          overlap_min_identity = overlap_min_identity,
          min_cluster_purity = min_cluster_purity,
          read_id = prepared$read_id
        )
      } else {
        .cluster_and_consensus_reads(
          seqs = prepared$seqs,
          gene_id = gid,
          index_pair_id = pair_id,
          method = method,
          min_cluster_reads = min_cluster_reads,
          max_clusters = max_clusters,
          min_cluster_identity = min_cluster_identity,
          cluster_backend = cluster_backend,
          min_cluster_purity = min_cluster_purity,
          read_id = prepared$read_id
        )
      }
      clusters <- assembled$clusters
      if (!is.null(assembled$unassigned) && nrow(assembled$unassigned) > 0L) {
        pair_unassigned[[length(pair_unassigned) + 1L]] <- assembled$unassigned
      }
      if (!is.null(assembled$members) && nrow(assembled$members) > 0L) {
        pair_members[[length(pair_members) + 1L]] <- assembled$members
      }
      if (is.null(clusters) || nrow(clusters) < 1L) {
        n_genes_no_clusters <- n_genes_no_clusters + 1L
        next
      }
      pair_clusters[[gid]] <- clusters
      n_assigned <- n_assigned + sum(clusters$n_reads)
      n_clusters_total <- n_clusters_total + nrow(clusters)
    }

    unassigned_df <- if (length(pair_unassigned) > 0L) {
      do.call(rbind, pair_unassigned)
    } else {
      .empty_unassigned_to_cluster()
    }
    if (!is.null(unassigned_df) && nrow(unassigned_df) > 0L) {
      rownames(unassigned_df) <- NULL
      unassigned_df <- .relabel_unit_col(unassigned_df)
      utils::write.table(
        unassigned_df,
        file = file.path(pair_dir, "unassigned_to_cluster.tsv"),
        sep = "\t", row.names = FALSE, quote = FALSE
      )
    }

    if (length(pair_clusters) < 1L) {
      skip_reason <- if (n_genes_low_reads > 0L) "low_gene_reads" else "no_clusters"
      .write_sample_skip(pair_dir, pair_id, n_in, n_genes_low_reads, skip_reason)
      summary_rows[[i]] <- .make_summary_row(pair_id, n_in, 0L, 0L, skip_reason)
      next
    }

    clusters <- do.call(rbind, pair_clusters)
    rownames(clusters) <- NULL
    dna <- Biostrings::DNAStringSet(clusters$seq)
    clusters <- .relabel_unit_col(clusters)
    names(dna) <- sprintf(
      "%s|cluster_%d;size=%d;index_pair=%s",
      clusters$gene_id, clusters$cluster_id, clusters$n_reads, clusters$index_pair_id
    )
    Biostrings::writeXStringSet(dna, filepath = file.path(pair_dir, "consensus.fasta"))
    if (method %in% c("both", "cluster") && length(pair_members) > 0L) {
      members_df <- do.call(rbind, pair_members)
      rownames(members_df) <- NULL
      mem_dna <- Biostrings::DNAStringSet(members_df$seq)
      names(mem_dna) <- paste(members_df$read_id, members_df$cluster_id)
      Biostrings::writeXStringSet(
        mem_dna,
        filepath = file.path(pair_dir, "clusters.fasta")
      )
      combined_members[[pair_id]] <- members_df
    }
    count_cols <- c(
      "index_pair_id", "gene_id", "cluster_id",
      "n_reads", "fraction", "fraction_bucket",
      "n_reads_gene", "method", "cluster_purity"
    )
    count_cols <- count_cols[count_cols %in% names(clusters)]
    utils::write.table(
      clusters[, count_cols, drop = FALSE],
      file = file.path(pair_dir, "cluster_counts.tsv"),
      sep = "\t", row.names = FALSE, quote = FALSE
    )
    stats_df <- .make_amplicon_stats(
      index_pair_id = pair_id,
      n_reads_in = n_in,
      n_reads_assigned_gene = n_assigned,
      n_genes_detected = length(unique(clusters$gene_id)),
      n_clusters_total = n_clusters_total,
      n_skipped_low_count = n_genes_low_reads,
      skip_reason = ""
    )
    utils::write.table(
      stats_df,
      file = file.path(pair_dir, "stats.tsv"),
      sep = "\t", row.names = FALSE, quote = FALSE
    )
    summary_rows[[i]] <- .make_summary_row(
      pair_id, n_in, n_clusters_total, length(unique(clusters$gene_id)), ""
    )
    combined_counts[[pair_id]] <- clusters
  }

  summary_df <- do.call("rbind", summary_rows[!vapply(summary_rows, is.null, logical(1))])
  if (is.null(summary_df)) {
    summary_df <- data.frame(
      index_pair_id = character(),
      n_reads = integer(),
      n_clusters = integer(),
      n_genes = integer(),
      skip_reason = character(),
      stringsAsFactors = FALSE
    )
  }
  utils::write.table(
    summary_df,
    file = file.path(out_dir, "summary_by_sample.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )

  combined_table <- if (length(combined_counts)) {
    .relabel_unit_col(do.call("rbind", combined_counts))
  } else {
    data.frame(
      index_pair_id = character(),
      gene_id = character(),
      cluster_id = integer(),
      seq = character(),
      n_reads = integer(),
      n_reads_gene = integer(),
      fraction = numeric(),
      fraction_bucket = numeric(),
      method = character(),
      stringsAsFactors = FALSE
    )
  }

  tools_needed <- c("vsearch", "abpoa")
  versions <- .tool_versions(tools_needed)
  pkg_ver <- tryCatch(
    as.character(utils::packageVersion("miaoseq")),
    error = function(e) NA_character_
  )
  write(
    paste0(
      "r_version=", R.version.string, "\n",
      "miaoseq_version=", pkg_ver, "\n",
      "cluster_backend=", cluster_backend, "\n",
      "assembly_backend=", assembly_backend, "\n",
      "min_cluster_identity=", min_cluster_identity, "\n",
      "max_clusters=", max_clusters, "\n",
      "min_cluster_reads=", min_cluster_reads, "\n",
      "min_reads=", min_reads, "\n",
      "min_cluster_purity=",
      if (is.null(min_cluster_purity)) "NULL" else min_cluster_purity, "\n",
      "strict_end_trim=", strict_end_trim, "\n",
      "n_samples=", length(index_pair_ids), "\n",
      "n_clusters=", nrow(combined_table), "\n",
      paste(paste0("version_", names(versions), "=", versions), collapse = "\n"),
      "\n"
    ),
    file = file.path(out_dir, "run_stats.txt")
  )

  list(
    samples = index_pair_ids,
    out_dir = out_dir,
    table = combined_table,
    gene_assignments = ga
  )
}

#' @export
doAmpliconResolve <- doAssembleAmplicons

################################################################################
# Internal helpers
################################################################################

#' @keywords internal
.require_primer_list <- function(primer_list) {
  if (is.null(primer_list) ||
      (is.character(primer_list) &&
         (length(primer_list) < 1L || !nzchar(as.character(primer_list[[1L]]))))) {
    stop(
      "primer_list is required. NULL / no_ref / unknown-gene assignment is not supported.",
      call. = FALSE
    )
  }
  .abs_path(primer_list[[1L]], mustWork = TRUE)
}

#' @keywords internal
.make_summary_row <- function(index_pair_id, n_reads, n_clusters, n_genes, skip_reason) {
  data.frame(
    index_pair_id = index_pair_id,
    n_reads = n_reads,
    n_clusters = n_clusters,
    n_genes = n_genes,
    skip_reason = skip_reason,
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
.write_sample_skip <- function(pair_dir, index_pair_id, n_in,
                               n_skipped_low_count = 0L,
                               skip_reason) {
  stats_df <- .make_amplicon_stats(
    index_pair_id = index_pair_id,
    n_reads_in = n_in,
    n_reads_assigned_gene = 0L,
    n_genes_detected = 0L,
    n_clusters_total = 0L,
    n_skipped_low_count = n_skipped_low_count,
    skip_reason = skip_reason
  )
  utils::write.table(
    stats_df,
    file = file.path(pair_dir, "stats.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  invisible(stats_df)
}

#' @keywords internal
.make_amplicon_stats <- function(index_pair_id,
                                 n_reads_in,
                                 n_reads_assigned_gene,
                                 n_genes_detected,
                                 n_clusters_total,
                                 n_skipped_low_count = 0L,
                                 skip_reason) {
  data.frame(
    index_pair_id = index_pair_id,
    n_reads_in = n_reads_in,
    n_reads_assigned_gene = n_reads_assigned_gene,
    n_reads_unassigned_gene = n_reads_in - n_reads_assigned_gene,
    n_genes_detected = n_genes_detected,
    n_clusters_total = n_clusters_total,
    n_skipped_low_count = as.integer(n_skipped_low_count),
    skip_reason = skip_reason,
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
.empty_gene_assign_df <- function() {
  data.frame(
    read_id = character(),
    index_pair_id = character(),
    gene_id = character(),
    assign_status = character(),
    flipped = logical(),
    total_edit = integer(),
    f_edit = integer(),
    r_edit = integer(),
    oriented_seq = character(),
    oriented_qual = character(),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
.as_gene_assign_df <- function(x) {
  if (is.list(x) && !is.null(x$gene_assignments) &&
      is.data.frame(x$gene_assignments)) {
    out <- x$gene_assignments
  } else if (is.list(x) && !is.null(x$table)) {
    out <- x$table
  } else if (is.character(x) && length(x) == 1L && file.exists(x)) {
    p <- x
    if (dir.exists(p)) {
      p <- file.path(p, "gene_assignments.tsv")
    }
    if (!file.exists(p)) {
      stop("gene_assign path does not contain gene_assignments.tsv: ", x)
    }
    out <- utils::read.delim(p, stringsAsFactors = FALSE)
  } else if (is.data.frame(x)) {
    out <- x
  } else {
    stop("gene_assign must be GeneAssignResult, path to gene_assignments.tsv, or data.frame.")
  }
  out <- .ensure_index_pair_col(out)
  req <- c("read_id", "index_pair_id", "gene_id", "assign_status", "oriented_seq")
  if (!all(req %in% names(out))) {
    stop("gene_assign must contain: ", paste(req, collapse = ", "))
  }
  if (!"oriented_qual" %in% names(out)) {
    out$oriented_qual <- rep("", nrow(out))
  }
  out
}

#' @keywords internal
.parse_primer_pairs <- function(primer_list) {
  primers <- .read_list_csv(
    primer_list,
    min_cols = 2L,
    header_tokens = c("primer_id", "primer", "id", "sequence", "seq", "gene")
  )
  if (ncol(primers) < 2L) {
    stop("primer_list must have at least two columns: primer_id, sequence.")
  }
  names(primers)[1:2] <- c("primer_id", "seq")
  primers$seq <- toupper(gsub("[^ACGT]", "", primers$seq))

  f_idx <- grepl("_F$", primers$primer_id)
  r_idx <- grepl("_R$", primers$primer_id)
  f_primers <- primers[f_idx, , drop = FALSE]
  r_primers <- primers[r_idx, , drop = FALSE]
  f_primers$gene_id <- sub("_F$", "", f_primers$primer_id)
  r_primers$gene_id <- sub("_R$", "", r_primers$primer_id)

  pairs <- merge(
    f_primers[, c("gene_id", "seq")],
    r_primers[, c("gene_id", "seq")],
    by = "gene_id",
    suffixes = c(".f", ".r")
  )
  if (nrow(pairs) < 1L) {
    stop("primer_list must contain at least one _F/_R primer pair.")
  }
  missing_r <- setdiff(f_primers$gene_id, r_primers$gene_id)
  missing_f <- setdiff(r_primers$gene_id, f_primers$gene_id)
  if (length(missing_r) || length(missing_f)) {
    stop(
      "primer_list has unpaired genes: missing _R for ",
      paste(missing_r, collapse = ", "),
      "; missing _F for ",
      paste(missing_f, collapse = ", ")
    )
  }

  pairs <- pairs[order(pairs$gene_id), , drop = FALSE]
  rownames(pairs) <- NULL
  pairs
}

#' @keywords internal
.assign_all_sample_reads <- function(pair_reads, primers,
                                     end_window, max_primer_edit,
                                     n_core = 1L) {
  pair_ids <- names(pair_reads)
  if (length(pair_ids) < 1L) return(pair_reads)

  n_per <- vapply(pair_reads, function(x) {
    if (is.null(x) || !is.data.frame(x)) 0L else nrow(x)
  }, integer(1))

  all_ids <- character(sum(n_per))
  all_seqs <- character(sum(n_per))
  all_quals <- character(sum(n_per))
  offset <- 0L
  for (i in seq_along(pair_ids)) {
    n <- n_per[[i]]
    if (n < 1L) next
    idx <- offset + seq_len(n)
    part <- pair_reads[[pair_ids[[i]]]]
    all_ids[idx] <- as.character(part$read_id)
    all_seqs[idx] <- as.character(part$seq)
    all_quals[idx] <- if ("qual" %in% names(part)) {
      as.character(part$qual)
    } else {
      rep("", n)
    }
    offset <- offset + n
  }

  gene_id <- rep(NA_character_, length(all_seqs))
  flipped <- rep(FALSE, length(all_seqs))
  total_edit <- rep(NA_integer_, length(all_seqs))
  f_edit <- rep(NA_integer_, length(all_seqs))
  r_edit <- rep(NA_integer_, length(all_seqs))
  assign_status <- rep(NA_character_, length(all_seqs))

  if (!is.null(primers) && length(all_seqs) > 0L) {
    res <- assign_genes_primers_cpp(
      seqs = all_seqs,
      gene_ids = primers$gene_id,
      f_primers = primers$seq.f,
      r_primers = primers$seq.r,
      end_window = as.integer(end_window),
      max_primer_edit = as.integer(max_primer_edit),
      n_core = as.integer(max(1L, n_core))
    )
    gene_id <- res$gene_id
    flipped <- res$flipped
    total_edit <- res$total_edit
    f_edit <- res$f_edit
    r_edit <- res$r_edit
    assign_status <- res$status
  }

  offset <- 0L
  for (i in seq_along(pair_ids)) {
    n <- n_per[[i]]
    if (n < 1L) {
      pair_reads[[pair_ids[[i]]]] <- data.frame(
        read_id = character(),
        seq = character(),
        qual = character(),
        gene_id = character(),
        flipped = logical(),
        total_edit = integer(),
        f_edit = integer(),
        r_edit = integer(),
        assign_status = character(),
        stringsAsFactors = FALSE
      )
      next
    }
    idx <- offset + seq_len(n)
    pair_reads[[pair_ids[[i]]]] <- data.frame(
      read_id = all_ids[idx],
      seq = all_seqs[idx],
      qual = all_quals[idx],
      gene_id = gene_id[idx],
      flipped = flipped[idx],
      total_edit = total_edit[idx],
      f_edit = f_edit[idx],
      r_edit = r_edit[idx],
      assign_status = assign_status[idx],
      stringsAsFactors = FALSE
    )
    offset <- offset + n
  }
  pair_reads
}

#' @keywords internal
.assign_reads_to_genes <- function(reads, primers,
                                   end_window, max_primer_edit,
                                   n_core = 1L) {
  # Compatibility wrapper: single-sample path via batch helper.
  tmp <- list(x = reads)
  names(tmp) <- "x"
  out <- .assign_all_sample_reads(
    pair_reads = tmp,
    primers = primers,
    end_window = end_window,
    max_primer_edit = max_primer_edit,
    n_core = n_core
  )[["x"]]
  out
}

#' @keywords internal
.orient_read_seqs <- function(seqs, flipped) {
  if (length(seqs) < 1L) return(character())
  as.character(reverse_complement_seqs_cpp(seqs, flipped))
}

#' @keywords internal
.orient_reads_with_qual <- function(seqs, quals, flipped) {
  if (length(seqs) < 1L) {
    return(list(seq = character(), qual = character()))
  }
  flipped <- as.logical(flipped)
  if (length(flipped) == 1L) flipped <- rep(flipped, length(seqs))
  out_seq <- as.character(reverse_complement_seqs_cpp(seqs, flipped))
  if (is.null(quals) || length(quals) != length(seqs)) {
    quals <- rep("", length(seqs))
  }
  out_qual <- as.character(reverse_quals_cpp(as.character(quals), flipped))
  list(seq = out_seq, qual = out_qual)
}

#' @keywords internal
.prepare_gene_seqs_for_cluster <- function(gene_reads, primers, gene_id,
                                           max_primer_edit, n_core,
                                           strict_end_trim = TRUE) {
  seqs <- as.character(gene_reads$oriented_seq)
  read_ids <- as.character(gene_reads$read_id)
  if (length(seqs) > 0L && (all(is.na(seqs) | !nzchar(seqs)))) {
    stop(
      "oriented_seq is empty for gene ", gene_id,
      ". doAssignGenes(store_sequences=FALSE) cannot be passed to Assemble; ",
      "re-run assignment with store_sequences=TRUE.",
      call. = FALSE
    )
  }

  if (is.null(primers)) {
    keep <- nzchar(seqs) & !is.na(seqs)
    return(list(
      seqs = seqs[keep],
      read_id = read_ids[keep]
    ))
  }

  fp <- primers$seq.f[primers$gene_id == gene_id]
  rp <- primers$seq.r[primers$gene_id == gene_id]
  trimmed <- trim_amplicon_insert_cpp(
    seqs,
    f_primer = fp,
    r_primer = rp,
    max_edit = as.integer(max_primer_edit),
    n_core = as.integer(max(1L, n_core))
  )
  ok <- !is.na(trimmed$start) & nzchar(trimmed$seq)
  use_seqs <- as.character(trimmed$seq)

  for (i in seq_along(seqs)) {
    if (!isTRUE(ok[[i]])) {
      if (isTRUE(strict_end_trim)) {
        use_seqs[[i]] <- ""
      } else {
        use_seqs[[i]] <- seqs[[i]]
      }
    }
  }

  keep <- nzchar(use_seqs)
  list(
    seqs = use_seqs[keep],
    read_id = read_ids[keep]
  )
}

#' @keywords internal
.resolve_all_sample_reads <- function(index_pair_ids, assignments, fastq,
                                      sample_fastq_dir = NULL) {
  result <- stats::setNames(vector("list", length(index_pair_ids)), index_pair_ids)
  need_fastq <- character()

  for (pair_id in index_pair_ids) {
    reads <- .read_sample_reads_if_present(pair_id, sample_fastq_dir)
    if (!is.null(reads)) {
      result[[pair_id]] <- reads
    } else {
      need_fastq <- c(need_fastq, pair_id)
      result[[pair_id]] <- data.frame(
        read_id = character(),
        seq = character(),
        qual = character(),
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(need_fastq) > 0L) {
    if (is.null(fastq) || length(fastq) < 1L) {
      stop("Either sample_fastq_dir (with per-sample FASTQ) or fastq must be provided.")
    }
    buckets <- .bucket_reads_from_fastq(
      fastq = fastq,
      assignments = assignments,
      target_pairs = need_fastq
    )
    for (pair_id in need_fastq) {
      result[[pair_id]] <- buckets[[pair_id]] %||% data.frame(
        read_id = character(),
        seq = character(),
        qual = character(),
        stringsAsFactors = FALSE
      )
    }
  }

  result
}

#' @keywords internal
.read_sample_reads_if_present <- function(index_pair_id, sample_fastq_dir) {
  if (is.null(sample_fastq_dir)) return(NULL)
  sample_fastq_dir <- .abs_path(sample_fastq_dir[[1L]], mustWork = FALSE)
  safe_id <- .safe_filename(index_pair_id[[1L]])
  cand <- file.path(sample_fastq_dir, paste0(safe_id, ".fq.gz"))
  if (!file.exists(cand)) {
    cand <- file.path(sample_fastq_dir, paste0(safe_id, ".fq"))
  }
  if (file.exists(cand)) {
    return(.read_all_reads_from_fastq(.abs_path(cand, mustWork = TRUE)))
  }
  NULL
}

#' @keywords internal
.bucket_reads_from_fastq <- function(fastq, assignments, target_pairs) {
  raw <- bucket_fastq_assignments_cpp(
    fastq_files = .abs_path(fastq, mustWork = TRUE),
    read_ids = as.character(assignments$read_id),
    sample_ids = as.character(assignments$index_pair_id),
    target_samples = as.character(target_pairs)
  )
  out <- stats::setNames(vector("list", length(target_pairs)), target_pairs)
  for (pair_id in target_pairs) {
    part <- raw[[pair_id]]
    out[[pair_id]] <- data.frame(
      read_id = as.character(part$read_id),
      seq = as.character(part$seq),
      qual = as.character(part$qual),
      stringsAsFactors = FALSE
    )
  }
  out
}

#' @keywords internal
.read_all_reads_from_fastq <- function(fastq) {
  raw <- read_fastq_seqs_cpp(.abs_path(fastq[[1]], mustWork = TRUE))
  data.frame(
    read_id = as.character(raw$read_id),
    seq = as.character(raw$seq),
    qual = as.character(raw$qual),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Capture external tool versions for provenance logging
#' @param tools Character vector of tool names to query (e.g. "vsearch", "abpoa").
#' @return Named character vector of version strings (NA if tool unavailable).
#' @keywords internal
.tool_versions <- function(tools = c("vsearch", "abpoa", "blastn", "mmseqs")) {
  vapply(tools, function(t) {
    bin <- Sys.which(t)
    if (!nzchar(bin)) return(NA_character_)
    tryCatch({
      out <- suppressWarnings(system2(bin, "--version", stdout = TRUE, stderr = TRUE))
      out <- out[nzchar(out)]
      if (length(out) < 1L) return("unknown")
      trimws(out[[1L]])
    }, error = function(e) "unknown")
  }, character(1))
}

#' @keywords internal
.assign_genes_by_amplicon_ref <- function(seqs, amplicon_refs, max_edit,
                                          n_core = 1L) {
  stop(
    "Amplicon-ref gene-assignment fallback has been removed. ",
    "Use primer_list for gene assignment.",
    call. = FALSE
  )
}

#' @keywords internal
.processable_assign_statuses <- function() {
  c("assigned", "length_outlier", "trim_fail")
}

#' @keywords internal
.is_processable_gene_row <- function(gene_id, assign_status) {
  if (is.null(assign_status)) {
    stop(
      "assign_status column is missing. Re-run doAssignGenes() or add the column.",
      call. = FALSE
    )
  }
  !is.na(gene_id) & nzchar(as.character(gene_id)) &
    as.character(assign_status) %in% .processable_assign_statuses()
}

#' @keywords internal
.expected_insert_lengths <- function(amplicon_widths, primers) {
  genes <- names(amplicon_widths)
  out <- as.integer(amplicon_widths)
  names(out) <- genes
  if (is.null(primers) || nrow(primers) < 1L) {
    return(out)
  }
  for (i in seq_along(genes)) {
    gid <- genes[[i]]
    fp <- primers$seq.f[primers$gene_id == gid]
    rp <- primers$seq.r[primers$gene_id == gid]
    f_len <- if (length(fp) >= 1L) nchar(fp[[1]]) else 0L
    r_len <- if (length(rp) >= 1L) nchar(rp[[1]]) else 0L
    out[[i]] <- as.integer(amplicon_widths[[i]] - f_len - r_len)
  }
  out
}

#' @keywords internal
.label_length_outliers <- function(reads,
                                   primers,
                                   expected_lengths,
                                   length_tolerance,
                                   max_primer_edit,
                                   n_core = 1L) {
  if (nrow(reads) < 1L) return(reads)
  assigned <- which(
    !is.na(reads$gene_id) &
      nzchar(as.character(reads$gene_id)) &
      reads$assign_status == "assigned"
  )
  if (length(assigned) < 1L) return(reads)
  by_gene <- split(assigned, as.character(reads$gene_id[assigned]))
  for (gid in names(by_gene)) {
    exp_len <- expected_lengths[[gid]]
    if (is.null(exp_len) || is.na(exp_len) || exp_len < 1L) next
    fp <- primers$seq.f[primers$gene_id == gid]
    rp <- primers$seq.r[primers$gene_id == gid]
    if (length(fp) < 1L || length(rp) < 1L) next
    idx <- by_gene[[gid]]
    trimmed <- trim_amplicon_insert_cpp(
      as.character(reads$oriented_seq[idx]),
      f_primer = fp[[1]],
      r_primer = rp[[1]],
      max_edit = as.integer(max_primer_edit),
      n_core = as.integer(max(1L, n_core))
    )
    tseq <- as.character(trimmed$seq)
    obs <- nchar(tseq)
    trim_fail <- is.na(trimmed$start) | !nzchar(tseq)
    len_bad <- !trim_fail & (
      obs < exp_len * (1 - length_tolerance) |
        obs > exp_len * (1 + length_tolerance)
    )
    if (any(trim_fail)) {
      reads$assign_status[idx[trim_fail]] <- "trim_fail"
    }
    if (any(len_bad)) {
      reads$assign_status[idx[len_bad]] <- "length_outlier"
    }
  }
  reads
}

#' @keywords internal
.consensus_via_abpoa <- function(members) {
  if (length(members) < 1L) return("")
  if (length(members) == 1L) return(members[[1]])
  abpoa <- Sys.which("abpoa")
  if (!nzchar(abpoa)) {
    stop("abpoa binary not found on PATH.")
  }
  tmp <- tempfile(fileext = ".fa")
  on.exit(unlink(tmp), add = TRUE)
  con <- file(tmp, open = "wt")
  for (i in seq_along(members)) {
    writeLines(sprintf(">r%d", i), con)
    writeLines(members[[i]], con)
  }
  close(con)
  out <- system2(abpoa, args = c("-r", "0", "-m", "0", tmp),
                 stdout = TRUE, stderr = FALSE)
  status <- attr(out, "status")
  if (!is.null(status) && !identical(as.integer(status), 0L)) {
    stop("abpoa failed (exit ", status, ")")
  }
  seq_lines <- out[!grepl("^>", out) & nzchar(out)]
  if (length(seq_lines) < 1L) {
    stop("abpoa produced no consensus sequence.")
  }
  paste(seq_lines, collapse = "")
}

#' @keywords internal
.vsearch_id_threshold <- function(min_cluster_identity) {
  id <- as.numeric(min_cluster_identity)
  # Pass through API range (0, 1]; no extra floor (low-id sweeps need < 0.5).
  max(1e-6, min(1.0, id))
}

#' Convert absolute edit budget to vsearch-style identity using median length.
#' Used by Reassess when `min_identity` is not supplied.
#' @keywords internal
.vsearch_id_from_edit <- function(seqs, max_edit) {
  lens <- nchar(as.character(seqs))
  lens <- lens[is.finite(lens) & lens > 0]
  med <- if (length(lens)) stats::median(lens) else 1500
  if (!is.finite(med) || med < 1) med <- 1500
  id <- 1 - as.numeric(max_edit) / as.numeric(med)
  max(0.50, min(0.99, id))
}

#' @keywords internal
.cluster_ids_vsearch <- function(seqs, min_cluster_identity, max_clusters) {
  n <- length(seqs)
  if (n < 1L) return(integer())
  vsearch <- Sys.which("vsearch")
  if (!nzchar(vsearch)) stop("vsearch binary not found on PATH.")

  tmp_dir <- tempfile("vsearch_")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
  fa <- file.path(tmp_dir, "in.fa")
  uc <- file.path(tmp_dir, "out.uc")

  con <- file(fa, open = "wt")
  for (i in seq_len(n)) {
    writeLines(sprintf(">s%d", i - 1L), con)
    writeLines(seqs[[i]], con)
  }
  close(con)

  id <- .vsearch_id_threshold(min_cluster_identity)
  args <- c(
    "--cluster_fast", fa,
    "--id", format(id, digits = 6, trim = TRUE),
    "--strand", "both",
    "--iddef", "2",
    "--uc", uc
  )
  out <- system2(vsearch, args = args, stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status")
  if (!is.null(status) && !identical(as.integer(status), 0L)) {
    stop("vsearch failed: ", paste(out, collapse = "\n"))
  }
  if (!file.exists(uc)) {
    stop("vsearch did not produce UC output.")
  }

  query_to_centroid <- character(n)
  names(query_to_centroid) <- sprintf("s%d", seq_len(n) - 1L)
  uc_lines <- readLines(uc, warn = FALSE)
  for (line in uc_lines) {
    if (!nzchar(line)) next
    parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(parts) < 9L) next
    rec <- parts[[1]]
    if (rec == "S") {
      query_to_centroid[parts[[9]]] <- parts[[9]]
    } else if (rec == "H" && length(parts) >= 10L) {
      query_to_centroid[parts[[9]]] <- parts[[10]]
    }
  }

  centroids <- unique(query_to_centroid[nzchar(query_to_centroid)])
  if (length(centroids) < 1L) {
    return(rep(NA_integer_, n))
  }

  # Rank centroids by cluster size; optionally keep top max_clusters.
  sizes <- vapply(centroids, function(cen) sum(query_to_centroid == cen, na.rm = TRUE), integer(1))
  ord <- order(sizes, decreasing = TRUE)
  if (is.infinite(max_clusters)) {
    keep <- centroids[ord]
  } else {
    keep <- centroids[ord][seq_len(min(length(centroids), as.integer(max_clusters)))]
  }
  cen_to_id <- stats::setNames(seq_along(keep), keep)

  cluster_ids <- rep(NA_integer_, n)
  for (i in seq_len(n)) {
    q <- sprintf("s%d", i - 1L)
    cen <- query_to_centroid[[q]]
    if (is.na(cen) || !nzchar(cen) || !cen %in% names(cen_to_id)) next
    cluster_ids[[i]] <- cen_to_id[[cen]]
  }
  cluster_ids
}

#' @keywords internal
.consensus_for_members <- function(members, method) {
  do_cons <- method %in% c("both", "consensus")
  if (!do_cons) return(list(seq = members[[1]], method = "greedy"))
  cons <- .consensus_via_abpoa(members)
  list(seq = cons, method = "abpoa")
}

#' @keywords internal
.exact_modal_purity <- function(members) {
  members <- as.character(members)
  members <- members[nzchar(members) & !is.na(members)]
  if (length(members) < 1L) return(NA_real_)
  as.numeric(max(table(members))) / length(members)
}

#' @keywords internal
.filter_clusters_by_purity <- function(clusters, min_cluster_purity) {
  if (is.null(min_cluster_purity) || nrow(clusters) < 1L) {
    return(clusters)
  }
  if (!"cluster_purity" %in% names(clusters)) {
    return(clusters)
  }
  keep <- !is.na(clusters$cluster_purity) &
    clusters$cluster_purity >= as.numeric(min_cluster_purity)
  out <- clusters[keep, , drop = FALSE]
  if (nrow(out) < 1L) {
    return(out)
  }
  out$cluster_id <- seq_len(nrow(out))
  n_gene <- sum(out$n_reads)
  out$n_reads_gene <- as.integer(n_gene)
  out$fraction <- out$n_reads / n_gene
  rownames(out) <- NULL
  out
}

#' @keywords internal
.clusters_from_ids <- function(seqs, cluster_ids, gene_id, index_pair_id,
                                 method, min_cluster_reads,
                                 min_cluster_purity = NULL,
                                 read_id = NULL) {
  n_in <- length(seqs)
  status <- rep(NA_character_, n_in)
  status[is.na(cluster_ids)] <- "max_clusters_overflow"

  empty_members <- data.frame(
    index_pair_id = character(),
    gene_id = character(),
    cluster_id = integer(),
    read_id = character(),
    seq = character(),
    stringsAsFactors = FALSE
  )
  valid <- !is.na(cluster_ids)
  if (!any(valid)) {
    return(list(
      clusters = data.frame(),
      members = empty_members,
      unassigned = .membership_to_unassigned(
        index_pair_id, gene_id, read_id, status
      )
    ))
  }

  seqs_v <- as.character(seqs[valid])
  cluster_ids_v <- cluster_ids[valid]
  read_ids_v <- if (!is.null(read_id)) {
    as.character(read_id[valid])
  } else {
    paste0("r", which(valid))
  }
  valid_idx <- which(valid)

  tab <- sort(table(cluster_ids_v), decreasing = TRUE)
  below <- names(tab)[tab < as.integer(min_cluster_reads)]
  if (length(below) > 0L) {
    for (cid in below) {
      hit <- valid_idx[cluster_ids_v == as.integer(cid)]
      status[hit] <- "below_min_cluster_reads"
    }
  }
  tab <- tab[tab >= as.integer(min_cluster_reads)]
  if (length(tab) < 1L) {
    return(list(
      clusters = data.frame(),
      members = empty_members,
      unassigned = .membership_to_unassigned(
        index_pair_id, gene_id, read_id, status
      )
    ))
  }

  n_gene_reads <- sum(tab)
  out <- vector("list", length(tab))
  member_rows <- vector("list", length(tab))
  retained_orig_cid <- integer(0)

  for (j in seq_along(tab)) {
    cid <- as.integer(names(tab)[[j]])
    idx <- cluster_ids_v == cid
    members <- seqs_v[idx]
    member_ids <- read_ids_v[idx]
    n_reads <- as.integer(tab[[j]])
    purity <- .exact_modal_purity(members)
    cons <- .consensus_for_members(members = members, method = method)
    out[[j]] <- data.frame(
      index_pair_id = index_pair_id,
      gene_id = gene_id,
      cluster_id = j,
      seq = cons$seq,
      n_reads = n_reads,
      n_reads_gene = as.integer(n_gene_reads),
      fraction = n_reads / n_gene_reads,
      fraction_bucket = n_reads / n_in,
      method = cons$method,
      cluster_purity = purity,
      stringsAsFactors = FALSE,
      .orig_cid = cid
    )
    member_rows[[j]] <- data.frame(
      index_pair_id = index_pair_id,
      gene_id = gene_id,
      cluster_id = j,
      read_id = member_ids,
      seq = members,
      .orig_cid = cid,
      stringsAsFactors = FALSE
    )
    retained_orig_cid <- c(retained_orig_cid, cid)
  }
  clusters <- do.call(rbind, out)
  members_df <- do.call(rbind, member_rows)
  clusters <- .filter_clusters_by_purity(clusters, min_cluster_purity)

  pass_orig <- if (nrow(clusters) > 0L && ".orig_cid" %in% names(clusters)) {
    as.integer(clusters$.orig_cid)
  } else {
    integer(0)
  }
  dropped_cid <- setdiff(retained_orig_cid, pass_orig)
  if (length(dropped_cid) > 0L) {
    for (cid in dropped_cid) {
      hit <- valid_idx[cluster_ids_v == cid]
      status[hit] <- "low_cluster_purity"
    }
  }

  if (nrow(clusters) > 0L) {
    orig_to_new <- stats::setNames(
      as.integer(clusters$cluster_id),
      as.character(clusters$.orig_cid)
    )
    members_df <- members_df[
      as.character(members_df$.orig_cid) %in% names(orig_to_new),
      ,
      drop = FALSE
    ]
    members_df$cluster_id <- as.integer(
      orig_to_new[as.character(members_df$.orig_cid)]
    )
    members_df$.orig_cid <- NULL
    clusters$.orig_cid <- NULL
    rownames(members_df) <- NULL
  } else {
    members_df <- empty_members
  }

  list(
    clusters = clusters,
    members = members_df,
    unassigned = .membership_to_unassigned(
      index_pair_id, gene_id, read_id, status
    )
  )
}

#' @keywords internal
.membership_to_unassigned <- function(index_pair_id, gene_id, read_id, status) {
  if (is.null(read_id) || length(read_id) < 1L) {
    return(.empty_unassigned_to_cluster())
  }
  keep <- !is.na(status)
  if (!any(keep)) {
    return(.empty_unassigned_to_cluster())
  }
  data.frame(
    index_pair_id = as.character(index_pair_id)[[1]],
    gene_id = as.character(gene_id)[[1]],
    read_id = as.character(read_id)[keep],
    reason = as.character(status)[keep],
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
.empty_unassigned_to_cluster <- function() {
  data.frame(
    index_pair_id = character(),
    gene_id = character(),
    read_id = character(),
    reason = character(),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
.unassigned_rows <- function(index_pair_id, gene_id, read_id, reason) {
  if (length(read_id) < 1L) {
    return(.empty_unassigned_to_cluster())
  }
  data.frame(
    index_pair_id = as.character(index_pair_id)[[1]],
    gene_id = as.character(gene_id)[[1]],
    read_id = as.character(read_id),
    reason = as.character(reason)[[1]],
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
.cluster_and_consensus_reads <- function(seqs, gene_id, index_pair_id, method,
                                         min_cluster_reads, max_clusters,
                                         min_cluster_identity,
                                         cluster_backend = "vsearch",
                                         min_cluster_purity = NULL,
                                         read_id = NULL) {
  empty <- list(
    clusters = data.frame(),
    members = data.frame(
      index_pair_id = character(), gene_id = character(),
      cluster_id = integer(), read_id = character(), seq = character(),
      stringsAsFactors = FALSE
    ),
    unassigned = .empty_unassigned_to_cluster()
  )
  if (length(seqs) < 1L) {
    return(empty)
  }
  .cluster_and_consensus_core(
    seqs = seqs,
    gene_id = gene_id,
    index_pair_id = index_pair_id,
    method = method,
    min_cluster_reads = min_cluster_reads,
    max_clusters = max_clusters,
    min_cluster_identity = min_cluster_identity,
    cluster_backend = cluster_backend,
    min_cluster_purity = min_cluster_purity,
    read_id = read_id
  )
}

#' @keywords internal
.cluster_and_consensus_core <- function(seqs, gene_id, index_pair_id, method,
                                        min_cluster_reads, max_clusters,
                                        min_cluster_identity,
                                        cluster_backend = "vsearch",
                                        min_cluster_purity = NULL,
                                        read_id = NULL) {
  if (!identical(cluster_backend, "vsearch")) {
    stop("Only cluster_backend = \"vsearch\" is supported.")
  }
  cluster_ids <- .cluster_ids_vsearch(
    seqs = as.character(seqs),
    min_cluster_identity = min_cluster_identity,
    max_clusters = max_clusters
  )
  .clusters_from_ids(
    seqs = seqs,
    cluster_ids = cluster_ids,
    gene_id = gene_id,
    index_pair_id = index_pair_id,
    method = method,
    min_cluster_reads = min_cluster_reads,
    min_cluster_purity = min_cluster_purity,
    read_id = read_id
  )
}

#' @keywords internal
.overlap_graph_merge_edit <- function(min_cluster_identity, overlap_min_identity,
                                      seqs) {
  lens <- nchar(as.character(seqs))
  lens <- lens[lens > 0L]
  med_len <- if (length(lens)) stats::median(lens) else 200
  id <- min(as.numeric(min_cluster_identity), as.numeric(overlap_min_identity))
  as.integer(max(floor((1 - id) * med_len), 0L))
}

#' @keywords internal
.assemble_overlap_graph_reads <- function(seqs, gene_id, index_pair_id,
                                          min_cluster_reads, max_clusters,
                                          min_cluster_identity,
                                          overlap_min_identity = 0.90,
                                          min_cluster_purity = NULL,
                                          read_id = NULL) {
  seqs <- as.character(seqs)
  empty_members <- data.frame(
    index_pair_id = character(), gene_id = character(),
    cluster_id = integer(), read_id = character(), seq = character(),
    stringsAsFactors = FALSE
  )
  empty <- list(
    clusters = data.frame(),
    members = empty_members,
    unassigned = .empty_unassigned_to_cluster()
  )
  if (length(seqs) < 1L) {
    return(empty)
  }

  merge_edit <- .overlap_graph_merge_edit(
    min_cluster_identity = min_cluster_identity,
    overlap_min_identity = overlap_min_identity,
    seqs = seqs
  )

  tab <- sort(table(seqs), decreasing = TRUE)
  uniq <- names(tab)
  exact_counts <- as.integer(tab)

  kept_idx <- integer(0)
  for (i in seq_along(uniq)) {
    if (exact_counts[[i]] < min_cluster_reads) {
      next
    }
    dominated <- FALSE
    if (length(kept_idx) > 0L) {
      for (j in kept_idx) {
        d <- edlib_edit_distance_cpp(
          uniq[[i]],
          uniq[[j]],
          max_edit = merge_edit
        )
        if (!is.na(d) && d <= merge_edit) {
          dominated <- TRUE
          break
        }
      }
    }
    if (!dominated) {
      kept_idx <- c(kept_idx, i)
    }
    if (!is.infinite(max_clusters) &&
        length(kept_idx) >= as.integer(max_clusters)) {
      break
    }
  }
  if (length(kept_idx) < 1L) {
    return(list(
      clusters = data.frame(),
      members = empty_members,
      unassigned = .unassigned_rows(
        index_pair_id, gene_id, read_id, "below_min_cluster_reads"
      )
    ))
  }

  paths <- uniq[kept_idx]
  path_counts <- integer(length(paths))
  path_exact <- integer(length(paths))
  best_path <- rep(NA_integer_, length(seqs))
  read_ids <- if (!is.null(read_id)) {
    as.character(read_id)
  } else {
    paste0("r", seq_along(seqs))
  }
  for (si in seq_along(seqs)) {
    s <- seqs[[si]]
    best <- NA_integer_
    best_d <- merge_edit + 1L
    for (p in seq_along(paths)) {
      d <- edlib_edit_distance_cpp(s, paths[[p]], max_edit = merge_edit)
      if (!is.na(d) && d <= merge_edit && d < best_d) {
        best_d <- d
        best <- p
      }
    }
    if (!is.na(best)) {
      path_counts[[best]] <- path_counts[[best]] + 1L
      best_path[[si]] <- best
      if (identical(s, paths[[best]])) {
        path_exact[[best]] <- path_exact[[best]] + 1L
      }
    }
  }

  valid <- which(path_counts >= as.integer(min_cluster_reads))
  if (length(valid) < 1L) {
    return(list(
      clusters = data.frame(),
      members = empty_members,
      unassigned = .unassigned_rows(
        index_pair_id, gene_id, read_id, "below_min_cluster_reads"
      )
    ))
  }

  n_gene_reads <- sum(path_counts[valid])
  n_in <- length(seqs)
  out <- lapply(seq_along(valid), function(k) {
    j <- valid[[k]]
    n_reads <- as.integer(path_counts[[j]])
    purity <- if (n_reads > 0L) {
      as.numeric(path_exact[[j]]) / n_reads
    } else {
      NA_real_
    }
    data.frame(
      index_pair_id = as.character(index_pair_id)[[1]],
      gene_id = as.character(gene_id)[[1]],
      cluster_id = k,
      seq = paths[[j]],
      n_reads = n_reads,
      n_reads_gene = as.integer(n_gene_reads),
      fraction = path_counts[[j]] / n_gene_reads,
      fraction_bucket = path_counts[[j]] / n_in,
      method = "overlap_graph",
      cluster_purity = purity,
      stringsAsFactors = FALSE,
      .path_idx = j
    )
  })
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res <- .filter_clusters_by_purity(res, min_cluster_purity)

  kept_paths <- if (".path_idx" %in% names(res) && nrow(res) > 0L) {
    as.integer(res$.path_idx)
  } else {
    integer(0)
  }
  path_to_cid <- if (nrow(res) > 0L) {
    stats::setNames(as.integer(res$cluster_id), as.character(res$.path_idx))
  } else {
    integer(0)
  }
  if (".path_idx" %in% names(res)) {
    res$.path_idx <- NULL
  }

  status <- rep(NA_character_, length(seqs))
  mem_list <- list()
  for (si in seq_along(seqs)) {
    bp <- best_path[[si]]
    if (is.na(bp)) {
      status[[si]] <- "overlap_unmatched"
    } else if (!(bp %in% valid)) {
      status[[si]] <- "below_min_cluster_reads"
    } else if (!(bp %in% kept_paths)) {
      status[[si]] <- "low_cluster_purity"
    } else {
      cid <- as.integer(path_to_cid[[as.character(bp)]])
      mem_list[[length(mem_list) + 1L]] <- data.frame(
        index_pair_id = as.character(index_pair_id)[[1]],
        gene_id = as.character(gene_id)[[1]],
        cluster_id = cid,
        read_id = read_ids[[si]],
        seq = seqs[[si]],
        stringsAsFactors = FALSE
      )
    }
  }
  members_df <- if (length(mem_list)) {
    do.call(rbind, mem_list)
  } else {
    empty_members
  }

  list(
    clusters = res,
    members = members_df,
    unassigned = .membership_to_unassigned(
      index_pair_id, gene_id, read_id, status
    )
  )
}
