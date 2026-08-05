################################################################################
# Phase H: reassessment of assemblies (recall diagnostic; Path A)
################################################################################

#' Reassess assembled amplicon consensus and unassigned reads
#'
#' Diagnostic step after [doAssembleAmplicons()]. Computes (Q1) pairwise
#' consensus similarity (and optional near-duplicate groups for review) and
#' (Q2) whether U1 unassigned reads can be reassigned to an existing
#' consensus. Does **not** auto-merge or rewrite sequences: the user chooses
#' which consensus to keep from the pairwise report.
#'
#' Does **not** modify `amplicon/`. Results are written under `out_dir`
#' (typically `amplicon_reassess/`), with one subdirectory per `backend`
#' so edlib / VSEARCH / BLASTN / MMseqs2 can be compared.
#'
#' @param amplicon_out Result of [doAssembleAmplicons()] or path to `amplicon/`.
#' @param gene_assign Gene assignments with `oriented_seq` (same contract as
#'   assemble). Used to recover U1 sequences when
#'   `unassigned_to_cluster.tsv` is present, and as fallback when that file is
#'   missing.
#' @param out_dir Output directory (e.g. `.../amplicon_reassess`).
#' @param backend Character vector of backends to run. Allowed:
#'   `"edlib"`, `"vsearch"`, `"blastn"`, `"mmseqs"`. Default `"edlib"` provides
#'   an independent exact edit-distance evaluation of clustering results (a
#'   different perspective from vsearch %identity used in clustering). Multiple
#'   values produce parallel result trees and a `summary_compare.tsv`.
#' @param consensus_merge_max_edit Max edit distance (bp) used only as a
#'   **highlight threshold** for Q1 pairwise / near-duplicate grouping
#'   (`pass_threshold`). No sequences are merged. Identity-based backends use
#'   the derived `min_identity` unless that argument is set explicitly.
#' @param read_assign_max_edit Max edit distance (bp) for U1 → consensus
#'   reassignment diagnostics (Q2). Same identity derivation as above.
#' @param min_identity Optional identity threshold in \eqn{[0,1]} for
#'   `vsearch` / `blastn` / `mmseqs`. If `NULL`, derived from edit distance and
#'   median consensus length via [`.vsearch_id_from_edit()`].
#' @param layers Unassigned layers to include. Currently `"U1"` only.
#' @param primer_list Optional primers for insert trim when reconstructing U1
#'   sequences from `gene_assign`.
#' @param max_primer_edit Primer trim edit distance.
#' @param skip_missing_backend If `TRUE` (default), backends whose binary is
#'   missing on `PATH` are skipped with a warning. If `FALSE`, stop.
#' @param overwrite If `FALSE`, stop when `out_dir` is non-empty.
#' @param n_core Reserved / passed to primer trim.
#'
#' @return A list with `out_dir`, `summary`, `summary_compare`, `backends_run`,
#'   and `samples`.
#' @export
doReassessAssemblies <- function(
  amplicon_out,
  gene_assign = NULL,
  out_dir,
  backend = "edlib",
  consensus_merge_max_edit = 12L,
  read_assign_max_edit = 12L,
  min_identity = NULL,
  layers = c("U1"),
  primer_list = NULL,
  max_primer_edit = 5L,
  skip_missing_backend = TRUE,
  overwrite = FALSE,
  n_core = 1L
) {
  allowed <- c("edlib", "vsearch", "blastn", "mmseqs")
  backend <- unique(as.character(backend))
  if (length(backend) < 1L) {
    stop("backend must contain at least one of: ", paste(allowed, collapse = ", "))
  }
  bad <- setdiff(backend, allowed)
  if (length(bad) > 0L) {
    stop("Unknown backend(s): ", paste(bad, collapse = ", "))
  }

  layers <- unique(as.character(layers))
  if (!all(layers %in% c("U0", "U1"))) {
    stop("layers must be a subset of c(\"U0\", \"U1\").")
  }
  if (!"U1" %in% layers) {
    stop("layers must include \"U1\".")
  }
  if ("U0" %in% layers) {
    warning("U0 is not implemented; ignoring.")
  }

  amp_dir <- if (is.list(amplicon_out)) amplicon_out$out_dir else amplicon_out
  if (is.null(amp_dir) || !dir.exists(amp_dir)) {
    stop("amplicon_out directory does not exist: ", amp_dir)
  }
  if (!overwrite && dir.exists(out_dir) && length(list.files(out_dir)) > 0L) {
    stop("out_dir already exists and is not empty. Set overwrite = TRUE.")
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  backends_run <- character()
  for (b in backend) {
    ok <- .reassess_backend_available(b)
    if (isTRUE(ok$ok)) {
      backends_run <- c(backends_run, b)
    } else if (isTRUE(skip_missing_backend)) {
      warning("Skipping backend \"", b, "\": ", ok$reason)
    } else {
      stop("Backend \"", b, "\" unavailable: ", ok$reason)
    }
  }
  if (length(backends_run) < 1L) {
    stop("No usable backends after availability checks.")
  }

  ga <- if (!is.null(gene_assign)) {
    .as_gene_assign_df(gene_assign)
  } else if (is.list(amplicon_out) && !is.null(amplicon_out$gene_assignments)) {
    .as_gene_assign_df(amplicon_out$gene_assignments)
  } else {
    NULL
  }
  primers <- if (!is.null(primer_list)) .parse_primer_pairs(primer_list) else NULL

  pair_dirs <- list.dirs(amp_dir, recursive = FALSE, full.names = TRUE)
  pair_dirs <- pair_dirs[
    file.exists(file.path(pair_dirs, "cluster_counts.tsv")) |
      file.exists(file.path(pair_dirs, "consensus.fasta"))
  ]
  if (length(pair_dirs) < 1L) {
    stop("No per-sample assemble outputs found under: ", amp_dir)
  }

  summary_rows <- list()
  for (b in backends_run) {
    b_root <- file.path(out_dir, b)
    dir.create(b_root, recursive = TRUE, showWarnings = FALSE)
    for (pair_dir_path in pair_dirs) {
      pair_id <- basename(pair_dir_path)
      summary_rows[[length(summary_rows) + 1L]] <- .reassess_one_sample(
        index_pair_id = pair_id,
        amplicon_pair_dir = pair_dir_path,
        out_pair_dir = file.path(b_root, pair_id),
        gene_assign = ga,
        primers = primers,
        max_primer_edit = as.integer(max_primer_edit),
        consensus_merge_max_edit = as.integer(consensus_merge_max_edit),
        read_assign_max_edit = as.integer(read_assign_max_edit),
        min_identity = min_identity,
        backend = b,
        n_core = as.integer(n_core)
      )
    }
  }

  summary_df <- do.call("rbind", summary_rows)
  rownames(summary_df) <- NULL
  utils::write.table(
    summary_df,
    file = file.path(out_dir, "summary_by_sample.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  compare_df <- .reassess_summary_compare(summary_df)
  utils::write.table(
    compare_df,
    file = file.path(out_dir, "summary_compare.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  versions <- .tool_versions(backends_run)
  pkg_ver <- tryCatch(
    as.character(utils::packageVersion("miaoseq")),
    error = function(e) NA_character_
  )
  write(
    paste0(
      "r_version=", R.version.string, "\n",
      "miaoseq_version=", pkg_ver, "\n",
      "backends_requested=", paste(backend, collapse = ","), "\n",
      "backends_run=", paste(backends_run, collapse = ","), "\n",
      "consensus_merge_max_edit=", consensus_merge_max_edit, "\n",
      "read_assign_max_edit=", read_assign_max_edit, "\n",
      "min_identity=", if (is.null(min_identity)) "auto" else min_identity, "\n",
      "n_samples=", length(pair_dirs), "\n",
      paste(paste0("version_", names(versions), "=", versions), collapse = "\n"),
      "\n"
    ),
    file = file.path(out_dir, "run_stats.txt")
  )

  list(
    samples = unique(summary_df$index_pair_id),
    out_dir = out_dir,
    summary = summary_df,
    summary_compare = compare_df,
    backends_run = backends_run
  )
}

#' @keywords internal
.reassess_backend_available <- function(backend) {
  if (identical(backend, "edlib")) {
    return(list(ok = TRUE, reason = ""))
  }
  bins <- switch(
    backend,
    vsearch = "vsearch",
    blastn = c("blastn", "makeblastdb"),
    mmseqs = "mmseqs",
    character()
  )
  missing <- bins[!nzchar(Sys.which(bins))]
  if (length(missing) > 0L) {
    return(list(
      ok = FALSE,
      reason = paste0("missing binary: ", paste(missing, collapse = ", "))
    ))
  }
  list(ok = TRUE, reason = "")
}

#' @keywords internal
.reassess_summary_compare <- function(summary_df) {
  if (nrow(summary_df) < 1L) {
    return(data.frame(
      index_pair_id = character(),
      stringsAsFactors = FALSE
    ))
  }
  backends <- unique(summary_df$backend)
  base <- unique(summary_df[, c("index_pair_id", "n_consensus", "n_unassigned_U1"), drop = FALSE])
  out <- base
  for (b in backends) {
    sub <- summary_df[summary_df$backend == b, , drop = FALSE]
    tmp <- data.frame(
      index_pair_id = sub$index_pair_id,
      stringsAsFactors = FALSE
    )
    tmp[[paste0(b, "_n_near_duplicate_groups")]] <- sub$n_near_duplicate_groups
    tmp[[paste0(b, "_reassignable_frac")]] <- sub$reassignable_frac
    tmp[[paste0(b, "_n_reassignable_U1")]] <- sub$n_reassignable_U1
    tmp[[paste0(b, "_min_identity")]] <- sub$min_identity
    out <- merge(out, tmp, by = "index_pair_id", all.x = TRUE)
  }
  out[order(out$index_pair_id), , drop = FALSE]
}

#' @keywords internal
.reassess_one_sample <- function(index_pair_id,
                                 amplicon_pair_dir,
                                 out_pair_dir,
                                 gene_assign,
                                 primers,
                                 max_primer_edit,
                                 consensus_merge_max_edit,
                                 read_assign_max_edit,
                                 min_identity,
                                 backend,
                                 n_core) {
  dir.create(out_pair_dir, recursive = TRUE, showWarnings = FALSE)

  consensus <- .load_sample_consensus_table(amplicon_pair_dir, index_pair_id)
  n_consensus <- nrow(consensus)

  id_merge <- if (!is.null(min_identity)) {
    as.numeric(min_identity)
  } else {
    .vsearch_id_from_edit(
      as.character(consensus$seq),
      as.integer(consensus_merge_max_edit)
    )
  }
  id_assign <- if (!is.null(min_identity)) {
    as.numeric(min_identity)
  } else {
    .vsearch_id_from_edit(
      as.character(consensus$seq),
      as.integer(read_assign_max_edit)
    )
  }

  pairwise <- .consensus_pairwise(
    consensus = consensus,
    backend = backend,
    max_edit = consensus_merge_max_edit,
    min_identity = id_merge
  )
  utils::write.table(
    pairwise,
    file = file.path(out_pair_dir, "consensus_pairwise.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  # Soft grouping of near-duplicates for review only (no merged sequences).
  near_dup <- .consensus_near_duplicate_groups(
    consensus = consensus,
    pairwise = pairwise
  )
  utils::write.table(
    near_dup,
    file = file.path(out_pair_dir, "consensus_near_duplicate_groups.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  u1 <- .load_or_infer_u1_reads(
    index_pair_id = index_pair_id,
    amplicon_pair_dir = amplicon_pair_dir,
    gene_assign = gene_assign,
    consensus = consensus,
    primers = primers,
    max_primer_edit = max_primer_edit,
    read_assign_max_edit = read_assign_max_edit,
    n_core = n_core
  )
  assign_tbl <- .assign_u1_to_consensus(
    u1 = u1,
    consensus = consensus,
    backend = backend,
    max_edit = read_assign_max_edit,
    min_identity = id_assign
  )
  utils::write.table(
    assign_tbl,
    file = file.path(out_pair_dir, "unassigned_to_consensus.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  n_u1 <- nrow(assign_tbl)
  n_reassignable <- if (n_u1 > 0L) sum(assign_tbl$pass_threshold, na.rm = TRUE) else 0L
  reassignable_frac <- if (n_u1 > 0L) {
    n_reassignable / n_u1
  } else {
    NA_real_
  }

  data.frame(
    index_pair_id = index_pair_id,
    backend = backend,
    n_consensus = as.integer(n_consensus),
    n_near_duplicate_groups = as.integer(
      if (nrow(near_dup) < 1L) 0L else length(unique(near_dup$group_id))
    ),
    n_unassigned_U1 = as.integer(n_u1),
    n_reassignable_U1 = as.integer(n_reassignable),
    reassignable_frac = as.numeric(reassignable_frac),
    consensus_merge_max_edit = as.integer(consensus_merge_max_edit),
    read_assign_max_edit = as.integer(read_assign_max_edit),
    min_identity = as.numeric(id_merge),
    stringsAsFactors = FALSE
  )
}

################################################################################
# Loaders / shared geometry (unchanged logic)
################################################################################

#' @keywords internal
.load_sample_consensus_table <- function(sample_dir, index_pair_id) {
  counts_path <- file.path(sample_dir, "cluster_counts.tsv")
  fa_path <- file.path(sample_dir, "consensus.fasta")
  if (file.exists(counts_path)) {
    counts <- .ensure_index_pair_col(utils::read.delim(counts_path, stringsAsFactors = FALSE))
  } else {
    counts <- data.frame(
      index_pair_id = character(),
      gene_id = character(),
      cluster_id = integer(),
      n_reads = integer(),
      stringsAsFactors = FALSE
    )
  }
  if (!file.exists(fa_path)) {
    counts$seq <- character(nrow(counts))
    return(counts)
  }
  dna <- Biostrings::readDNAStringSet(fa_path)
  hdr <- names(dna)
  gene_id <- sub("\\|.*", "", hdr)
  cluster_id <- as.integer(sub("^.*cluster_([0-9]+).*$", "\\1", hdr))
  seq_map <- data.frame(
    gene_id = gene_id,
    cluster_id = cluster_id,
    seq = as.character(dna),
    stringsAsFactors = FALSE
  )
  if (nrow(counts) < 1L) {
    seq_map$index_pair_id <- index_pair_id
    seq_map$n_reads <- NA_integer_
    return(seq_map)
  }
  if (!"seq" %in% names(counts)) {
    counts <- merge(counts, seq_map, by = c("gene_id", "cluster_id"), all.x = TRUE)
  }
  if (!"index_pair_id" %in% names(counts) ||
      any(is.na(counts$index_pair_id) | !nzchar(as.character(counts$index_pair_id)))) {
    counts$index_pair_id <- index_pair_id
  }
  counts
}

#' @keywords internal
.empty_pairwise <- function() {
  data.frame(
    gene_id = character(),
    cluster_i = integer(),
    cluster_j = integer(),
    tool = character(),
    metric = character(),
    value = numeric(),
    pident = numeric(),
    edit_distance = integer(),
    pass_threshold = logical(),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
.empty_assign <- function() {
  data.frame(
    index_pair_id = character(),
    gene_id = character(),
    read_id = character(),
    reason = character(),
    layer = character(),
    best_cluster_id = integer(),
    tool = character(),
    metric = character(),
    value = numeric(),
    pident = numeric(),
    edit_distance = integer(),
    pass_threshold = logical(),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
.write_named_fasta <- function(ids, seqs, path) {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  for (i in seq_along(ids)) {
    writeLines(paste0(">", ids[[i]]), con)
    writeLines(as.character(seqs[[i]]), con)
  }
  invisible(path)
}

#' @keywords internal
.consensus_id_label <- function(gene_id, cluster_id) {
  paste0(gene_id, "__c", as.integer(cluster_id))
}

#' @keywords internal
.parse_consensus_id_label <- function(label) {
  gene_id <- sub("__c[0-9]+$", "", label)
  cluster_id <- as.integer(sub("^.*__c", "", label))
  list(gene_id = gene_id, cluster_id = cluster_id)
}

################################################################################
# Backend dispatch
################################################################################

#' @keywords internal
.consensus_pairwise <- function(consensus, backend, max_edit, min_identity) {
  if (nrow(consensus) < 2L) {
    return(.empty_pairwise())
  }
  switch(
    backend,
    edlib = .consensus_pairwise_edlib(consensus, max_edit),
    vsearch = .consensus_pairwise_search_tool(
      consensus, max_edit, min_identity, tool = "vsearch"
    ),
    blastn = .consensus_pairwise_search_tool(
      consensus, max_edit, min_identity, tool = "blastn"
    ),
    mmseqs = .consensus_pairwise_search_tool(
      consensus, max_edit, min_identity, tool = "mmseqs"
    ),
    stop("Unknown backend: ", backend)
  )
}

#' @keywords internal
.assign_u1_to_consensus <- function(u1, consensus, backend, max_edit,
                                    min_identity) {
  if (nrow(u1) < 1L) {
    return(.empty_assign())
  }
  switch(
    backend,
    edlib = .assign_u1_to_consensus_edlib(u1, consensus, max_edit),
    vsearch = .assign_u1_search_tool(
      u1, consensus, max_edit, min_identity, tool = "vsearch"
    ),
    blastn = .assign_u1_search_tool(
      u1, consensus, max_edit, min_identity, tool = "blastn"
    ),
    mmseqs = .assign_u1_search_tool(
      u1, consensus, max_edit, min_identity, tool = "mmseqs"
    ),
    stop("Unknown backend: ", backend)
  )
}

################################################################################
# edlib
################################################################################

#' @keywords internal
.consensus_pairwise_edlib <- function(consensus, max_edit) {
  if (nrow(consensus) < 1L) return(.empty_pairwise())
  rows <- list()
  for (gid in unique(consensus$gene_id)) {
    sub <- consensus[consensus$gene_id == gid, , drop = FALSE]
    if (nrow(sub) < 2L) next
    for (a in seq_len(nrow(sub) - 1L)) {
      for (b in seq.int(a + 1L, nrow(sub))) {
        sa <- as.character(sub$seq[[a]])
        sb <- as.character(sub$seq[[b]])
        d <- edlib_edit_distance_cpp(sa, sb, max_edit = as.integer(max_edit))
        pass <- !is.na(d) && d <= as.integer(max_edit)
        aln_len <- max(nchar(sa), nchar(sb))
        pident <- if (is.na(d) || aln_len < 1L) {
          NA_real_
        } else {
          max(0, 1 - as.numeric(d) / aln_len) * 100
        }
        rows[[length(rows) + 1L]] <- data.frame(
          gene_id = gid,
          cluster_i = as.integer(sub$cluster_id[[a]]),
          cluster_j = as.integer(sub$cluster_id[[b]]),
          tool = "edlib",
          metric = "edit_distance",
          value = if (is.na(d)) NA_real_ else as.numeric(d),
          pident = pident,
          edit_distance = if (is.na(d)) NA_integer_ else as.integer(d),
          pass_threshold = pass,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (!length(rows)) return(.empty_pairwise())
  do.call(rbind, rows)
}

#' @keywords internal
.assign_u1_to_consensus_edlib <- function(u1, consensus, max_edit) {
  if (nrow(u1) < 1L) return(.empty_assign())
  rows <- vector("list", nrow(u1))
  for (i in seq_len(nrow(u1))) {
    gid <- u1$gene_id[[i]]
    cons <- consensus[consensus$gene_id == gid, , drop = FALSE]
    best_id <- NA_integer_
    best_d <- NA_integer_
    best_pident <- NA_real_
    if (nrow(cons) > 0L && !is.na(u1$seq[[i]]) && nzchar(u1$seq[[i]])) {
      for (j in seq_len(nrow(cons))) {
        d <- edlib_edit_distance_cpp(
          u1$seq[[i]],
          cons$seq[[j]],
          max_edit = as.integer(max_edit)
        )
        if (!is.na(d) && (is.na(best_d) || d < best_d)) {
          best_d <- as.integer(d)
          best_id <- as.integer(cons$cluster_id[[j]])
          aln_len <- max(nchar(u1$seq[[i]]), nchar(cons$seq[[j]]))
          best_pident <- if (aln_len < 1L) NA_real_ else max(0, 1 - d / aln_len) * 100
        }
      }
    }
    pass <- !is.na(best_d) && best_d <= as.integer(max_edit)
    rows[[i]] <- data.frame(
      index_pair_id = u1$index_pair_id[[i]],
      gene_id = gid,
      read_id = u1$read_id[[i]],
      reason = u1$reason[[i]],
      layer = u1$layer[[i]],
      best_cluster_id = best_id,
      tool = "edlib",
      metric = "edit_distance",
      value = if (is.na(best_d)) NA_real_ else as.numeric(best_d),
      pident = best_pident,
      edit_distance = best_d,
      pass_threshold = pass,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

################################################################################
# External search backends (vsearch / blastn / mmseqs)
################################################################################

#' @keywords internal
.consensus_pairwise_search_tool <- function(consensus, max_edit, min_identity,
                                            tool) {
  rows <- list()
  for (gid in unique(consensus$gene_id)) {
    sub <- consensus[consensus$gene_id == gid, , drop = FALSE]
    if (nrow(sub) < 2L) next
    ids <- .consensus_id_label(sub$gene_id, sub$cluster_id)
    hits <- .search_all_pairs(
      ids = ids,
      seqs = as.character(sub$seq),
      tool = tool,
      min_identity = min_identity
    )
    if (nrow(hits) < 1L) next
    for (k in seq_len(nrow(hits))) {
      qi <- .parse_consensus_id_label(hits$query[[k]])
      tj <- .parse_consensus_id_label(hits$target[[k]])
      if (!identical(qi$gene_id, gid) || !identical(tj$gene_id, gid)) next
      if (is.na(qi$cluster_id) || is.na(tj$cluster_id)) next
      if (qi$cluster_id >= tj$cluster_id) next
      pident <- as.numeric(hits$pident[[k]])
      aln_len <- as.integer(hits$alnlen[[k]])
      approx_edit <- if (!is.na(pident) && !is.na(aln_len) && aln_len > 0L) {
        as.integer(round((1 - pident / 100) * aln_len))
      } else {
        NA_integer_
      }
      pass <- !is.na(pident) && (pident / 100) >= as.numeric(min_identity)
      rows[[length(rows) + 1L]] <- data.frame(
        gene_id = gid,
        cluster_i = qi$cluster_id,
        cluster_j = tj$cluster_id,
        tool = tool,
        metric = "pident",
        value = pident,
        pident = pident,
        edit_distance = approx_edit,
        pass_threshold = pass,
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(rows)) return(.empty_pairwise())
  # Keep best (highest pident) for each unordered pair.
  out <- do.call(rbind, rows)
  out <- out[order(out$gene_id, out$cluster_i, out$cluster_j, -out$pident), , drop = FALSE]
  out <- out[!duplicated(paste(out$gene_id, out$cluster_i, out$cluster_j, sep = "\t")), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' @keywords internal
.assign_u1_search_tool <- function(u1, consensus, max_edit, min_identity, tool) {
  if (nrow(u1) < 1L || nrow(consensus) < 1L) return(.empty_assign())

  # Search per gene to avoid cross-gene best hits.
  parts <- list()
  for (gid in unique(u1$gene_id)) {
    u_sub <- u1[u1$gene_id == gid & !is.na(u1$seq) & nzchar(u1$seq), , drop = FALSE]
    c_sub <- consensus[consensus$gene_id == gid, , drop = FALSE]
    if (nrow(u_sub) < 1L || nrow(c_sub) < 1L) {
      for (i in which(u1$gene_id == gid)) {
        parts[[length(parts) + 1L]] <- data.frame(
          index_pair_id = u1$index_pair_id[[i]],
          gene_id = gid,
          read_id = u1$read_id[[i]],
          reason = u1$reason[[i]],
          layer = u1$layer[[i]],
          best_cluster_id = NA_integer_,
          tool = tool,
          metric = "pident",
          value = NA_real_,
          pident = NA_real_,
          edit_distance = NA_integer_,
          pass_threshold = FALSE,
          stringsAsFactors = FALSE
        )
      }
      next
    }
    q_ids <- as.character(u_sub$read_id)
    t_ids <- .consensus_id_label(c_sub$gene_id, c_sub$cluster_id)
    hits <- .search_queries_vs_targets(
      query_ids = q_ids,
      query_seqs = as.character(u_sub$seq),
      target_ids = t_ids,
      target_seqs = as.character(c_sub$seq),
      tool = tool,
      min_identity = min_identity
    )
    best <- .best_hit_per_query(hits)
    for (i in seq_len(nrow(u_sub))) {
      rid <- u_sub$read_id[[i]]
      h <- best[[rid]]
      if (is.null(h)) {
        parts[[length(parts) + 1L]] <- data.frame(
          index_pair_id = u_sub$index_pair_id[[i]],
          gene_id = gid,
          read_id = rid,
          reason = u_sub$reason[[i]],
          layer = u_sub$layer[[i]],
          best_cluster_id = NA_integer_,
          tool = tool,
          metric = "pident",
          value = NA_real_,
          pident = NA_real_,
          edit_distance = NA_integer_,
          pass_threshold = FALSE,
          stringsAsFactors = FALSE
        )
      } else {
        parsed <- .parse_consensus_id_label(h$target)
        pident <- as.numeric(h$pident)
        aln_len <- as.integer(h$alnlen)
        approx_edit <- if (!is.na(pident) && !is.na(aln_len) && aln_len > 0L) {
          as.integer(round((1 - pident / 100) * aln_len))
        } else {
          NA_integer_
        }
        pass <- !is.na(pident) && (pident / 100) >= as.numeric(min_identity)
        parts[[length(parts) + 1L]] <- data.frame(
          index_pair_id = u_sub$index_pair_id[[i]],
          gene_id = gid,
          read_id = rid,
          reason = u_sub$reason[[i]],
          layer = u_sub$layer[[i]],
          best_cluster_id = parsed$cluster_id,
          tool = tool,
          metric = "pident",
          value = pident,
          pident = pident,
          edit_distance = approx_edit,
          pass_threshold = pass,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (!length(parts)) return(.empty_assign())
  do.call(rbind, parts)
}

#' @keywords internal
.best_hit_per_query <- function(hits) {
  out <- list()
  if (nrow(hits) < 1L) return(out)
  hits <- hits[order(hits$query, -hits$pident, hits$evalue), , drop = FALSE]
  for (i in seq_len(nrow(hits))) {
    q <- hits$query[[i]]
    if (is.null(out[[q]])) {
      out[[q]] <- list(
        target = hits$target[[i]],
        pident = hits$pident[[i]],
        alnlen = hits$alnlen[[i]],
        evalue = hits$evalue[[i]]
      )
    }
  }
  out
}

#' @keywords internal
.search_all_pairs <- function(ids, seqs, tool, min_identity) {
  .search_queries_vs_targets(
    query_ids = ids,
    query_seqs = seqs,
    target_ids = ids,
    target_seqs = seqs,
    tool = tool,
    min_identity = min_identity
  )
}

#' @keywords internal
.search_queries_vs_targets <- function(query_ids, query_seqs, target_ids,
                                       target_seqs, tool, min_identity) {
  empty <- data.frame(
    query = character(),
    target = character(),
    pident = numeric(),
    alnlen = integer(),
    evalue = numeric(),
    stringsAsFactors = FALSE
  )
  if (length(query_ids) < 1L || length(target_ids) < 1L) return(empty)

  tmp <- tempfile(paste0("reassess_", tool, "_"))
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  qfa <- file.path(tmp, "query.fa")
  tfa <- file.path(tmp, "target.fa")
  .write_named_fasta(query_ids, query_seqs, qfa)
  .write_named_fasta(target_ids, target_seqs, tfa)

  switch(
    tool,
    vsearch = .run_vsearch_usearch_global(qfa, tfa, tmp, min_identity),
    blastn = .run_blastn_search(qfa, tfa, tmp),
    mmseqs = .run_mmseqs_easy_search(qfa, tfa, tmp),
    stop("Unknown search tool: ", tool)
  )
}

#' Run an external binary without forcing a shell rewrite of argv.
#'
#' `system2(..., stdout = TRUE, stderr = TRUE)` joins args with spaces via
#' `/bin/sh`, which breaks values that contain spaces (e.g. BLAST `-outfmt`).
#' Capture stderr to a temp file instead so argv stays intact.
#'
#' @keywords internal
.system2_safe <- function(command, args) {
  err_file <- tempfile("miaoseq_sys2_err_")
  on.exit(unlink(err_file), add = TRUE)
  out <- system2(command, args = args, stdout = TRUE, stderr = err_file)
  status <- attr(out, "status")
  err <- if (file.exists(err_file) && file.info(err_file)$size > 0) {
    paste(readLines(err_file, warn = FALSE), collapse = "\n")
  } else {
    ""
  }
  list(
    stdout = out,
    stderr = err,
    status = if (is.null(status)) 0L else as.integer(status)
  )
}

#' @keywords internal
.run_vsearch_usearch_global <- function(query_fa, target_fa, tmp, min_identity) {
  vsearch <- Sys.which("vsearch")
  out <- file.path(tmp, "vsearch.tsv")
  id <- max(0.5, min(0.99, as.numeric(min_identity)))
  args <- c(
    "--usearch_global", query_fa,
    "--db", target_fa,
    "--id", format(id, digits = 6, trim = TRUE),
    "--iddef", "2",
    "--strand", "plus",
    "--maxaccepts", "0",
    "--maxrejects", "0",
    "--userout", out,
    "--userfields", "query+target+id+alnlen+mism+gaps"
  )
  res <- .system2_safe(vsearch, args)
  if (!identical(res$status, 0L)) {
    stop(
      "vsearch --usearch_global failed: ",
      paste(c(res$stdout, res$stderr), collapse = "\n")
    )
  }
  if (!file.exists(out) || file.info(out)$size < 1) {
    return(data.frame(
      query = character(), target = character(), pident = numeric(),
      alnlen = integer(), evalue = numeric(), stringsAsFactors = FALSE
    ))
  }
  tab <- utils::read.delim(
    out,
    header = FALSE,
    stringsAsFactors = FALSE,
    col.names = c("query", "target", "pident", "alnlen", "mism", "gaps")
  )
  data.frame(
    query = as.character(tab$query),
    target = as.character(tab$target),
    pident = as.numeric(tab$pident),
    alnlen = as.integer(tab$alnlen),
    evalue = NA_real_,
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
.run_blastn_search <- function(query_fa, target_fa, tmp) {
  blastn <- Sys.which("blastn")
  makeblastdb <- Sys.which("makeblastdb")
  db <- file.path(tmp, "blastdb")
  out <- file.path(tmp, "blast.tsv")
  mk <- .system2_safe(
    makeblastdb,
    c("-in", target_fa, "-dbtype", "nucl", "-out", db)
  )
  if (!identical(mk$status, 0L)) {
    stop(
      "makeblastdb failed: ",
      paste(c(mk$stdout, mk$stderr), collapse = "\n")
    )
  }
  args <- c(
    "-query", query_fa,
    "-db", db,
    "-task", "blastn",
    # Avoid spaces in argv: Unix system2 joins args with spaces (no shell quoting).
    # Format 6 defaults: qseqid sseqid pident length mismatch gapopen
    #                     qstart qend sstart send evalue bitscore
    "-outfmt", "6",
    "-max_hsps", "1",
    "-max_target_seqs", "50",
    "-out", out
  )
  res <- .system2_safe(blastn, args)
  if (!identical(res$status, 0L)) {
    stop(
      "blastn failed: ",
      paste(c(res$stdout, res$stderr), collapse = "\n")
    )
  }
  if (!file.exists(out) || file.info(out)$size < 1) {
    return(data.frame(
      query = character(), target = character(), pident = numeric(),
      alnlen = integer(), evalue = numeric(), stringsAsFactors = FALSE
    ))
  }
  tab <- utils::read.delim(
    out,
    header = FALSE,
    stringsAsFactors = FALSE,
    col.names = c(
      "query", "target", "pident", "alnlen", "mismatch", "gapopen",
      "qstart", "qend", "sstart", "send", "evalue", "bitscore"
    )
  )
  data.frame(
    query = as.character(tab$query),
    target = as.character(tab$target),
    pident = as.numeric(tab$pident),
    alnlen = as.integer(tab$alnlen),
    evalue = as.numeric(tab$evalue),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
.run_mmseqs_easy_search <- function(query_fa, target_fa, tmp) {
  mmseqs <- Sys.which("mmseqs")
  out <- file.path(tmp, "mmseqs.m8")
  tmp_mm <- file.path(tmp, "mmtmp")
  dir.create(tmp_mm, showWarnings = FALSE)
  args <- c(
    "easy-search",
    query_fa,
    target_fa,
    out,
    tmp_mm,
    "--format-output",
    "query,target,fident,alnlen,mismatch,gapopen,evalue,bits",
    "-s", "7.5",
    "--search-type", "3"
  )
  res <- .system2_safe(mmseqs, args)
  if (!identical(res$status, 0L)) {
    stop(
      "mmseqs easy-search failed: ",
      paste(c(res$stdout, res$stderr), collapse = "\n")
    )
  }
  if (!file.exists(out) || file.info(out)$size < 1) {
    return(data.frame(
      query = character(), target = character(), pident = numeric(),
      alnlen = integer(), evalue = numeric(), stringsAsFactors = FALSE
    ))
  }
  tab <- utils::read.delim(
    out,
    header = FALSE,
    stringsAsFactors = FALSE,
    col.names = c(
      "query", "target", "fident", "alnlen", "mismatch", "gapopen",
      "evalue", "bits"
    )
  )
  # MMseqs fident is typically fraction in [0,1]; convert to percent.
  pident <- as.numeric(tab$fident)
  if (length(pident) > 0L && isTRUE(max(pident, na.rm = TRUE) <= 1.0001)) {
    pident <- pident * 100
  }
  data.frame(
    query = as.character(tab$query),
    target = as.character(tab$target),
    pident = pident,
    alnlen = as.integer(tab$alnlen),
    evalue = as.numeric(tab$evalue),
    stringsAsFactors = FALSE
  )
}

################################################################################
# Near-duplicate grouping (diagnostic only) + U1 loading
################################################################################

#' Connected components of consensus pairs that pass the Q1 highlight threshold.
#'
#' Groups are for **review only** (near-duplicate candidates). No sequences are
#' merged or rewritten.
#'
#' @keywords internal
.consensus_near_duplicate_groups <- function(consensus, pairwise) {
  empty_groups <- data.frame(
    index_pair_id = character(),
    gene_id = character(),
    group_id = integer(),
    member_cluster_ids = character(),
    sum_n_reads = integer(),
    highest_n_reads_cluster_id = integer(),
    stringsAsFactors = FALSE
  )
  if (nrow(consensus) < 1L) return(empty_groups)

  group_rows <- list()
  next_group <- 1L

  for (gid in unique(consensus$gene_id)) {
    sub <- consensus[consensus$gene_id == gid, , drop = FALSE]
    cids <- as.integer(sub$cluster_id)
    lab <- stats::setNames(seq_along(cids), as.character(cids))
    edges <- pairwise[
      pairwise$gene_id == gid & pairwise$pass_threshold %in% TRUE,
      ,
      drop = FALSE
    ]
    if (nrow(edges) > 0L) {
      for (e in seq_len(nrow(edges))) {
        ia <- match(edges$cluster_i[[e]], cids)
        ib <- match(edges$cluster_j[[e]], cids)
        if (is.na(ia) || is.na(ib)) next
        la <- lab[[ia]]
        lb <- lab[[ib]]
        if (!identical(la, lb)) {
          lab[lab == lb] <- la
        }
      }
    }
    for (ul in unique(lab)) {
      member_idx <- which(lab == ul)
      members <- sort(cids[member_idx])
      member_rows <- sub[sub$cluster_id %in% members, , drop = FALSE]
      ord <- order(
        -as.integer(member_rows$n_reads),
        as.integer(member_rows$cluster_id)
      )
      top_row <- member_rows[ord[[1]], , drop = FALSE]
      sum_n <- sum(as.integer(member_rows$n_reads), na.rm = TRUE)
      group_rows[[length(group_rows) + 1L]] <- data.frame(
        index_pair_id = as.character(top_row$index_pair_id[[1]]),
        gene_id = gid,
        group_id = next_group,
        member_cluster_ids = paste(members, collapse = ","),
        sum_n_reads = as.integer(sum_n),
        highest_n_reads_cluster_id = as.integer(top_row$cluster_id[[1]]),
        stringsAsFactors = FALSE
      )
      next_group <- next_group + 1L
    }
  }

  groups <- do.call(rbind, group_rows)
  rownames(groups) <- NULL
  groups
}

#' @keywords internal
.load_or_infer_u1_reads <- function(index_pair_id,
                                    amplicon_pair_dir,
                                    gene_assign,
                                    consensus,
                                    primers,
                                    max_primer_edit,
                                    read_assign_max_edit,
                                    n_core) {
  empty <- data.frame(
    index_pair_id = character(),
    gene_id = character(),
    read_id = character(),
    reason = character(),
    seq = character(),
    layer = character(),
    stringsAsFactors = FALSE
  )
  u_path <- file.path(amplicon_pair_dir, "unassigned_to_cluster.tsv")
  if (file.exists(u_path) && !is.null(gene_assign)) {
    u <- .ensure_index_pair_col(utils::read.delim(u_path, stringsAsFactors = FALSE))
    u <- u[u$index_pair_id == index_pair_id, , drop = FALSE]
    if (nrow(u) < 1L) return(empty)
    ga <- gene_assign[
      gene_assign$index_pair_id == index_pair_id &
        gene_assign$read_id %in% u$read_id,
      ,
      drop = FALSE
    ]
    seq_map <- .oriented_inserts_for_reads(ga, primers, max_primer_edit, n_core)
    out <- merge(u, seq_map, by = c("gene_id", "read_id"), all.x = TRUE)
    out$layer <- "U1"
    out$index_pair_id <- index_pair_id
    return(out[, c("index_pair_id", "gene_id", "read_id", "reason", "seq", "layer")])
  }

  if (is.null(gene_assign) || nrow(consensus) < 1L) return(empty)
  ga <- gene_assign[
    gene_assign$index_pair_id == index_pair_id &
      .is_processable_gene_row(gene_assign$gene_id, gene_assign$assign_status),
    ,
    drop = FALSE
  ]
  if (nrow(ga) < 1L) return(empty)
  seq_map <- .oriented_inserts_for_reads(ga, primers, max_primer_edit, n_core)
  rows <- list()
  for (i in seq_len(nrow(seq_map))) {
    gid <- seq_map$gene_id[[i]]
    cons <- consensus[consensus$gene_id == gid, , drop = FALSE]
    if (nrow(cons) < 1L) {
      rows[[length(rows) + 1L]] <- data.frame(
        index_pair_id = index_pair_id,
        gene_id = gid,
        read_id = seq_map$read_id[[i]],
        reason = "inferred_no_consensus",
        seq = seq_map$seq[[i]],
        layer = "U1",
        stringsAsFactors = FALSE
      )
      next
    }
    best <- Inf
    for (j in seq_len(nrow(cons))) {
      d <- edlib_edit_distance_cpp(
        seq_map$seq[[i]],
        cons$seq[[j]],
        max_edit = as.integer(read_assign_max_edit)
      )
      if (!is.na(d) && d < best) best <- d
    }
    if (!is.finite(best) || best > as.integer(read_assign_max_edit)) {
      rows[[length(rows) + 1L]] <- data.frame(
        index_pair_id = index_pair_id,
        gene_id = gid,
        read_id = seq_map$read_id[[i]],
        reason = "inferred_far_from_consensus",
        seq = seq_map$seq[[i]],
        layer = "U1",
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(rows)) return(empty)
  do.call(rbind, rows)
}

#' @keywords internal
.oriented_inserts_for_reads <- function(gene_reads, primers, max_primer_edit,
                                        n_core) {
  if (nrow(gene_reads) < 1L) {
    return(data.frame(
      gene_id = character(),
      read_id = character(),
      seq = character(),
      stringsAsFactors = FALSE
    ))
  }
  parts <- lapply(split(gene_reads, gene_reads$gene_id), function(gr) {
    prepared <- .prepare_gene_seqs_for_cluster(
      gene_reads = gr,
      primers = primers,
      gene_id = gr$gene_id[[1]],
      max_primer_edit = max_primer_edit,
      n_core = n_core,
      strict_end_trim = TRUE
    )
    data.frame(
      gene_id = gr$gene_id[[1]],
      read_id = prepared$read_id,
      seq = prepared$seqs,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, parts)
}
