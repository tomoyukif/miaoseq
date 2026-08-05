################################################################################
# Demultiplexing (C++ edlib core: suffix → barcode; R wrapper for I/O)
################################################################################

#' Expand "~" and normalize for C++ file I/O (gzopen/ofstream do not expand "~").
#' @keywords internal
.abs_path <- function(path, mustWork = FALSE) {
    if (is.null(path)) {
        return(NULL)
    }
    path <- path.expand(as.character(path))
    if (!length(path)) {
        return(path)
    }
    normalizePath(path, winslash = "/", mustWork = mustWork)
}

#' Demultiplex reads by dual barcodes
#'
#' Classifies multiplexed FASTQ reads into samples using a two-step procedure
#' implemented in C++ (edlib HW for adapter **suffix** anchors; barcode matching
#' via a precomputed mutant dictionary). FASTQ is streamed in chunks (no full
#' Biostrings load). Assignment tables are written incrementally; when
#' `split_reads = TRUE`, per-sample FASTQ files are written in the same pass.
#'
#' @param fastq Character vector of FASTQ paths (`.fq` / `.fastq`, optionally `.gz`).
#' @param demult_dir Output directory for demultiplex results.
#' @param index_list Path to a 5-column headerless CSV:
#'   `index_pair_id`, forward index ID, forward index sequence,
#'   reverse index ID, reverse index sequence.
#' @param sample_list Optional path to a CSV with at least
#'   `index_pair_id` and `sample_name` (headerless or with header).
#'   If omitted, `sample_id` equals `index_pair_id`.
#' @param n_core Number of OpenMP threads for the C++ core (`1` disables).
#' @param end_window Bases from each read end to search for anchors.
#' @param max_anchor_edit Maximum edit distance for suffix detection.
#' @param max_barcode_edit Maximum edit distance for barcode mutant dictionary.
#' @param allow_revcomp If `TRUE`, also search reverse-complement windows.
#' @param allow_single_end If `TRUE`, assign reads when only one end (F or R)
#'   yields a unique barcode hit. Allowed only when every F and every R barcode
#'   ID appears in exactly one index pair (non-combinatorial layout). On
#'   combinatorial plates (shared F/R IDs across wells), the run stops at
#'   startup. Conflicting hits on both ends are rejected as `ambiguous_ends`.
#' @param split_reads If `TRUE`, write per-sample FASTQ under
#'   `demult_dir/by_sample/` in the same streaming pass.
#' @param compress If `TRUE` and `split_reads = TRUE`, write `.fq.gz`.
#' @param chunk_size Reads per C++ batch (streaming to limit peak memory).
#' @param stats_unassign If `TRUE`, write `stats_unassigned.tsv` summarizing
#'   unassigned read counts by `reason`.
#' @param return_tables If `TRUE` (default), read `assignments.tsv` /
#'   `unassigned.tsv` back into data.frames for the return value. Set
#'   `FALSE` to keep peak RAM low (tables remain on disk).
#'
#' @return A list with `assignments`, `summary`, `unassigned`, `demult_dir`,
#'   and `n_edlib` (total edlib calls in the C++ core). When
#'   `stats_unassign = TRUE`, also includes `stats_unassigned`. When
#'   `return_tables = FALSE`, `assignments` / `unassigned` are empty
#'   placeholders and paths are in `assignments_tsv` / `unassigned_tsv`.
#'
#' @export
#' @useDynLib miaoseq, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats setNames
#' @importFrom utils read.csv write.table
doDemultiplex <- function(fastq,
                          demult_dir,
                          index_list,
                          sample_list = NULL,
                          n_core = 1,
                          end_window = 120,
                          max_anchor_edit = 10,
                          max_barcode_edit = 2,
                          allow_revcomp = TRUE,
                          allow_single_end = FALSE,
                          split_reads = FALSE,
                          compress = TRUE,
                          chunk_size = 20000,
                          stats_unassign = FALSE,
                          return_tables = TRUE) {
    # C++ ofstream/gzopen do not expand "~"; normalize before streaming I/O.
    demult_dir <- .abs_path(demult_dir[[1L]], mustWork = FALSE)
    index_list <- .abs_path(index_list[[1L]], mustWork = TRUE)
    if (!is.null(sample_list)) {
        sample_list <- .abs_path(sample_list[[1L]], mustWork = TRUE)
    }
    if (!dir.exists(demult_dir)) {
        dir.create(demult_dir, recursive = TRUE, showWarnings = FALSE)
    }
    demult_dir <- .abs_path(demult_dir, mustWork = TRUE)
    fastq <- .abs_path(fastq, mustWork = TRUE)
    if (length(fastq) < 1) {
        stop("All fastq paths must exist.")
    }

    layout <- .parse_index_layout(index_list)
    .assert_single_end_allowed(layout, allow_single_end)
    sample_map <- .load_sample_map(sample_list, layout$sample_map)
    design <- .validate_barcode_design(layout, max_barcode_edit)

    # diagnostic dict build (also warms logic; core rebuilds maps internally)
    f_diag <- build_barcode_dict_cpp(unname(layout$f_barcodes),
                                     names(layout$f_barcodes),
                                     as.integer(max_barcode_edit))
    r_diag <- build_barcode_dict_cpp(unname(layout$r_barcodes),
                                     names(layout$r_barcodes),
                                     as.integer(max_barcode_edit))

    .write_index_layout(layout, file.path(demult_dir, "index_layout.tsv"))
    .conflict_rows <- function(side, diag) {
        n <- length(diag$conflict_mutants)
        if (!n) {
            return(data.frame(
                side = character(),
                mutant = character(),
                parents = character(),
                stringsAsFactors = FALSE
            ))
        }
        data.frame(
            side = rep(side, n),
            mutant = diag$conflict_mutants,
            parents = diag$conflict_parents,
            stringsAsFactors = FALSE
        )
    }
    write.table(
        rbind(
            .conflict_rows("F", f_diag),
            .conflict_rows("R", r_diag)
        ),
        file.path(demult_dir, "barcode_conflicts.tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE
    )
    write.table(design, file.path(demult_dir, "design_check.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)

    if (isTRUE(split_reads)) {
        dir.create(file.path(demult_dir, "by_sample"),
                   recursive = TRUE, showWarnings = FALSE)
    }

    message("Demultiplexing ", length(fastq), " FASTQ file(s) (streaming)...")
    res <- demux_fastq_stream_cpp(
        fastq_files = fastq,
        demult_dir = demult_dir,
        f_suffix = layout$f_suffix,
        r_suffix = layout$r_suffix,
        f_barcodes = unname(layout$f_barcodes),
        f_barcode_names = names(layout$f_barcodes),
        r_barcodes = unname(layout$r_barcodes),
        r_barcode_names = names(layout$r_barcodes),
        pair_f_names = sample_map$f_index_id,
        pair_r_names = sample_map$r_index_id,
        pair_ids = sample_map$index_pair_id,
        sample_ids = sample_map$sample_id,
        end_window = as.integer(end_window),
        max_anchor_edit = as.integer(max_anchor_edit),
        max_barcode_edit = as.integer(max_barcode_edit),
        allow_revcomp = isTRUE(allow_revcomp),
        allow_single_end = isTRUE(allow_single_end),
        n_core = as.integer(n_core),
        chunk_size = as.integer(chunk_size),
        split_reads = isTRUE(split_reads),
        compress = isTRUE(compress),
        include_unassigned_fastq = FALSE
    )

    summary_df <- data.frame(
        sample_id = as.character(res$summary_sample_id),
        index_pair_id = as.character(res$summary_index_pair_id),
        n_reads = as.integer(res$summary_n_reads),
        n_complete = as.integer(res$summary_n_complete),
        n_fuzzy = as.integer(res$summary_n_fuzzy),
        n_single_end = as.integer(res$summary_n_single_end),
        stringsAsFactors = FALSE
    )
    summary_df <- summary_df[summary_df$n_reads > 0L, , drop = FALSE]
    rownames(summary_df) <- NULL
    write.table(summary_df, file.path(demult_dir, "summary_by_sample.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)

    stats_unassigned_df <- NULL
    if (isTRUE(stats_unassign)) {
        rn <- as.character(res$reason_names)
        rc <- as.integer(res$reason_n)
        if (length(rn)) {
            n_total <- sum(rc)
            stats_unassigned_df <- data.frame(
                scope = "overall",
                sample_id = "",
                reason = rn,
                n = rc,
                fraction = as.numeric(rc) / n_total,
                stringsAsFactors = FALSE
            )
            stats_unassigned_df <-
                stats_unassigned_df[order(-stats_unassigned_df$n), , drop = FALSE]
            rownames(stats_unassigned_df) <- NULL
        } else {
            stats_unassigned_df <- .stats_unassign_reason_table(
                character(), scope = "overall"
            )
        }
        write.table(
            stats_unassigned_df,
            file.path(demult_dir, "stats_unassigned.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE
        )
    }

    assignments <- .empty_assignments()
    unassigned <- .empty_unassigned()
    if (isTRUE(return_tables)) {
        assignments <- utils::read.delim(
            res$assignments_tsv, stringsAsFactors = FALSE
        )
        unassigned <- utils::read.delim(
            res$unassigned_tsv, stringsAsFactors = FALSE
        )
    }

    message(
        "Demultiplex done: ",
        format(as.integer(res$n_assigned), big.mark = ","), " assigned / ",
        format(as.integer(res$n_records), big.mark = ","), " reads"
    )

    list(
        assignments = assignments,
        summary = summary_df,
        unassigned = unassigned,
        demult_dir = demult_dir,
        n_edlib = res$n_edlib,
        stats_unassigned = stats_unassigned_df,
        assignments_tsv = res$assignments_tsv,
        unassigned_tsv = res$unassigned_tsv,
        n_records = res$n_records,
        n_assigned = res$n_assigned,
        n_unassigned = res$n_unassigned
    )
}

#' Split a FASTQ by demultiplex assignments
#'
#' Writes one FASTQ per `sample_id` using `assignments` from [doDemultiplex()].
#' Can be run after `doDemultiplex(..., split_reads = FALSE)`. Uses a C++
#' streaming reader/writer (not line-by-line R I/O).
#'
#' @param fastq Character vector of source FASTQ paths.
#' @param assignments `assignments.tsv` path or a data.frame with
#'   `read_id` and `sample_id` columns.
#' @param out_dir Output directory for per-sample FASTQ files.
#' @param compress If `TRUE`, write `.fq.gz`.
#' @param include_unassigned If `TRUE`, also write `unassigned.fq(.gz)`.
#' @param unassigned Optional unassigned table / TSV path (required when
#'   `include_unassigned = TRUE` unless `assignments` is a TSV path next to
#'   `unassigned.tsv`).
#' @param n_core Unused (reserved).
#'
#' @return Invisibly, a named integer vector of reads written per sample.
#'
#' @export
#' @importFrom utils read.delim write.table
splitDemultiplexReads <- function(fastq,
                                  assignments,
                                  out_dir,
                                  compress = TRUE,
                                  include_unassigned = FALSE,
                                  unassigned = NULL,
                                  n_core = 1) {
    invisible(.split_fastq_by_assignment(
        fastq = fastq,
        assignments = assignments,
        out_dir = out_dir,
        compress = compress,
        include_unassigned = include_unassigned,
        unassigned = unassigned
    ))
}

################################################################################
# Layout / diagnostics (R)
################################################################################

#' @keywords internal
.parse_index_layout <- function(index_list) {
    index_df <- utils::read.csv(index_list, header = FALSE, stringsAsFactors = FALSE)
    if (ncol(index_df) < 5) stop("index_list must have 5 columns.")
    names(index_df)[1:5] <- c(
        "index_pair_id", "f_index_id", "f_seq", "r_index_id", "r_seq"
    )
    index_df$f_seq <- toupper(gsub("[^ACGT]", "", index_df$f_seq))
    index_df$r_seq <- toupper(gsub("[^ACGT]", "", index_df$r_seq))

    f_parts <- .decompose_index_seqs(index_df$f_seq[!duplicated(index_df$f_index_id)])
    r_parts <- .decompose_index_seqs(index_df$r_seq[!duplicated(index_df$r_index_id)])

    f_ids <- index_df$f_index_id[!duplicated(index_df$f_index_id)]
    r_ids <- index_df$r_index_id[!duplicated(index_df$r_index_id)]
    f_barcodes <- setNames(f_parts$barcodes, f_ids)
    r_barcodes <- setNames(r_parts$barcodes, r_ids)

    sample_map <- data.frame(
        index_pair_id = index_df$index_pair_id,
        f_index_id = index_df$f_index_id,
        r_index_id = index_df$r_index_id,
        stringsAsFactors = FALSE
    )
    sample_map$pair_key <- paste(sample_map$f_index_id, sample_map$r_index_id, sep = "\t")

    list(
        f_prefix = f_parts$prefix,
        f_suffix = f_parts$suffix,
        f_barcode_len = nchar(f_parts$barcodes[[1]]),
        r_prefix = r_parts$prefix,
        r_suffix = r_parts$suffix,
        r_barcode_len = nchar(r_parts$barcodes[[1]]),
        f_barcodes = f_barcodes,
        r_barcodes = r_barcodes,
        sample_map = sample_map
    )
}

#' @keywords internal
.decompose_index_seqs <- function(seqs) {
    seqs <- as.character(seqs)
    if (length(seqs) < 1) stop("No index sequences to decompose.")
    prefix_len <- .common_prefix_len(seqs)
    suffix_len <- .common_suffix_len(seqs)
    if (prefix_len < 1 || suffix_len < 1) {
        stop("Could not detect common prefix/suffix in index sequences.")
    }
    mid_lens <- nchar(seqs) - prefix_len - suffix_len
    if (any(mid_lens <= 0) || length(unique(mid_lens)) != 1) {
        stop("Index sequences do not share a uniform barcode middle region.")
    }
    barcodes <- substr(seqs, prefix_len + 1, prefix_len + mid_lens[1])
    list(
        prefix = substr(seqs[[1]], 1, prefix_len),
        suffix = substr(seqs[[1]], nchar(seqs[[1]]) - suffix_len + 1, nchar(seqs[[1]])),
        barcodes = barcodes
    )
}

#' @keywords internal
.common_prefix_len <- function(seqs) {
    if (length(seqs) == 1) return(0L)
    chars <- strsplit(seqs, "", fixed = TRUE)
    n <- min(vapply(chars, length, integer(1)))
    k <- 0L
    for (i in seq_len(n)) {
        bases <- vapply(chars, `[[`, character(1), i)
        if (length(unique(bases)) == 1) k <- i else break
    }
    k
}

#' @keywords internal
.common_suffix_len <- function(seqs) {
    if (length(seqs) == 1) return(0L)
    chars <- lapply(strsplit(seqs, "", fixed = TRUE), rev)
    n <- min(vapply(chars, length, integer(1)))
    k <- 0L
    for (i in seq_len(n)) {
        bases <- vapply(chars, `[[`, character(1), i)
        if (length(unique(bases)) == 1) k <- i else break
    }
    k
}

#' @keywords internal
#' Stop early when allow_single_end is incompatible with combinatorial indexes.
.assert_single_end_allowed <- function(layout, allow_single_end) {
    if (!isTRUE(allow_single_end)) {
        return(invisible(NULL))
    }
    sm <- layout$sample_map
    f_reuse <- any(duplicated(sm$f_index_id))
    r_reuse <- any(duplicated(sm$r_index_id))
    if (!f_reuse && !r_reuse) {
        return(invisible(NULL))
    }
    f_max <- max(as.integer(table(sm$f_index_id)))
    r_max <- max(as.integer(table(sm$r_index_id)))
    msg <- paste0(
        "allow_single_end=TRUE is not allowed for combinatorial index layouts ",
        "where the same F and/or R barcode ID is reused across multiple wells. ",
        "Detected max F reuse=", f_max, ", max R reuse=", r_max, ". ",
        "A single-end barcode hit cannot uniquely identify a sample on such ",
        "plates; refusing to run to avoid silent cross-sample assignment. ",
        "Use allow_single_end=FALSE, or supply a 1:1 (non-reused) index_list."
    )
    warning(msg, call. = FALSE)
    stop(msg, call. = FALSE)
}

#' @keywords internal
.validate_barcode_design <- function(layout, max_barcode_edit) {
    # Light check without edlibR: Hamming on equal-length barcodes as lower bound proxy
    f_min <- .min_hamming(layout$f_barcodes)
    r_min <- .min_hamming(layout$r_barcodes)
    need <- 2L * as.integer(max_barcode_edit) + 1L
    data.frame(
        side = c("F", "R"),
        n_barcodes = c(length(layout$f_barcodes), length(layout$r_barcodes)),
        min_hamming = c(f_min, r_min),
        recommended_min_edit = c(need, need),
        note = c(
            "Hamming lower-bound; C++ dict applies collision removal",
            "Hamming lower-bound; C++ dict applies collision removal"
        ),
        stringsAsFactors = FALSE
    )
}

#' @keywords internal
.min_hamming <- function(seqs) {
    seqs <- unique(as.character(seqs))
    if (length(seqs) < 2) return(NA_integer_)
    best <- Inf
    for (i in seq_len(length(seqs) - 1L)) {
        for (j in (i + 1L):length(seqs)) {
            a <- strsplit(seqs[i], "", fixed = TRUE)[[1]]
            b <- strsplit(seqs[j], "", fixed = TRUE)[[1]]
            n <- min(length(a), length(b))
            d <- sum(a[seq_len(n)] != b[seq_len(n)]) + abs(length(a) - length(b))
            if (d < best) best <- d
        }
    }
    as.integer(best)
}

#' @keywords internal
.load_sample_map <- function(sample_list, sample_map) {
    sample_map$sample_id <- sample_map$index_pair_id
    if (is.null(sample_list)) return(sample_map)
    sl <- utils::read.csv(sample_list, header = FALSE, stringsAsFactors = FALSE)
    if (ncol(sl) < 2) stop("sample_list needs at least 2 columns.")
    if (identical(tolower(sl[1, 1]), "index_pair_id") ||
        identical(tolower(sl[1, 1]), "index_pair")) {
        sl <- sl[-1, , drop = FALSE]
    }
    names(sl)[1:2] <- c("index_pair_id", "sample_name")
    hit <- match(sample_map$index_pair_id, sl$index_pair_id)
    sample_map$sample_id[!is.na(hit)] <- sl$sample_name[hit[!is.na(hit)]]
    sample_map
}

#' @keywords internal
.write_index_layout <- function(layout, path) {
    df <- data.frame(
        side = c("F", "R"),
        prefix = c(layout$f_prefix, layout$r_prefix),
        suffix = c(layout$f_suffix, layout$r_suffix),
        barcode_len = c(layout$f_barcode_len, layout$r_barcode_len),
        n_barcodes = c(length(layout$f_barcodes), length(layout$r_barcodes)),
        stringsAsFactors = FALSE
    )
    write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

################################################################################
# C++ driver (legacy in-memory path kept for demux_reads_cpp unit use)
################################################################################

#' @keywords internal
.demultiplex_fastq_cpp <- function(fq, layout, sample_map, n_core,
                                   end_window, max_anchor_edit,
                                   max_barcode_edit,
                                   allow_revcomp, allow_single_end, chunk_size) {
    # Prefer streaming API via doDemultiplex(); this helper remains for
    # callers that already hold sequences in memory.
    reads <- .read_fastq_seqs(fq)
    n <- length(reads)
    if (n < 1) {
        return(list(assignments = .empty_assignments(),
                    unassigned = .empty_unassigned(),
                    n_edlib = 0))
    }
    ids <- names(reads)
    seqs <- as.character(reads)

    starts <- seq(1L, n, by = as.integer(chunk_size))
    assign_rows <- list()
    un_rows <- list()
    n_edlib <- 0

    for (s in starts) {
        e <- min(s + as.integer(chunk_size) - 1L, n)
        res <- demux_reads_cpp(
            seqs = seqs[s:e],
            ids = ids[s:e],
            f_suffix = layout$f_suffix,
            r_suffix = layout$r_suffix,
            f_barcodes = unname(layout$f_barcodes),
            f_barcode_names = names(layout$f_barcodes),
            r_barcodes = unname(layout$r_barcodes),
            r_barcode_names = names(layout$r_barcodes),
            pair_f_names = sample_map$f_index_id,
            pair_r_names = sample_map$r_index_id,
            pair_ids = sample_map$index_pair_id,
            sample_ids = sample_map$sample_id,
            end_window = as.integer(end_window),
            max_anchor_edit = as.integer(max_anchor_edit),
            max_barcode_edit = as.integer(max_barcode_edit),
            allow_revcomp = isTRUE(allow_revcomp),
            allow_single_end = isTRUE(allow_single_end),
            n_core = as.integer(n_core)
        )
        n_edlib <- n_edlib + res$n_edlib
        assigned <- which(res$status == 1L)
        failed <- which(res$status == 0L)
        if (length(assigned)) {
            assign_rows[[length(assign_rows) + 1L]] <- data.frame(
                read_id = res$read_id[assigned],
                index_pair_id = res$index_pair_id[assigned],
                f_index_id = res$f_index_id[assigned],
                r_index_id = res$r_index_id[assigned],
                sample_id = res$sample_id[assigned],
                source_file = fq,
                barcode_edit_f = res$barcode_edit_f[assigned],
                barcode_edit_r = res$barcode_edit_r[assigned],
                match_class = res$match_class[assigned],
                assign_mode = res$assign_mode[assigned],
                stringsAsFactors = FALSE
            )
        }
        if (length(failed)) {
            un_rows[[length(un_rows) + 1L]] <- data.frame(
                read_id = res$read_id[failed],
                reason = res$reason[failed],
                source_file = fq,
                stringsAsFactors = FALSE
            )
        }
    }

    list(
        assignments = .rbind_or_empty(assign_rows, .empty_assignments()),
        unassigned = .rbind_or_empty(un_rows, .empty_unassigned()),
        n_edlib = n_edlib
    )
}

#' @keywords internal
.read_fastq_seqs <- function(fastq) {
    ext <- tolower(fastq)
    if (grepl("\\.(fq|fastq)(\\.gz)?$", ext) || grepl("\\.gz$", ext)) {
        ss <- Biostrings::readDNAStringSet(fastq, format = "fastq")
    } else {
        ss <- Biostrings::readDNAStringSet(fastq)
    }
    names(ss) <- vapply(strsplit(names(ss), "\\s+"), `[[`, character(1), 1L)
    ss
}

################################################################################
# Output helpers
################################################################################

#' @keywords internal
.empty_assignments <- function() {
    data.frame(
        read_id = character(),
        index_pair_id = character(),
        f_index_id = character(),
        r_index_id = character(),
        sample_id = character(),
        source_file = character(),
        barcode_edit_f = integer(),
        barcode_edit_r = integer(),
        match_class = character(),
        assign_mode = character(),
        stringsAsFactors = FALSE
    )
}

#' @keywords internal
.empty_unassigned <- function() {
    data.frame(
        read_id = character(),
        reason = character(),
        source_file = character(),
        stringsAsFactors = FALSE
    )
}

#' @keywords internal
.rbind_or_empty <- function(parts, empty) {
    parts <- Filter(function(x) !is.null(x) && nrow(x) > 0, parts)
    if (!length(parts)) return(empty)
    out <- do.call(rbind, parts)
    rownames(out) <- NULL
    out
}

#' @keywords internal
.stats_unassign_reason_table <- function(reasons, scope, sample_id = NA_character_) {
    reasons <- as.character(reasons)
    reasons <- reasons[!is.na(reasons) & nzchar(reasons)]
    if (!length(reasons)) {
        return(data.frame(
            scope = character(),
            sample_id = character(),
            reason = character(),
            n = integer(),
            fraction = numeric(),
            stringsAsFactors = FALSE
        ))
    }
    tab <- sort(table(reasons), decreasing = TRUE)
    n_total <- sum(tab)
    data.frame(
        scope = scope,
        sample_id = if (is.na(sample_id)) "" else sample_id,
        reason = names(tab),
        n = as.integer(tab),
        fraction = as.numeric(tab) / n_total,
        stringsAsFactors = FALSE,
        row.names = NULL
    )
}

#' @keywords internal
.stats_unassign_demultiplex <- function(unassigned) {
    parts <- list(
        .stats_unassign_reason_table(unassigned$reason, scope = "overall")
    )
    if ("source_file" %in% names(unassigned) && nrow(unassigned) > 0L) {
        by_file <- split(unassigned$reason, unassigned$source_file)
        parts <- c(
            parts,
            lapply(names(by_file), function(src) {
                .stats_unassign_reason_table(by_file[[src]], scope = "source_file", sample_id = src)
            })
        )
    }
    do.call(rbind, parts)
}

#' @keywords internal
.stats_unassign_amplicon <- function(gene_assignments) {
    if (nrow(gene_assignments) < 1L) {
        return(.stats_unassign_reason_table(character(), scope = "overall"))
    }
    ok_status <- "assigned"
    is_unassigned <- is.na(gene_assignments$gene_id) |
        gene_assignments$gene_id == "" |
        !(gene_assignments$assign_status %in% ok_status)
    un_df <- gene_assignments[is_unassigned, , drop = FALSE]
    parts <- list(
        .stats_unassign_reason_table(un_df$assign_status, scope = "overall")
    )
    if (nrow(un_df) > 0L && "sample_id" %in% names(un_df)) {
        by_sample <- split(un_df$assign_status, un_df$sample_id)
        parts <- c(
            parts,
            lapply(names(by_sample), function(sid) {
                .stats_unassign_reason_table(by_sample[[sid]], scope = "sample", sample_id = sid)
            })
        )
    }
    do.call(rbind, parts)
}

#' @keywords internal
.summarize_assignments <- function(assignments) {
    if (nrow(assignments) < 1) {
        return(data.frame(
            sample_id = character(),
            index_pair_id = character(),
            n_reads = integer(),
            n_complete = integer(),
            n_fuzzy = integer(),
            n_single_end = integer(),
            stringsAsFactors = FALSE
        ))
    }
    spl <- split(assignments, assignments$sample_id)
    do.call(rbind, lapply(spl, function(df) {
        data.frame(
            sample_id = df$sample_id[1],
            index_pair_id = df$index_pair_id[1],
            n_reads = nrow(df),
            n_complete = sum(df$match_class == "complete_match"),
            n_fuzzy = sum(df$match_class == "fuzzy_match"),
            n_single_end = sum(df$assign_mode %in% c("single_f", "single_r")),
            stringsAsFactors = FALSE
        )
    }))
}

#' @keywords internal
.write_demultiplex_tables <- function(demult_dir, assignments, summary_df, unassigned) {
    write.table(assignments, file.path(demult_dir, "assignments.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(summary_df, file.path(demult_dir, "summary_by_sample.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(unassigned, file.path(demult_dir, "unassigned.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
}

################################################################################
# split FASTQ
################################################################################

#' @keywords internal
.split_fastq_by_assignment <- function(fastq,
                                       assignments,
                                       out_dir,
                                       compress = TRUE,
                                       include_unassigned = FALSE,
                                       unassigned = NULL) {
    out_dir <- .abs_path(out_dir[[1L]], mustWork = FALSE)
    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }
    out_dir <- .abs_path(out_dir, mustWork = TRUE)
    fastq <- .abs_path(fastq, mustWork = TRUE)

    assign_tsv <- NULL
    un_tsv <- ""
    tmp_assign <- NULL
    tmp_un <- NULL
    on.exit({
        if (!is.null(tmp_assign) && file.exists(tmp_assign)) unlink(tmp_assign)
        if (!is.null(tmp_un) && file.exists(tmp_un)) unlink(tmp_un)
    }, add = TRUE)

    if (is.character(assignments) && length(assignments) == 1L &&
        file.exists(.abs_path(assignments, mustWork = FALSE))) {
        assign_tsv <- .abs_path(assignments, mustWork = TRUE)
        if (isTRUE(include_unassigned)) {
            if (is.character(unassigned) && length(unassigned) == 1L &&
                file.exists(.abs_path(unassigned, mustWork = FALSE))) {
                un_tsv <- .abs_path(unassigned, mustWork = TRUE)
            } else {
                cand <- file.path(dirname(assign_tsv), "unassigned.tsv")
                if (file.exists(cand)) {
                    un_tsv <- .abs_path(cand, mustWork = TRUE)
                }
            }
        }
    } else {
        assignments <- .as_assignment_df(assignments)
        if (!all(c("read_id", "sample_id") %in% names(assignments))) {
            stop("assignments must contain read_id and sample_id columns.")
        }
        tmp_assign <- tempfile(fileext = ".tsv")
        utils::write.table(
            data.frame(
                read_id = assignments$read_id,
                index_pair_id = if ("index_pair_id" %in% names(assignments))
                    assignments$index_pair_id else assignments$sample_id,
                f_index_id = if ("f_index_id" %in% names(assignments))
                    assignments$f_index_id else "",
                r_index_id = if ("r_index_id" %in% names(assignments))
                    assignments$r_index_id else "",
                sample_id = assignments$sample_id,
                stringsAsFactors = FALSE
            ),
            tmp_assign, sep = "\t", quote = FALSE, row.names = FALSE
        )
        assign_tsv <- tmp_assign
        if (isTRUE(include_unassigned) && !is.null(unassigned)) {
            unassigned <- .as_assignment_df(unassigned)
            tmp_un <- tempfile(fileext = ".tsv")
            utils::write.table(
                data.frame(
                    read_id = unassigned$read_id,
                    reason = if ("reason" %in% names(unassigned))
                        unassigned$reason else "",
                    stringsAsFactors = FALSE
                ),
                tmp_un, sep = "\t", quote = FALSE, row.names = FALSE
            )
            un_tsv <- tmp_un
        }
    }

    res <- split_fastq_by_assignment_cpp(
        fastq_files = fastq,
        assignments_tsv = assign_tsv,
        out_dir = out_dir,
        compress = isTRUE(compress),
        include_unassigned = isTRUE(include_unassigned),
        unassigned_tsv = un_tsv
    )
    counts <- setNames(as.integer(res$n_reads), as.character(res$sample_id))
    invisible(counts)
}

#' @keywords internal
.as_assignment_df <- function(x) {
    if (is.character(x) && length(x) == 1L && file.exists(x)) {
        utils::read.delim(x, stringsAsFactors = FALSE)
    } else {
        as.data.frame(x, stringsAsFactors = FALSE)
    }
}

#' @keywords internal
.safe_filename <- function(x) {
    gsub("[^A-Za-z0-9._-]+", "_", x)
}

#' @keywords internal
.stream_fastq_write <- function(fastq, id2sample, cons, counts, un_ids, un_con) {
    .Deprecated(msg = "Internal R FASTQ splitter replaced by C++ streaming.")
    counts
}
