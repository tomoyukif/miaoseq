################################################################################
# Demultiplexing (C++ edlib core: suffix → barcode; R wrapper for I/O)
################################################################################

#' Demultiplex reads by dual barcodes
#'
#' Classifies multiplexed FASTQ reads into samples using a two-step procedure
#' implemented in C++ (edlib HW for adapter **suffix** anchors; barcode matching
#' via a precomputed mutant dictionary). Optional prefix checks are off by
#' default on the hot path. Primary outputs are assignment tables; sample FASTQ
#' files are written only when `split_reads = TRUE` (via [splitDemultiplexReads()]).
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
#' @param max_prefix_edit Maximum edit distance for optional prefix checks.
#' @param max_barcode_edit Maximum edit distance for barcode mutant dictionary.
#' @param require_prefix If `TRUE`, reject reads without a detectable prefix.
#' @param check_prefix If `TRUE`, run optional prefix labeling on every hit
#'   (slower). Ignored when `require_prefix = TRUE` (prefix is always checked).
#' @param require_unique_pair Kept for API compatibility (pair map is unique by construction).
#' @param allow_revcomp If `TRUE`, also search reverse-complement windows.
#' @param split_reads If `TRUE`, call [splitDemultiplexReads()] after assignment.
#' @param compress Passed to [splitDemultiplexReads()] when `split_reads = TRUE`.
#' @param chunk_size Reads per C++ batch (streaming to limit peak memory).
#'
#' @return A list with `assignments`, `summary`, `unassigned`, `demult_dir`,
#'   and `n_edlib` (total edlib calls in the C++ core).
#'
#' @export
#' @useDynLib miaoseq, .registration = TRUE
#' @importFrom Biostrings readDNAStringSet
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
                          max_prefix_edit = 5,
                          max_barcode_edit = 2,
                          require_prefix = FALSE,
                          check_prefix = FALSE,
                          require_unique_pair = TRUE,
                          allow_revcomp = TRUE,
                          split_reads = FALSE,
                          compress = TRUE,
                          chunk_size = 20000) {
    if (!dir.exists(demult_dir)) {
        dir.create(demult_dir, recursive = TRUE, showWarnings = FALSE)
    }
    fastq <- as.character(fastq)
    if (length(fastq) < 1 || any(!file.exists(fastq))) {
        stop("All fastq paths must exist.")
    }

    layout <- .parse_index_layout(index_list)
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
    write.table(
        rbind(
            data.frame(side = "F",
                       mutant = f_diag$conflict_mutants,
                       parents = f_diag$conflict_parents,
                       stringsAsFactors = FALSE),
            data.frame(side = "R",
                       mutant = r_diag$conflict_mutants,
                       parents = r_diag$conflict_parents,
                       stringsAsFactors = FALSE)
        ),
        file.path(demult_dir, "barcode_conflicts.tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE
    )
    write.table(design, file.path(demult_dir, "design_check.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)

    assign_parts <- list()
    un_parts <- list()
    n_edlib_total <- 0

    for (fq in fastq) {
        message("Demultiplexing: ", fq)
        part <- .demultiplex_fastq_cpp(
            fq = fq,
            layout = layout,
            sample_map = sample_map,
            n_core = n_core,
            end_window = end_window,
            max_anchor_edit = max_anchor_edit,
            max_prefix_edit = max_prefix_edit,
            max_barcode_edit = max_barcode_edit,
            require_prefix = require_prefix,
            check_prefix = check_prefix || require_prefix,
            allow_revcomp = allow_revcomp,
            chunk_size = chunk_size
        )
        assign_parts[[fq]] <- part$assignments
        un_parts[[fq]] <- part$unassigned
        n_edlib_total <- n_edlib_total + part$n_edlib
    }

    assignments <- .rbind_or_empty(assign_parts, .empty_assignments())
    unassigned <- .rbind_or_empty(un_parts, .empty_unassigned())
    summary_df <- .summarize_assignments(assignments)
    .write_demultiplex_tables(demult_dir, assignments, summary_df, unassigned)

    if (isTRUE(split_reads) && nrow(assignments) > 0) {
        splitDemultiplexReads(
            fastq = fastq,
            assignments = assignments,
            out_dir = file.path(demult_dir, "by_sample"),
            compress = compress
        )
    }

    list(
        assignments = assignments,
        summary = summary_df,
        unassigned = unassigned,
        demult_dir = demult_dir,
        n_edlib = n_edlib_total
    )
}

#' Split a FASTQ by demultiplex assignments
#'
#' Writes one FASTQ per `sample_id` using `assignments` from [doDemultiplex()].
#' Can be run after `doDemultiplex(..., split_reads = FALSE)`.
#'
#' @param fastq Character vector of source FASTQ paths.
#' @param assignments `assignments.tsv` path or a data.frame with
#'   `read_id` and `sample_id` columns.
#' @param out_dir Output directory for per-sample FASTQ files.
#' @param compress If `TRUE`, write `.fq.gz`.
#' @param include_unassigned If `TRUE`, also write `unassigned.fq(.gz)`.
#' @param unassigned Optional unassigned table / TSV path.
#' @param n_core Unused (reserved).
#'
#' @return Invisibly, a named integer vector of reads written per sample.
#'
#' @export
#' @importFrom utils read.delim
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
# C++ driver
################################################################################

#' @keywords internal
.demultiplex_fastq_cpp <- function(fq, layout, sample_map, n_core,
                                   end_window, max_anchor_edit, max_prefix_edit,
                                   max_barcode_edit, require_prefix, check_prefix,
                                   allow_revcomp, chunk_size) {
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
            f_prefix = layout$f_prefix,
            r_prefix = layout$r_prefix,
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
            max_prefix_edit = as.integer(max_prefix_edit),
            max_barcode_edit = as.integer(max_barcode_edit),
            require_prefix = isTRUE(require_prefix),
            check_prefix = isTRUE(check_prefix),
            allow_revcomp = isTRUE(allow_revcomp),
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
                anchor_status = res$anchor_status[assigned],
                prefix_edit_f = ifelse(res$prefix_edit_f[assigned] < 0, NA_integer_,
                                       res$prefix_edit_f[assigned]),
                prefix_edit_r = ifelse(res$prefix_edit_r[assigned] < 0, NA_integer_,
                                       res$prefix_edit_r[assigned]),
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
        anchor_status = character(),
        prefix_edit_f = integer(),
        prefix_edit_r = integer(),
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
.summarize_assignments <- function(assignments) {
    if (nrow(assignments) < 1) {
        return(data.frame(
            sample_id = character(),
            index_pair_id = character(),
            n_reads = integer(),
            n_complete = integer(),
            n_fuzzy = integer(),
            n_high_confidence = integer(),
            n_partial_anchor = integer(),
            n_rescued = integer(),
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
            n_high_confidence = sum(df$anchor_status == "high_confidence"),
            n_partial_anchor = sum(df$anchor_status == "partial_anchor"),
            n_rescued = sum(df$anchor_status == "rescued"),
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
    assignments <- .as_assignment_df(assignments)
    if (!all(c("read_id", "sample_id") %in% names(assignments))) {
        stop("assignments must contain read_id and sample_id columns.")
    }
    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }

    id2sample <- setNames(assignments$sample_id, assignments$read_id)
    un_ids <- character()
    if (isTRUE(include_unassigned) && !is.null(unassigned)) {
        unassigned <- .as_assignment_df(unassigned)
        un_ids <- unique(unassigned$read_id)
    }

    samples <- sort(unique(assignments$sample_id))
    ext <- if (compress) ".fq.gz" else ".fq"
    cons <- list()
    counts <- setNames(integer(length(samples)), samples)
    for (s in samples) {
        path <- file.path(out_dir, paste0(.safe_filename(s), ext))
        cons[[s]] <- if (compress) gzfile(path, open = "wt") else file(path, open = "wt")
    }
    un_con <- NULL
    if (length(un_ids)) {
        path <- file.path(out_dir, paste0("unassigned", ext))
        un_con <- if (compress) gzfile(path, open = "wt") else file(path, open = "wt")
    }
    on.exit({
        for (con in cons) close(con)
        if (!is.null(un_con)) close(un_con)
    }, add = TRUE)

    for (fq in fastq) {
        counts <- .stream_fastq_write(fq, id2sample, cons, counts, un_ids, un_con)
    }
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
    con <- if (grepl("\\.gz$", fastq, ignore.case = TRUE)) {
        gzfile(fastq, open = "rt")
    } else {
        file(fastq, open = "rt")
    }
    on.exit(close(con), add = TRUE)

    repeat {
        h <- readLines(con, n = 1L)
        if (length(h) == 0L) break
        s <- readLines(con, n = 1L)
        p <- readLines(con, n = 1L)
        q <- readLines(con, n = 1L)
        if (length(s) == 0L || length(p) == 0L || length(q) == 0L) break
        rid <- sub("^@", "", strsplit(h, "\\s+")[[1]][1])
        sample <- unname(id2sample[rid])
        if (length(sample) == 1L && !is.na(sample) && !is.null(cons[[sample]])) {
            writeLines(c(h, s, p, q), cons[[sample]])
            counts[[sample]] <- counts[[sample]] + 1L
        } else if (!is.null(un_con) && rid %in% un_ids) {
            writeLines(c(h, s, p, q), un_con)
        }
    }
    counts
}
