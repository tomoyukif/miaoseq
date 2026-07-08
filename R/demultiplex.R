################################################################################
# Demultiplexing (edlib: suffix → barcode → optional prefix)
################################################################################

#' Demultiplex reads by dual barcodes (edlib)
#'
#' Classifies multiplexed FASTQ reads into samples using a two-step procedure:
#' find common adapter **suffix** anchors (loose), extract the nearby **barcode**
#' (strict dictionary / edit distance), then optionally verify the outer **prefix**.
#' Primary outputs are assignment tables; sample-specific FASTQ files are written
#' only when `split_reads = TRUE` (via [splitDemultiplexReads()]).
#'
#' @param fastq Character vector of FASTQ paths (`.fq` / `.fastq`, optionally `.gz`).
#' @param demult_dir Output directory for demultiplex results.
#' @param index_list Path to a 5-column headerless CSV:
#'   `index_pair_id`, forward index ID, forward index sequence,
#'   reverse index ID, reverse index sequence.
#' @param sample_list Optional path to a CSV with at least
#'   `index_pair_id` and `sample_name` (headerless or with header).
#'   If omitted, `sample_id` equals `index_pair_id`.
#' @param n_core Number of CPU cores for parallel chunk processing.
#' @param end_window Bases from each read end to search for anchors.
#' @param max_anchor_edit Maximum edit distance for suffix detection.
#' @param max_prefix_edit Maximum edit distance for optional prefix checks.
#' @param max_barcode_edit Maximum edit distance for barcode matching.
#' @param require_prefix If `TRUE`, reject reads without a detectable prefix.
#' @param prefix_rescue If `TRUE`, re-check ambiguous hits with prefix geometry.
#' @param require_unique_pair If `TRUE`, reject non-unique valid F/R pairs.
#' @param allow_revcomp If `TRUE`, also search reverse-complement windows.
#' @param split_reads If `TRUE`, call [splitDemultiplexReads()] after assignment.
#' @param compress Passed to [splitDemultiplexReads()] when `split_reads = TRUE`.
#' @param chunk_size Number of reads per processing chunk.
#'
#' @return A list with `assignments`, `summary`, `unassigned`, and `demult_dir`.
#'
#' @export
#' @importFrom Biostrings DNAStringSet reverseComplement readDNAStringSet
#' @importFrom Biostrings writeXStringSet QualityScaledDNAStringSet
#' @importFrom BiocGenerics width
#' @importFrom edlibR align
#' @importFrom parallel mclapply
#' @importFrom stats setNames
#' @importFrom utils read.csv write.table
doDemultiplex <- function(fastq,
                          demult_dir,
                          index_list,
                          sample_list = NULL,
                          n_core = 1,
                          end_window = 60,
                          max_anchor_edit = 10,
                          max_prefix_edit = 5,
                          max_barcode_edit = 2,
                          require_prefix = FALSE,
                          prefix_rescue = TRUE,
                          require_unique_pair = TRUE,
                          allow_revcomp = TRUE,
                          split_reads = FALSE,
                          compress = TRUE,
                          chunk_size = 5000) {
    if (!dir.exists(demult_dir)) {
        dir.create(demult_dir, recursive = TRUE, showWarnings = FALSE)
    }
    fastq <- as.character(fastq)
    if (length(fastq) < 1 || any(!file.exists(fastq))) {
        stop("All fastq paths must exist.")
    }

    layout <- .parse_index_layout(index_list)
    f_dict <- .build_barcode_dict(layout$f_barcodes, max_barcode_edit)
    r_dict <- .build_barcode_dict(layout$r_barcodes, max_barcode_edit)
    design <- .validate_barcode_design(layout, max_barcode_edit)
    sample_map <- .load_sample_map(sample_list, layout$sample_map)

    .write_index_layout(layout, file.path(demult_dir, "index_layout.tsv"))
    write.table(rbind(cbind(side = "F", f_dict$conflicts),
                      cbind(side = "R", r_dict$conflicts)),
                file.path(demult_dir, "barcode_conflicts.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(design, file.path(demult_dir, "design_check.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)

    ctx <- list(
        layout = layout,
        f_dict = f_dict,
        r_dict = r_dict,
        sample_map = sample_map,
        end_window = as.integer(end_window),
        max_anchor_edit = as.integer(max_anchor_edit),
        max_prefix_edit = as.integer(max_prefix_edit),
        max_barcode_edit = as.integer(max_barcode_edit),
        require_prefix = isTRUE(require_prefix),
        prefix_rescue = isTRUE(prefix_rescue),
        require_unique_pair = isTRUE(require_unique_pair),
        allow_revcomp = isTRUE(allow_revcomp)
    )

    results <- list()
    for (fq in fastq) {
        message("Demultiplexing: ", fq)
        results[[fq]] <- .demultiplex_fastq_file(
            fq, ctx, n_core = n_core, chunk_size = chunk_size
        )
    }
    assignments <- do.call(rbind, lapply(results, `[[`, "assignments"))
    unassigned <- do.call(rbind, lapply(results, `[[`, "unassigned"))
    if (is.null(assignments)) {
        assignments <- .empty_assignments()
    }
    if (is.null(unassigned)) {
        unassigned <- .empty_unassigned()
    }
    rownames(assignments) <- NULL
    rownames(unassigned) <- NULL

    summary_df <- .summarize_assignments(assignments)
    .write_demultiplex_tables(demult_dir, assignments, summary_df, unassigned)

    if (isTRUE(split_reads) && nrow(assignments) > 0) {
        splitDemultiplexReads(
            fastq = fastq,
            assignments = assignments,
            out_dir = file.path(demult_dir, "by_sample"),
            compress = compress,
            include_unassigned = FALSE,
            n_core = n_core
        )
    }

    list(
        assignments = assignments,
        summary = summary_df,
        unassigned = unassigned,
        demult_dir = demult_dir
    )
}

#' Split a FASTQ by demultiplex assignments
#'
#' Writes one FASTQ per `sample_id` using `assignments` from [doDemultiplex()].
#' Can be run after `doDemultiplex(..., split_reads = FALSE)`.
#'
#' @param fastq Character vector of source FASTQ paths.
#' @param assignments `assignments.tsv` path or a data.frame with
#'   `read_id` and `sample_id` columns (and ideally `source_file`).
#' @param out_dir Output directory for per-sample FASTQ files.
#' @param compress If `TRUE`, write `.fq.gz`.
#' @param include_unassigned If `TRUE`, also write `unassigned.fq(.gz)`.
#' @param unassigned Optional unassigned table / TSV path.
#' @param n_core Unused currently (reserved); splitting is single-pass I/O.
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
# Phase A: layout / dictionary
################################################################################

#' @keywords internal
.parse_index_layout <- function(index_list) {
    index_df <- utils::read.csv(index_list, header = FALSE, stringsAsFactors = FALSE)
    if (ncol(index_df) < 5) {
        stop("index_list must have 5 columns.")
    }
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
.build_barcode_dict <- function(barcodes, max_barcode_edit = 2L) {
    ids <- names(barcodes)
    if (is.null(ids) || anyNA(ids) || any(ids == "")) {
        stop("barcodes must be a named character vector.")
    }
    barcodes <- setNames(as.character(unname(barcodes)), ids)
    max_barcode_edit <- as.integer(max_barcode_edit)

    # mutant_seq -> character vector of parent ids
    owners <- new.env(hash = TRUE, parent = emptyenv())

    add_owner <- function(seq, id) {
        if (exists(seq, envir = owners, inherits = FALSE)) {
            owners[[seq]] <- c(owners[[seq]], id)
        } else {
            owners[[seq]] <- id
        }
    }

    for (i in seq_along(barcodes)) {
        parent <- barcodes[[i]]
        id <- ids[[i]]
        add_owner(parent, id)
        if (max_barcode_edit < 1L) next
        level1 <- .mutate_barcode(parent)
        for (m in level1) add_owner(m, id)
        if (max_barcode_edit >= 2L) {
            for (m1 in level1) {
                for (m2 in .mutate_barcode(m1)) {
                    if (!identical(m2, parent)) add_owner(m2, id)
                }
            }
        }
    }

    keys <- ls(owners, all.names = TRUE)
    mutant_map <- vector("list", length(keys))
    names(mutant_map) <- keys
    conflict_mut <- character()
    conflict_par <- character()
    keep <- logical(length(keys))
    for (i in seq_along(keys)) {
        k <- keys[[i]]
        u <- unique(owners[[k]])
        if (length(u) == 1L) {
            mutant_map[[i]] <- u
            keep[i] <- TRUE
        } else {
            conflict_mut <- c(conflict_mut, k)
            conflict_par <- c(conflict_par, paste(u, collapse = ","))
            keep[i] <- FALSE
        }
    }
    mutant_map <- mutant_map[keep]

    list(
        barcodes = barcodes,
        mutant_map = mutant_map,
        conflicts = data.frame(
            mutant = conflict_mut,
            parents = conflict_par,
            stringsAsFactors = FALSE
        )
    )
}

#' @keywords internal
.mutate_barcode <- function(seq) {
    bases <- c("A", "C", "G", "T")
    chars <- strsplit(seq, "", fixed = TRUE)[[1]]
    n <- length(chars)
    out <- new.env(hash = TRUE, parent = emptyenv())
    put <- function(x) out[[x]] <- TRUE

    for (i in seq_len(n)) {
        orig <- chars[i]
        for (b in bases) {
            if (b != orig) {
                chars[i] <- b
                put(paste(chars, collapse = ""))
            }
        }
        chars[i] <- orig
    }
    if (n > 1L) {
        for (i in seq_len(n)) {
            put(paste(chars[-i], collapse = ""))
        }
    }
    for (i in 0:n) {
        left <- if (i == 0L) "" else paste(chars[seq_len(i)], collapse = "")
        right <- if (i == n) "" else paste(chars[(i + 1L):n], collapse = "")
        for (b in bases) put(paste0(left, b, right))
    }
    ls(out, all.names = TRUE)
}

#' @keywords internal
.validate_barcode_design <- function(layout, max_barcode_edit) {
    f_min <- .min_pairwise_edit(layout$f_barcodes)
    r_min <- .min_pairwise_edit(layout$r_barcodes)
    need <- 2L * as.integer(max_barcode_edit) + 1L
    data.frame(
        side = c("F", "R"),
        n_barcodes = c(length(layout$f_barcodes), length(layout$r_barcodes)),
        min_edit_distance = c(f_min, r_min),
        recommended_min = c(need, need),
        ok = c(f_min >= need, r_min >= need),
        note = c(
            if (f_min < need) "Use collision removal; consider lowering max_barcode_edit" else "OK",
            if (r_min < need) "Use collision removal; consider lowering max_barcode_edit" else "OK"
        ),
        stringsAsFactors = FALSE
    )
}

#' @keywords internal
.min_pairwise_edit <- function(seqs) {
    seqs <- unique(as.character(seqs))
    if (length(seqs) < 2) return(NA_integer_)
    best <- Inf
    for (i in seq_len(length(seqs) - 1L)) {
        for (j in (i + 1L):length(seqs)) {
            d <- edlibR::align(seqs[i], seqs[j], mode = "NW", task = "distance")$editDistance
            if (!is.na(d) && d < best) best <- d
        }
    }
    as.integer(best)
}

#' @keywords internal
.load_sample_map <- function(sample_list, sample_map) {
    sample_map$sample_id <- sample_map$index_pair_id
    if (is.null(sample_list)) return(sample_map)
    sl <- utils::read.csv(sample_list, header = FALSE, stringsAsFactors = FALSE)
    if (ncol(sl) < 2) stop("sample_list needs at least 2 columns: index_pair_id, sample_name.")
    # If first row looks like a header, drop it
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
# Phase B: per-read scoring
################################################################################

#' @keywords internal
.demultiplex_fastq_file <- function(fastq, ctx, n_core = 1, chunk_size = 5000) {
    reads <- .read_fastq_seqs(fastq)
    n <- length(reads)
    if (n < 1) {
        return(list(assignments = .empty_assignments(),
                    unassigned = .empty_unassigned()))
    }
    ids <- names(reads)
    seqs <- as.character(reads)
    starts <- seq(1L, n, by = chunk_size)
    chunks <- lapply(starts, function(s) {
        e <- min(s + chunk_size - 1L, n)
        list(ids = ids[s:e], seqs = seqs[s:e])
    })

    worker <- function(chunk) {
        .demultiplex_chunk(chunk$ids, chunk$seqs, source_file = fastq, ctx = ctx)
    }

    if (n_core > 1 && length(chunks) > 1) {
        parts <- parallel::mclapply(chunks, worker, mc.cores = n_core)
    } else {
        parts <- lapply(chunks, worker)
    }

    assignments <- do.call(rbind, lapply(parts, `[[`, "assignments"))
    unassigned <- do.call(rbind, lapply(parts, `[[`, "unassigned"))
    list(assignments = assignments, unassigned = unassigned)
}

#' @keywords internal
.read_fastq_seqs <- function(fastq) {
    # Quality is ignored for demux; keep sequence + IDs.
    ext <- tolower(fastq)
    if (grepl("\\.gz$", ext)) {
        con <- gzfile(fastq, open = "rt")
        on.exit(close(con), add = TRUE)
        # Biostrings handles gz paths directly
        ss <- Biostrings::readDNAStringSet(fastq, format = "fastq")
    } else if (grepl("\\.(fq|fastq)$", ext)) {
        ss <- Biostrings::readDNAStringSet(fastq, format = "fastq")
    } else {
        # allow FASTA for testing
        ss <- Biostrings::readDNAStringSet(fastq)
    }
    names(ss) <- vapply(strsplit(names(ss), "\\s+"), `[[`, character(1), 1L)
    ss
}

#' @keywords internal
.demultiplex_chunk <- function(ids, seqs, source_file, ctx) {
    assign_rows <- list()
    un_rows <- list()
    for (i in seq_along(ids)) {
        res <- .score_read_indices(ids[[i]], seqs[[i]], source_file, ctx)
        if (identical(res$status, "assigned")) {
            assign_rows[[length(assign_rows) + 1L]] <- res$row
        } else {
            un_rows[[length(un_rows) + 1L]] <- data.frame(
                read_id = ids[[i]],
                reason = res$reason,
                source_file = source_file,
                stringsAsFactors = FALSE
            )
        }
    }
    list(
        assignments = if (length(assign_rows)) do.call(rbind, assign_rows) else .empty_assignments(),
        unassigned = if (length(un_rows)) do.call(rbind, un_rows) else .empty_unassigned()
    )
}

#' @keywords internal
.score_read_indices <- function(read_id, seq, source_file, ctx) {
    seq <- toupper(gsub("[^ACGT]", "N", seq))
    n <- nchar(seq)
    w <- min(ctx$end_window, n)
    front <- substr(seq, 1, w)
    rear <- substr(seq, n - w + 1, n)

    f_hit <- .match_end(
        window = front,
        suffix = ctx$layout$f_suffix,
        prefix = ctx$layout$f_prefix,
        barcode_len = ctx$layout$f_barcode_len,
        dict = ctx$f_dict,
        ctx = ctx,
        barcode_before_suffix = TRUE
    )
    r_hit <- .match_end(
        window = rear,
        suffix = ctx$layout$r_suffix,
        prefix = ctx$layout$r_prefix,
        barcode_len = ctx$layout$r_barcode_len,
        dict = ctx$r_dict,
        ctx = ctx,
        barcode_before_suffix = FALSE
    )

    # Orientation flip: F may be on rear / R on front (adapter strand order unchanged)
    if ((is.null(f_hit) || is.null(r_hit)) && ctx$allow_revcomp) {
        f_hit2 <- .match_end(
            window = rear,
            suffix = ctx$layout$f_suffix,
            prefix = ctx$layout$f_prefix,
            barcode_len = ctx$layout$f_barcode_len,
            dict = ctx$f_dict,
            ctx = ctx,
            barcode_before_suffix = TRUE
        )
        r_hit2 <- .match_end(
            window = front,
            suffix = ctx$layout$r_suffix,
            prefix = ctx$layout$r_prefix,
            barcode_len = ctx$layout$r_barcode_len,
            dict = ctx$r_dict,
            ctx = ctx,
            barcode_before_suffix = FALSE
        )
        if (!is.null(f_hit2) && !is.null(r_hit2)) {
            f_hit <- f_hit2
            r_hit <- r_hit2
        }
    }

    if (is.null(f_hit) && is.null(r_hit)) {
        return(list(status = "unassigned", reason = "no_suffix"))
    }
    if (is.null(f_hit) || is.null(r_hit)) {
        return(list(status = "unassigned", reason = "barcode_fail"))
    }
    if (identical(f_hit$reason, "ambiguous") || identical(r_hit$reason, "ambiguous")) {
        if (ctx$prefix_rescue) {
            rescued <- .prefix_rescue_pair(f_hit, r_hit, ctx)
            if (is.null(rescued)) {
                return(list(status = "unassigned", reason = "rescue_fail"))
            }
            f_hit <- rescued$f
            r_hit <- rescued$r
            f_hit$rescued <- TRUE
            r_hit$rescued <- TRUE
        } else {
            return(list(status = "unassigned", reason = "ambiguous_pair"))
        }
    }
    if (is.null(f_hit$id) || is.null(r_hit$id)) {
        return(list(status = "unassigned", reason = "barcode_fail"))
    }

    if (ctx$require_prefix &&
        (is.na(f_hit$prefix_edit) || is.na(r_hit$prefix_edit))) {
        return(list(status = "unassigned", reason = "no_prefix"))
    }

    pair_key <- paste(f_hit$id, r_hit$id, sep = "\t")
    sm <- ctx$sample_map
    hits <- which(sm$pair_key == pair_key)
    if (length(hits) == 0) {
        return(list(status = "unassigned", reason = "invalid_pair"))
    }
    if (ctx$require_unique_pair && length(hits) > 1) {
        return(list(status = "unassigned", reason = "ambiguous_pair"))
    }
    hit <- hits[[1]]

    match_class <- if (f_hit$barcode_edit == 0L && r_hit$barcode_edit == 0L) {
        "complete_match"
    } else {
        "fuzzy_match"
    }
    rescued <- isTRUE(f_hit$rescued) || isTRUE(r_hit$rescued)
    both_prefix <- !is.na(f_hit$prefix_edit) && !is.na(r_hit$prefix_edit)
    anchor_status <- if (rescued) {
        "rescued"
    } else if (both_prefix) {
        "high_confidence"
    } else {
        "partial_anchor"
    }

    row <- data.frame(
        read_id = read_id,
        index_pair_id = sm$index_pair_id[hit],
        f_index_id = f_hit$id,
        r_index_id = r_hit$id,
        sample_id = sm$sample_id[hit],
        source_file = source_file,
        barcode_edit_f = f_hit$barcode_edit,
        barcode_edit_r = r_hit$barcode_edit,
        match_class = match_class,
        anchor_status = anchor_status,
        prefix_edit_f = f_hit$prefix_edit,
        prefix_edit_r = r_hit$prefix_edit,
        stringsAsFactors = FALSE
    )
    list(status = "assigned", row = row)
}

#' @keywords internal
.find_suffix_anchor <- function(window, suffix, max_anchor_edit, allow_revcomp = TRUE) {
    window <- toupper(window)
    cands <- list()
    consider <- function(seq, flipped) {
        aln <- edlibR::align(suffix, seq, mode = "HW", task = "locations",
                             k = max_anchor_edit)
        if (is.null(aln$editDistance) || is.na(aln$editDistance) ||
            aln$editDistance < 0 || aln$editDistance > max_anchor_edit) {
            return(invisible(NULL))
        }
        locs <- aln$locations
        if (is.null(locs) || length(locs) < 1) return(invisible(NULL))
        for (loc in locs) {
            cands[[length(cands) + 1L]] <<- list(
                seq = seq,
                flipped = flipped,
                edit = as.integer(aln$editDistance),
                start = as.integer(loc[1] + 1L),
                end = as.integer(loc[2] + 1L)
            )
        }
        invisible(NULL)
    }
    consider(window, FALSE)
    if (allow_revcomp) {
        consider(.revcomp_char(window), TRUE)
    }
    if (!length(cands)) return(NULL)
    # Prefer lower edit distance, then enough room for a barcode before suffix.
    edits <- vapply(cands, `[[`, integer(1), "edit")
    cands <- cands[order(edits)]
    cands
}

#' @keywords internal
.match_end <- function(window, suffix, prefix, barcode_len, dict, ctx,
                       barcode_before_suffix = TRUE) {
    anchors <- .find_suffix_anchor(window, suffix, ctx$max_anchor_edit, ctx$allow_revcomp)
    if (is.null(anchors)) return(NULL)

    best <- NULL
    max_shift <- min(4L, as.integer(ctx$max_anchor_edit))
    for (anchor in anchors) {
        for (off in c(0L, unlist(lapply(seq_len(max_shift), function(x) c(-x, x))))) {
            bc <- .extract_barcode(anchor$seq, anchor$start, anchor$end,
                                   barcode_len, barcode_before_suffix, offset = off)
            if (is.null(bc)) next
            matched <- .match_barcode(bc, dict, ctx$max_barcode_edit)
            barcode_start <- if (barcode_before_suffix) {
                anchor$start - barcode_len + off
            } else {
                anchor$end + 1L + off
            }
            barcode_end <- barcode_start + barcode_len - 1L
            prefix_edit <- .check_prefix_optional(
                seq = anchor$seq,
                prefix = prefix,
                barcode_start = barcode_start,
                barcode_end = barcode_end,
                max_prefix_edit = ctx$max_prefix_edit,
                barcode_before_suffix = barcode_before_suffix
            )
            ambiguous <- isTRUE(matched$ambiguous)
            cand <- list(
                id = matched$id,
                barcode_edit = matched$edit,
                prefix_edit = prefix_edit,
                reason = if (is.null(matched$id) && !ambiguous) "barcode_fail" else
                    if (ambiguous) "ambiguous" else NA_character_,
                ambiguous = ambiguous,
                candidates = matched$candidates,
                barcode_seq = bc,
                anchor = anchor,
                rescued = FALSE
            )
            score <- c(
                as.integer(!is.null(cand$id) || ambiguous),
                -if (is.null(cand$barcode_edit) || is.na(cand$barcode_edit)) 99L else cand$barcode_edit,
                -abs(as.integer(off)),
                -anchor$edit
            )
            if (is.null(best) ||
                score[1] > best$.score[1] ||
                (score[1] == best$.score[1] && score[2] > best$.score[2]) ||
                (score[1] == best$.score[1] && score[2] == best$.score[2] &&
                 score[3] > best$.score[3]) ||
                (score[1] == best$.score[1] && score[2] == best$.score[2] &&
                 score[3] == best$.score[3] && score[4] > best$.score[4])) {
                best <- cand
                best$.score <- score
            }
            if (!is.null(cand$id) && !ambiguous &&
                identical(cand$barcode_edit, 0L) && identical(off, 0L)) {
                best$.score <- NULL
                return(best)
            }
        }
    }
    if (is.null(best)) return(NULL)
    best$.score <- NULL
    best
}

#' @keywords internal
.revcomp_char <- function(seq) {
    as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(seq)))
}

#' @keywords internal
.extract_barcode <- function(seq, suffix_start, suffix_end, barcode_len,
                             barcode_before_suffix = TRUE, offset = 0L) {
    if (barcode_before_suffix) {
        to <- suffix_start - 1L + offset
        from <- to - barcode_len + 1L
    } else {
        from <- suffix_end + 1L + offset
        to <- from + barcode_len - 1L
    }
    if (from < 1 || to > nchar(seq) || from > to) return(NULL)
    substr(seq, from, to)
}

#' @keywords internal
.match_barcode <- function(barcode_seq, dict, max_barcode_edit) {
    barcode_seq <- toupper(barcode_seq)
    if (!is.null(dict$mutant_map[[barcode_seq]])) {
        id <- dict$mutant_map[[barcode_seq]]
        edit <- edlibR::align(barcode_seq, dict$barcodes[[id]],
                              mode = "NW", task = "distance")$editDistance
        return(list(id = id, edit = as.integer(edit),
                    ambiguous = FALSE, candidates = id))
    }
    # fallback: NW against all barcodes
    edits <- vapply(dict$barcodes, function(bc) {
        edlibR::align(barcode_seq, bc, mode = "NW",
                      task = "distance", k = max_barcode_edit)$editDistance
    }, integer(1))
    edits[is.na(edits) | edits < 0] <- max_barcode_edit + 1L
    ok <- which(edits <= max_barcode_edit)
    if (length(ok) == 0) {
        return(list(id = NULL, edit = NA_integer_, ambiguous = FALSE,
                    candidates = character()))
    }
    min_e <- min(edits[ok])
    cand <- names(dict$barcodes)[ok][edits[ok] == min_e]
    if (length(cand) > 1) {
        return(list(id = NULL, edit = as.integer(min_e),
                    ambiguous = TRUE, candidates = cand))
    }
    list(id = cand[[1]], edit = as.integer(min_e),
         ambiguous = FALSE, candidates = cand)
}

#' @keywords internal
.check_prefix_optional <- function(seq, prefix, barcode_start, barcode_end,
                                   max_prefix_edit, barcode_before_suffix) {
    # Search a small window outside the barcode for prefix.
    if (barcode_before_suffix) {
        # prefix should be immediately before barcode
        from <- max(1L, barcode_start - nchar(prefix) - 4L)
        to <- max(1L, barcode_start - 1L)
    } else {
        # prefix after barcode toward outer end
        from <- min(nchar(seq), barcode_end + 1L)
        to <- min(nchar(seq), barcode_end + nchar(prefix) + 4L)
    }
    if (from > to) return(NA_integer_)
    region <- substr(seq, from, to)
    if (nchar(region) < 1) return(NA_integer_)
    aln <- edlibR::align(prefix, region, mode = "HW", task = "distance",
                         k = max_prefix_edit)
    if (is.null(aln$editDistance) || is.na(aln$editDistance) ||
        aln$editDistance < 0 || aln$editDistance > max_prefix_edit) {
        return(NA_integer_)
    }
    as.integer(aln$editDistance)
}

#' @keywords internal
.prefix_rescue_pair <- function(f_hit, r_hit, ctx) {
    # Keep candidate IDs that have detectable prefix when ambiguous.
    resolve_one <- function(hit, prefix, barcode_len, before) {
        if (!isTRUE(hit$ambiguous) && !is.null(hit$id)) return(hit)
        cands <- hit$candidates
        if (length(cands) < 1) return(NULL)
        keep <- cands
        # Prefer candidates whose expected barcode equals extracted (already true)
        # and whose prefix is present.
        if (is.na(hit$prefix_edit)) return(NULL)
        if (length(keep) == 1) {
            hit$id <- keep[[1]]
            hit$ambiguous <- FALSE
            hit$reason <- NA_character_
            return(hit)
        }
        NULL
    }
    f2 <- resolve_one(f_hit, ctx$layout$f_prefix, ctx$layout$f_barcode_len, TRUE)
    r2 <- resolve_one(r_hit, ctx$layout$r_prefix, ctx$layout$r_barcode_len, FALSE)
    if (is.null(f2) || is.null(r2) || is.null(f2$id) || is.null(r2$id)) return(NULL)
    # Among ambiguous candidate pairs, keep those present in sample_map
    f_cands <- if (length(f_hit$candidates)) f_hit$candidates else f2$id
    r_cands <- if (length(r_hit$candidates)) r_hit$candidates else r2$id
    keys <- as.vector(outer(f_cands, r_cands, paste, sep = "\t"))
    hit_i <- which(ctx$sample_map$pair_key %in% keys)
    if (length(hit_i) != 1L) return(NULL)
    sm <- ctx$sample_map[hit_i, ]
    f2$id <- sm$f_index_id
    r2$id <- sm$r_index_id
    f2$ambiguous <- FALSE
    r2$ambiguous <- FALSE
    list(f = f2, r = r2)
}

################################################################################
# Phase C: outputs
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
# Phase D: split FASTQ
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
    un_count <- 0L
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
