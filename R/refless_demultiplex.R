################################################################################
# Demultiplex v2: Index + Primer, edit-distance scoring, 4-tier classification
################################################################################

#' Align a short query to a read-end window and return scores
#' @keywords internal
.align_end_score <- function(query, subject){
    query <- as.character(Biostrings::DNAString(query))
    subject <- as.character(Biostrings::DNAString(subject))
    if(!nchar(query) || !nchar(subject)){
        return(c(score = 0, identity = 0, aln_len = 0, nm = Inf))
    }
    aln <- pwalign::pairwiseAlignment(
        pattern = query,
        subject = subject,
        type = "local",
        gapOpening = 1,
        gapExtension = 1
    )
    aln_pat <- as.character(pwalign::aligned(aln@pattern))
    aln_len <- nchar(gsub("-", "", aln_pat))
    qlen <- nchar(query)
    identity <- if(qlen > 0) aln_len / qlen else 0
    nm <- sum(aln_pat == "-") + sum(strsplit(aln_pat, "")[[1]] != strsplit(substr(query, 1, aln_len), "")[[1]])
    c(score = pwalign::score(aln), identity = identity, aln_len = aln_len, nm = nm)
}

#' Score all samples for one read
#' @keywords internal
.score_read_samples <- function(read_seq,
                                samples,
                                primer_f,
                                primer_r,
                                end_window){
    read_seq <- toupper(as.character(read_seq))
    n <- nchar(read_seq)
    if(n < 20){
        return(NULL)
    }
    w <- min(end_window, n)
    win_start <- substr(read_seq, 1, w)
    win_end <- substr(read_seq, n - w + 1, n)

    n_samp <- nrow(samples)
    score_f <- numeric(n_samp)
    score_r <- numeric(n_samp)
    id_f <- numeric(n_samp)
    id_r <- numeric(n_samp)

    for(i in seq_len(n_samp)){
        q_f <- paste0(samples$f_index[i], primer_f)
        q_r <- paste0(primer_r, samples$r_index[i])
        q_r_rc <- .revcomp(q_r)

        sf1 <- .align_end_score(q_f, win_start)
        sf2 <- .align_end_score(q_r_rc, win_start)
        if(sf1["identity"] >= sf2["identity"]){
            score_f[i] <- sf1["score"]
            id_f[i] <- sf1["identity"]
        } else {
            score_f[i] <- sf2["score"]
            id_f[i] <- sf2["identity"]
        }

        sr1 <- .align_end_score(q_r, win_end)
        sr2 <- .align_end_score(.revcomp(q_f), win_end)
        if(sr1["identity"] >= sr2["identity"]){
            score_r[i] <- sr1["score"]
            id_r[i] <- sr1["identity"]
        } else {
            score_r[i] <- sr2["score"]
            id_r[i] <- sr2["identity"]
        }
    }

    total <- score_f + score_r
    data.frame(
        sample_id = samples$sample_id,
        score_f = score_f,
        score_r = score_r,
        identity_f = id_f,
        identity_r = id_r,
        total_score = total,
        stringsAsFactors = FALSE
    )
}

#' Classify demultiplex assignment from per-read scores
#' @keywords internal
.classify_demultiplex <- function(scores,
                                  t_high = 0.85,
                                  t_low = 0.70,
                                  delta = 0.05,
                                  min_end_id_high = 0.80){
    scores <- scores[order(scores$total_score, decreasing = TRUE), ]
    best <- scores[1, ]
    second <- if(nrow(scores) > 1) scores[2, ] else scores[1, ]

    both_ends <- best$identity_f >= min_end_id_high & best$identity_r >= min_end_id_high
    mean_id <- (best$identity_f + best$identity_r) / 2
    margin <- if(nrow(scores) > 1){
        (best$total_score - second$total_score) / max(abs(best$total_score), 1e-6)
    } else {
        1
    }

    if(best$total_score <= 0 || mean_id < t_low){
        class <- "Unclassified"
        sample_id <- NA_character_
    } else if(margin < delta){
        class <- "Ambiguous"
        sample_id <- best$sample_id
    } else if(both_ends && mean_id >= t_high){
        class <- "High"
        sample_id <- best$sample_id
    } else if(mean_id >= t_low){
        class <- "Low"
        sample_id <- best$sample_id
    } else {
        class <- "Unclassified"
        sample_id <- NA_character_
    }

    list(
        sample_id = sample_id,
        class = class,
        best = best,
        second = second,
        margin = margin
    )
}

#' Demultiplex Nanopore reads using Index + Primer composite search
#'
#' Assigns reads to samples using local alignment of Index+Primer queries to read
#' ends. Supports four confidence classes: High, Low, Ambiguous, Unclassified.
#'
#' @param fastq_fn Path to FASTQ (plain or gz).
#' @param index_list CSV with columns: sample_id, f_index_id, f_index, r_index_id, r_index.
#' @param primer_list CSV with primer IDs ending in _F / _R.
#' @param out_dir Output directory.
#' @param end_window Number of bases at each read end to search.
#' @param t_high Minimum mean identity for High confidence.
#' @param t_low Minimum mean identity for Low confidence.
#' @param delta Minimum relative margin between best and second-best hit.
#' @param min_end_id_high Minimum per-end identity for High (both ends).
#' @param max_reads Optional cap on reads processed (for testing).
#' @param write_per_sample If TRUE, write per-sample FASTQ files.
#' @param n_core Not used yet; reserved for parallel processing.
#' @return Data frame of per-read assignments.
#' @export
doDemultiplex2 <- function(fastq_fn,
                           index_list,
                           primer_list,
                           out_dir,
                           end_window = 200,
                           t_high = 0.85,
                           t_low = 0.70,
                           delta = 0.05,
                           min_end_id_high = 0.80,
                           max_reads = NULL,
                           write_per_sample = TRUE,
                           n_core = 1){
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    samples <- .read_index_list(index_list)
    primers <- .read_primers(primer_list)

    read_fmt <- if(grepl("\\.(fq|fastq)(\\.gz)?$", fastq_fn, ignore.case = TRUE)) "fastq" else "fasta"
    reads <- Biostrings::readDNAStringSet(fastq_fn, format = read_fmt)
    if(!is.null(max_reads) && length(reads) > max_reads){
        reads <- reads[seq_len(max_reads)]
    }

    read_ids <- names(reads)
    if(is.null(read_ids)){
        read_ids <- paste0("read_", seq_along(reads))
    }

    process_one <- function(i){
        scores <- .score_read_samples(
            read_seq = reads[[i]],
            samples = samples,
            primer_f = primers$f,
            primer_r = primers$r,
            end_window = end_window
        )
        hit <- .classify_demultiplex(
            scores,
            t_high = t_high,
            t_low = t_low,
            delta = delta,
            min_end_id_high = min_end_id_high
        )
        data.frame(
            read_id = read_ids[i],
            sample_id = hit$sample_id,
            class = hit$class,
            margin = hit$margin,
            best_sample = hit$best$sample_id,
            score_f = hit$best$score_f,
            score_r = hit$best$score_r,
            identity_f = hit$best$identity_f,
            identity_r = hit$best$identity_r,
            total_score = hit$best$total_score,
            second_sample = hit$second$sample_id,
            stringsAsFactors = FALSE
        )
    }

    if(n_core > 1L && .Platform$OS.type != "windows"){
        results <- parallel::mclapply(seq_along(reads), process_one, mc.cores = n_core)
    } else {
        results <- lapply(seq_along(reads), process_one)
    }
    demult_out <- do.call(rbind, results)
    out_fn <- file.path(out_dir, "demultiplex_assignments.csv")
    write.csv(demult_out, out_fn, row.names = FALSE)

    if(write_per_sample){
        sample_dir <- file.path(out_dir, "per_sample")
        dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
        assigned <- demult_out$class %in% c("High", "Low") & !is.na(demult_out$sample_id)
        for(sid in unique(demult_out$sample_id[assigned])){
            idx <- which(assigned & demult_out$sample_id == sid)
            Biostrings::writeXStringSet(
                reads[idx],
                filepath = file.path(sample_dir, paste0(sid, ".fastq")),
                format = "fastq"
            )
        }
    }

    summary_tab <- as.data.frame(table(demult_out$class))
    names(summary_tab) <- c("class", "count")
    write.csv(summary_tab, file.path(out_dir, "demultiplex_summary.csv"), row.names = FALSE)

    invisible(demult_out)
}
