################################################################################
# Edit-calling module (consensus-based, no check_window)
################################################################################

#' Map short sequences to reference with minimap2 and parse PAF hits
#' @keywords internal
.map_short_seqs_minimap <- function(query_fa, ref_fa, out_dir, sif_path = NULL){
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    local_ref <- file.path(out_dir, "ref.fa")
    if(normalizePath(ref_fa) != normalizePath(local_ref)){
        file.copy(ref_fa, local_ref, overwrite = TRUE)
    }
    local_q <- file.path(out_dir, basename(query_fa))
    if(normalizePath(query_fa) != normalizePath(local_q)){
        file.copy(query_fa, local_q, overwrite = TRUE)
    }
    paf <- file.path(out_dir, paste0(sub("\\.[^.]+$", "", basename(query_fa)), ".paf"))
    cmd <- sprintf(
        "minimap2 -x map-ont %s %s > %s",
        shQuote(basename(local_ref)),
        shQuote(basename(local_q)),
        shQuote(basename(paf))
    )
    status <- .refless_run_shell(cmd, work_dir = out_dir, sif_path = sif_path)
    if(!is.null(status) && status != 0){
        stop("minimap2 mapping failed for ", query_fa)
    }
    if(!file.exists(paf) || !nzchar(readLines(paf, n = 1, warn = FALSE))){
        return(data.frame())
    }
    lines <- readLines(paf)
    parsed <- lapply(strsplit(lines, "\t"), function(x){
        if(length(x) < 12){
            return(NULL)
        }
        qname <- x[1]
        tname <- x[6]
        tstart <- as.integer(x[8])
        tend <- as.integer(x[9])
        nmatch <- as.integer(x[10])
        qlen <- as.integer(x[2])
        data.frame(
            qseqid = qname,
            sseqid = tname,
            sstart = tstart,
            send = tend,
            length = abs(tend - tstart) + 1L,
            qlen = qlen,
            pident = 100 * nmatch / qlen,
            stringsAsFactors = FALSE
        )
    })
    out <- do.call(rbind, parsed[!vapply(parsed, is.null, logical(1))])
    if(nrow(out) == 0){
        return(out)
    }
    out[order(out$qseqid, -out$pident), ]
}

#' Define reference targets for consensus-based edit calling
#'
#' Maps primers to the reference genome, defines amplicon coordinates, and
#' extracts wild-type sequences around gRNA cut sites.
#'
#' @param genome_fn Reference genome FASTA.
#' @param primer_list Primer CSV (_F / _R suffixes).
#' @param pam_list gRNA/PAM CSV: gene, chromosome number, cut position, optional guide_id.
#' @param out_dir Output directory (ref_targets/).
#' @param sif_path Optional Apptainer SIF path.
#' @return List with paths to targets.bed, intact_sequences.fa, target_metadata.csv.
#' @export
prepRefTargets <- function(genome_fn,
                           primer_list,
                           pam_list,
                           out_dir,
                           sif_path = NULL){
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    primers <- read.csv(primer_list, header = FALSE, stringsAsFactors = FALSE)
    f_primers <- Biostrings::DNAStringSet(primers$V2[grep("_F$", primers$V1)])
    names(f_primers) <- primers$V1[grep("_F$", primers$V1)]
    r_primers <- Biostrings::DNAStringSet(primers$V2[grep("_R$", primers$V1)])
    names(r_primers) <- primers$V1[grep("_R$", primers$V1)]

    map_dir <- file.path(out_dir, "primer_map")
    dir.create(map_dir, recursive = TRUE, showWarnings = FALSE)
    f_primer_fn <- file.path(map_dir, "f_primer.fa")
    r_primer_fn <- file.path(map_dir, "r_primer.fa")
    Biostrings::writeXStringSet(f_primers, f_primer_fn)
    Biostrings::writeXStringSet(r_primers, r_primer_fn)
    ref_fn <- file.path(map_dir, "ref_genome.fa")
    file.copy(genome_fn, ref_fn, overwrite = TRUE)

    f_out <- .map_short_seqs_minimap(f_primer_fn, ref_fn, file.path(map_dir, "f"), sif_path = sif_path)
    r_out <- .map_short_seqs_minimap(r_primer_fn, ref_fn, file.path(map_dir, "r"), sif_path = sif_path)
    if(nrow(f_out) == 0 || nrow(r_out) == 0){
        stop("Failed to map primers to reference genome")
    }

    f_out <- f_out[!duplicated(f_out$qseqid), ]
    r_out <- r_out[!duplicated(r_out$qseqid), ]
    f_hit <- f_out$length == Biostrings::width(f_primers[f_out$qseqid]) & f_out$pident >= 90
    r_hit <- r_out$length == Biostrings::width(r_primers[r_out$qseqid]) & r_out$pident >= 90
    f_out <- f_out[f_hit, , drop = FALSE]
    r_out <- r_out[r_hit, , drop = FALSE]

    f_out$qseqid <- sub("_F$", "", f_out$qseqid)
    r_out$qseqid <- sub("_R$", "", r_out$qseqid)
    amplicon_df <- dplyr::full_join(f_out, r_out, by = c("qseqid", "sseqid"))
    amplicon_df <- amplicon_df[, c("qseqid", "sseqid", "sstart.x", "send.x", "sstart.y", "send.y")]
    amplicon_df$start <- apply(amplicon_df[, -(1:2)], 1, min)
    amplicon_df$end <- apply(amplicon_df[, -(1:2)], 1, max)

    edit_site <- read.csv(pam_list, header = FALSE, stringsAsFactors = FALSE)
    chr_names <- paste0("chr", sprintf("%02d", edit_site$V2))
    cut_sites <- edit_site$V3
    target_genes <- edit_site$V1
    guide_ids <- if(ncol(edit_site) >= 4) edit_site$V4 else rep(NA_character_, nrow(edit_site))

    ref_genome <- Biostrings::readDNAStringSet(ref_fn)
    amplicon_gr <- GenomicRanges::GRanges(
        seqnames = amplicon_df$sseqid,
        ranges = IRanges::IRanges(start = amplicon_df$start, end = amplicon_df$end),
        id = amplicon_df$qseqid
    )
    amplicon_seq <- ref_genome[amplicon_gr]
    names(amplicon_seq) <- amplicon_gr$id

  # Intact reference: full amplicon per target gene (no check_window trimming)
    intact_list <- list()
    meta_rows <- list()
    for(i in seq_len(nrow(edit_site))){
        gene <- target_genes[i]
        chr <- chr_names[i]
        cut <- cut_sites[i]
        amp_i <- amplicon_df[amplicon_df$qseqid == gene, , drop = FALSE]
        if(nrow(amp_i) == 0){
            next
        }
        amp_chr <- amp_i$sseqid[1]
        amp_start <- amp_i$start[1]
        amp_end <- amp_i$end[1]
        target_name <- gene
        if(!is.na(guide_ids[i]) && nzchar(guide_ids[i])){
            target_name <- paste0(gene, "_", guide_ids[i])
        }
        gr <- GenomicRanges::GRanges(
            seqnames = amp_chr,
            ranges = IRanges::IRanges(start = amp_start, end = amp_end),
            target = target_name
        )
        intact_list[[target_name]] <- as.character(ref_genome[gr][[1]])
        meta_rows[[target_name]] <- data.frame(
            target_gene = target_name,
            gene = gene,
            guide_id = guide_ids[i],
            chromosome = chr,
            cut_site = cut,
            amplicon_chr = amp_chr,
            amplicon_start = amp_start,
            amplicon_end = amp_end,
            stringsAsFactors = FALSE
        )
    }
    intact_seq <- Biostrings::DNAStringSet(intact_list)
    meta <- do.call(rbind, meta_rows)

    bed_fn <- file.path(out_dir, "targets.bed")
    bed <- data.frame(
        chrom = meta$amplicon_chr,
        start = meta$amplicon_start - 1L,
        end = meta$amplicon_end,
        name = meta$target_gene,
        stringsAsFactors = FALSE
    )
    write.table(bed, bed_fn, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

    intact_fn <- file.path(out_dir, "intact_sequences.fa")
    Biostrings::writeXStringSet(intact_seq, intact_fn)
    meta_fn <- file.path(out_dir, "target_metadata.csv")
    write.csv(meta, meta_fn, row.names = FALSE)

    invisible(list(
        targets_bed = bed_fn,
        intact_sequences = intact_fn,
        target_metadata = meta_fn,
        amplicon_fasta = intact_fn,
        metadata = meta,
        intact_seq = intact_seq
    ))
}

#' Extract edit-site window sequence from a global alignment
#' @keywords internal
.extract_aligned_window <- function(aln, window_start, window_end){
    pat <- as.character(pwalign::aligned(aln@pattern))
    subj <- as.character(pwalign::aligned(aln@subject))
    ref_pos <- 0
    q_pos <- 0
    ref_chars <- character()
    q_chars <- character()
    for(i in seq_len(nchar(pat))){
        p <- substr(pat, i, i)
        s <- substr(subj, i, i)
        if(s != "-"){
            ref_pos <- ref_pos + 1
        }
        if(p != "-"){
            q_pos <- q_pos + 1
        }
        if(ref_pos >= window_start && ref_pos <= window_end){
            ref_chars <- c(ref_chars, if(s == "-") "" else s)
            q_chars <- c(q_chars, if(p == "-") "-" else p)
        }
        if(ref_pos > window_end){
            break
        }
    }
    list(
        ref_window = paste(ref_chars, collapse = ""),
        query_window = paste(q_chars, collapse = "")
    )
}

#' Align consensus sequences to reference targets
#'
#' @param consensus_df Consensus summary from doConsensus (must include consensus_file).
#' @param ref_targets Output list from prepRefTargets.
#' @param out_dir Output directory.
#' @param cut_site_margin bp around cut site for edit evaluation (default 3).
#' @return Data frame of alignment metrics per sample/cluster/target.
#' @export
alignConsensusToRef <- function(consensus_df,
                                ref_targets,
                                out_dir,
                                cut_site_margin = 3,
                                sample_id = NULL){
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    meta <- ref_targets$metadata
    intact <- ref_targets$intact_seq
    if(is.null(sample_id)){
        sample_id <- basename(dirname(dirname(out_dir)))
    }

    rows <- list()
    for(i in seq_len(nrow(consensus_df))){
        cid <- consensus_df$cluster_id[i]
        cons <- as.character(Biostrings::readDNAStringSet(consensus_df$consensus_file[i])[[1]])
        for(j in seq_len(nrow(meta))){
            target <- meta$target_gene[j]
            ref_seq <- as.character(intact[[target]])
            aln <- pwalign::pairwiseAlignment(
                pattern = cons,
                subject = ref_seq,
                type = "global",
                gapOpening = 1,
                gapExtension = 1
            )
            pat_aln <- as.character(pwalign::aligned(aln@pattern))
            sub_aln <- as.character(pwalign::aligned(aln@subject))
            matches <- sum(strsplit(pat_aln, "")[[1]] == strsplit(sub_aln, "")[[1]] & pat_aln != "-")
            aln_len <- sum(pat_aln != "-")
            identity <- if(aln_len > 0) matches / aln_len else 0

            cut_in_amp <- meta$cut_site[j] - meta$amplicon_start[j] + 1L
            w_start <- max(1L, cut_in_amp - cut_site_margin)
            w_end <- min(nchar(ref_seq), cut_in_amp + cut_site_margin)
            win <- .extract_aligned_window(aln, w_start, w_end)

            rows[[length(rows) + 1]] <- data.frame(
                sample_id = sample_id,
                cluster_id = cid,
                target_gene = target,
                alignment_score = pwalign::score(aln),
                identity = identity,
                ref_window = win$ref_window,
                consensus_window = win$query_window,
                cut_site_in_amplicon = cut_in_amp,
                support_reads = consensus_df$support_reads[i],
                fraction = consensus_df$fraction[i],
                stringsAsFactors = FALSE
            )
        }
    }
    out <- do.call(rbind, rows)
    write.csv(out, file.path(out_dir, "consensus_alignments.csv"), row.names = FALSE)
    invisible(out)
}

#' Classify refless genotype string into editViewer evaluation levels
#' @keywords internal
.refless_genotype_to_eval <- function(genotype){
    if(is.na(genotype) || genotype == "ref"){
        return("ref")
    }
    nums <- as.numeric(unlist(regmatches(genotype, gregexpr("[0-9]+", genotype))))
    inframe <- length(nums) > 0 && all(nums %% 3 == 0)
    if(grepl("^sub", genotype)){
        return("alt")
    }
    if(inframe){
        return("alt_inframe_homo")
    }
    "alt"
}

#' Build editViewer-compatible summary table from refless edit-calling results
#' @keywords internal
.build_refless_editcall_summary <- function(per_sample){
    rows <- list()
    for(sid in names(per_sample)){
        edits <- per_sample[[sid]]$edits
        if(is.null(edits) || nrow(edits) == 0){
            next
        }
        primary <- per_sample[[sid]]$selection$primary
        for(i in seq_len(nrow(edits))){
            is_primary <- FALSE
            if(!is.null(primary) && nrow(primary)){
                is_primary <- any(
                    primary$target_gene == edits$target_gene[i] &
                        primary$primary_cluster == edits$cluster_id[i]
                )
            }
            rows[[length(rows) + 1]] <- data.frame(
                sample_id = sid,
                index_pair_id = sid,
                cluster_id = edits$cluster_id[i],
                target_gene = edits$target_gene[i],
                genotype = edits$genotype[i],
                edit_pattern = edits$edit_pattern[i],
                vs_intact = edits$vs_intact[i],
                cluster_fraction = edits$cluster_fraction[i],
                cluster_support = edits$cluster_support[i],
                is_primary = is_primary,
                data_type = "genotype",
                stringsAsFactors = FALSE
            )
        }
    }
    if(length(rows) == 0){
        return(data.frame())
    }
    do.call(rbind, rows)
}

#' Classify edit pattern from consensus vs intact window
#' @keywords internal
.classify_edit_window <- function(consensus_window, ref_window, intact_ref){
    cw <- gsub("-", "", consensus_window)
    rw <- gsub("-", "", ref_window)
    if(cw == rw || cw == intact_ref){
        return(list(pattern = "intact", genotype = "ref", edit_signal = 0))
    }
    n_del <- length(unlist(gregexpr("-", consensus_window)))
    n_ins <- nchar(cw) - nchar(rw)
    if(n_del > 0 && n_ins > 0){
        genotype <- paste0("indel", n_ins, "-", n_del)
        pattern <- paste0("indel at gRNA site (", n_ins, " ins, ", n_del, " del)")
        signal <- n_del + abs(n_ins)
    } else if(n_del > 0){
        genotype <- paste0("del", n_del)
        pattern <- paste0("-", n_del, "bp deletion at gRNA site")
        signal <- n_del
    } else if(n_ins > 0){
        genotype <- paste0("ins", n_ins)
        pattern <- paste0("+", n_ins, "bp insertion at gRNA site")
        signal <- n_ins
    } else {
        genotype <- "sub"
        pattern <- "substitution at gRNA site"
        signal <- 1
    }
    list(pattern = pattern, genotype = genotype, edit_signal = signal)
}

#' Call edits from aligned consensus sequences
#'
#' @param align_df Output from alignConsensusToRef.
#' @param ref_targets Output from prepRefTargets.
#' @param confidence_df Optional confidence metrics from evalConfidence.
#' @param out_dir Output directory.
#' @return Data frame of edit calls per sample/cluster/target.
#' @export
callEditsFromConsensus <- function(align_df,
                                   ref_targets,
                                   confidence_df = NULL,
                                   out_dir){
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    intact <- ref_targets$intact_seq

    rows <- list()
    for(i in seq_len(nrow(align_df))){
        target <- align_df$target_gene[i]
        ref_seq <- as.character(intact[[target]])
        cls <- .classify_edit_window(
            align_df$consensus_window[i],
            align_df$ref_window[i],
            ref_seq
        )
        conf <- NA_real_
        if(!is.null(confidence_df)){
            m <- confidence_df$cluster_id == align_df$cluster_id[i]
            if(any(m)){
                conf <- confidence_df$mean_identity[m][1]
            }
        }
        rows[[i]] <- data.frame(
            sample_id = align_df$sample_id[i],
            cluster_id = align_df$cluster_id[i],
            target_gene = target,
            edit_pattern = cls$pattern,
            genotype = cls$genotype,
            vs_intact = if(cls$genotype == "ref") "wild-type" else "edited",
            alignment_identity = align_df$identity[i],
            cluster_support = align_df$support_reads[i],
            cluster_fraction = align_df$fraction[i],
            edit_signal_strength = cls$edit_signal,
            consensus_confidence = conf,
            stringsAsFactors = FALSE
        )
    }
    out <- do.call(rbind, rows)
    write.csv(out, file.path(out_dir, "editcall_per_cluster.csv"), row.names = FALSE)
    invisible(out)
}

#' Select representative edit cluster per sample
#'
#' @param edit_df Output from callEditsFromConsensus for one sample.
#' @param prefer_edited If TRUE, boost clusters with edit signal.
#' @param weights Named numeric weights for scoring components.
#' @return List with primary cluster selection and full scoring table.
#' @export
selectEditCluster <- function(edit_df,
                              prefer_edited = TRUE,
                              weights = c(
                                  log_fraction = 1.0,
                                  confidence = 0.5,
                                  edit_signal = 0.8,
                                  alignment = 0.3,
                                  small_no_edit_penalty = 0.5
                              )){
    if(nrow(edit_df) == 0){
        return(list(primary = NULL, scores = edit_df, summary = NULL))
    }

    score_one <- function(df){
        w <- weights
        df$score <- w["log_fraction"] * log(pmax(df$cluster_fraction, 1e-6))
        df$score <- df$score + w["confidence"] * ifelse(is.na(df$consensus_confidence), 0, df$consensus_confidence)
        df$score <- df$score + w["edit_signal"] * df$edit_signal_strength * if(prefer_edited) 1 else 0
        df$score <- df$score + w["alignment"] * df$alignment_identity
        penalty <- (df$cluster_fraction < 0.05) & (df$edit_signal_strength == 0)
        df$score <- df$score - w["small_no_edit_penalty"] * penalty
        df[order(-df$score), ]
    }

    pieces <- split(edit_df, edit_df$target_gene)
    summaries <- list()
  scores_all <- list()
    for(tg in names(pieces)){
        scored <- score_one(pieces[[tg]])
        primary <- scored[1, ]
        summaries[[tg]] <- data.frame(
            target_gene = tg,
            primary_cluster = primary$cluster_id,
            primary_score = primary$score,
            primary_edit_pattern = primary$edit_pattern,
            primary_vs_intact = primary$vs_intact,
            confidence = if(primary$alignment_identity >= 0.95) "high" else "moderate",
            stringsAsFactors = FALSE
        )
        scores_all[[tg]] <- scored
    }

    list(
        primary = do.call(rbind, summaries),
        scores = do.call(rbind, scores_all),
        summary = summaries
    )
}

#' Run consensus-based edit calling for pipeline results
#'
#' @param pipeline_results Output from miaoPipeline.
#' @param genome_fn Reference genome FASTA.
#' @param primer_list Primer CSV.
#' @param pam_list gRNA/PAM CSV.
#' @param out_dir Output root for edit calling.
#' @param sif_path Optional Apptainer SIF.
#' @param cut_site_margin bp margin around cut site.
#' @param prefer_edited Prefer edited clusters in selectEditCluster.
#' @return List with ref_targets, per-sample edit calls, and primary selections.
#' @export
doEditCalling <- function(pipeline_results,
                          genome_fn,
                          primer_list,
                          pam_list,
                          out_dir,
                          sif_path = NULL,
                          cut_site_margin = 3,
                          prefer_edited = TRUE){
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    ref_dir <- file.path(out_dir, "ref_targets")
    ref_targets <- prepRefTargets(
        genome_fn = genome_fn,
        primer_list = primer_list,
        pam_list = pam_list,
        out_dir = ref_dir,
        sif_path = sif_path
    )

    per_sample <- list()
    primary_all <- list()

    for(sid in names(pipeline_results$samples)){
        sample_res <- pipeline_results$samples[[sid]]
        sample_out <- file.path(out_dir, "samples", sid)
        align_dir <- file.path(sample_out, "align")
        edit_dir <- file.path(sample_out, "editcall")

        align_df <- alignConsensusToRef(
            consensus_df = sample_res$consensus,
            ref_targets = ref_targets,
            out_dir = align_dir,
            cut_site_margin = cut_site_margin,
            sample_id = sid
        )
        edit_df <- callEditsFromConsensus(
            align_df = align_df,
            ref_targets = ref_targets,
            confidence_df = sample_res$confidence,
            out_dir = edit_dir
        )
        sel <- selectEditCluster(edit_df, prefer_edited = prefer_edited)
        if(!is.null(sel$primary)){
            write.csv(sel$primary, file.path(sample_out, "primary_cluster.csv"), row.names = FALSE)
            write.csv(sel$scores, file.path(sample_out, "cluster_scores.csv"), row.names = FALSE)
            sel$primary$sample_id <- sid
            primary_all[[sid]] <- sel$primary
        }
        per_sample[[sid]] <- list(
            align = align_df,
            edits = edit_df,
            selection = sel
        )
    }

    result <- list(
        ref_targets = ref_targets,
        per_sample = per_sample,
        primary = if(length(primary_all)) do.call(rbind, primary_all) else NULL
    )
    saveRDS(result, file.path(out_dir, "doEditCalling_results.rds"))

    summary_df <- .build_refless_editcall_summary(per_sample)
    if(nrow(summary_df) > 0){
        editcall_dir <- file.path(out_dir, "editcall")
        dir.create(editcall_dir, showWarnings = FALSE, recursive = TRUE)
        write.csv(summary_df, file.path(editcall_dir, "editcall_summary.csv"), row.names = FALSE)
    }

    invisible(result)
}
