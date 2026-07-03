################################################################################
# Consensus generation and confidence evaluation
################################################################################

#' Majority consensus from aligned sequences (simple column-wise vote)
#' @keywords internal
.pileup_consensus <- function(seqs){
    if(length(seqs) == 0){
        return(NULL)
    }
    aln <- pwalign::pairwiseAlignment(
        pattern = seqs,
        subject = seqs[[1]],
        type = "global",
        gapOpening = 1,
        gapExtension = 1
    )
    mat <- do.call(rbind, strsplit(as.character(pwalign::aligned(aln)), ""))
    cons <- apply(mat, 2, function(col){
        tab <- sort(table(col[col != "-"]), decreasing = TRUE)
        names(tab)[1]
    })
    gsub("-", "", paste(cons, collapse = ""))
}

#' Per-position support rates from global alignment to consensus
#' @keywords internal
.pileup_support <- function(seqs, consensus){
    if(length(seqs) < 2){
        return(data.frame(position = integer(), base = character(), support = numeric()))
    }
    aln <- pwalign::pairwiseAlignment(
        pattern = seqs,
        subject = consensus,
        type = "global",
        gapOpening = 1,
        gapExtension = 1
    )
    mat <- do.call(rbind, strsplit(as.character(pwalign::aligned(aln)), ""))
    n <- ncol(mat)
    rows <- list()
    pos <- 0
    for(j in seq_len(n)){
        col <- mat[, j]
        if(any(col == "-")){
            next
        }
        pos <- pos + 1
        tab <- table(col)
        for(b in names(tab)){
            rows[[length(rows) + 1]] <- data.frame(
                position = pos,
                base = b,
                support = 100 * tab[b] / length(seqs),
                stringsAsFactors = FALSE
            )
        }
    }
    if(length(rows) == 0){
        return(data.frame(position = integer(), base = character(), support = numeric()))
    }
    do.call(rbind, rows)
}

#' Generate per-cluster consensus sequences
#'
#' @param reads_fasta FASTA/FASTQ of reads for one sample.
#' @param cluster_df Data frame from doCluster with read_id and cluster_id.
#' @param out_dir Output directory.
#' @param trim_primers Optional length to trim from each end after primer removal.
#' @return Data frame with cluster metadata and consensus paths.
#' @export
doConsensus <- function(reads_fasta, cluster_df, out_dir, trim_primers = NULL){
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    seqs <- .refless_read_seqs(reads_fasta)
    if(is.null(names(seqs))){
        names(seqs) <- paste0("read_", seq_along(seqs))
    }

    clusters <- split(cluster_df$read_id, cluster_df$cluster_id)
    total_reads <- nrow(cluster_df)
    results <- list()

    for(cid in names(clusters)){
        ids <- intersect(clusters[[cid]], names(seqs))
        if(length(ids) == 0){
            next
        }
        sub <- seqs[ids]
        if(!is.null(trim_primers) && trim_primers > 0 && all(Biostrings::width(sub) > 2 * trim_primers)){
            sub <- Biostrings::narrow(sub, start = trim_primers + 1, end = Biostrings::width(sub) - trim_primers)
        }
        cons <- .pileup_consensus(as.list(sub))
        cons_dss <- Biostrings::DNAStringSet(cons)
        names(cons_dss) <- cid
        cons_file <- file.path(out_dir, paste0(cid, "_consensus.fa"))
        Biostrings::writeXStringSet(cons_dss, cons_file)
        results[[cid]] <- data.frame(
            cluster_id = cid,
            support_reads = length(ids),
            fraction = length(ids) / total_reads,
            consensus_length = nchar(cons),
            consensus_file = cons_file,
            stringsAsFactors = FALSE
        )
    }
    out <- do.call(rbind, results)
    write.csv(out, file.path(out_dir, "consensus_summary.csv"), row.names = FALSE)
    invisible(out)
}

#' Evaluate confidence metrics for cluster consensus sequences
#'
#' @param reads_fasta FASTA/FASTQ for the sample.
#' @param cluster_df Cluster assignments from doCluster.
#' @param consensus_df Consensus summary from doConsensus.
#' @param out_dir Output directory.
#' @return Data frame with confidence metrics per cluster.
#' @export
evalConfidence <- function(reads_fasta, cluster_df, consensus_df, out_dir){
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    seqs <- .refless_read_seqs(reads_fasta)
    if(is.null(names(seqs))){
        names(seqs) <- paste0("read_", seq_along(seqs))
    }

    metrics <- list()
    for(i in seq_len(nrow(consensus_df))){
        cid <- consensus_df$cluster_id[i]
        cons <- as.character(Biostrings::readDNAStringSet(consensus_df$consensus_file)[[1]])
        ids <- cluster_df$read_id[cluster_df$cluster_id == cid]
        ids <- intersect(ids, names(seqs))
        sub <- seqs[ids]

        identities <- vapply(as.list(sub), function(s){
            aln <- pwalign::pairwiseAlignment(s, cons, type = "global", gapOpening = 1, gapExtension = 1)
            m <- as.character(pwalign::aligned(aln@pattern))
            sum(m != "-") / nchar(cons)
        }, numeric(1))

        support <- .pileup_support(as.list(sub), cons)
        min_support <- if(nrow(support) > 0) min(support$support[support$position > 0]) else NA_real_

        metrics[[cid]] <- data.frame(
            cluster_id = cid,
            support_reads = length(ids),
            fraction = consensus_df$fraction[i],
            consensus_length = nchar(cons),
            mean_identity = mean(identities, na.rm = TRUE),
            min_base_support = min_support,
            stringsAsFactors = FALSE
        )
        write.csv(support, file.path(out_dir, paste0(cid, "_base_support.csv")), row.names = FALSE)
    }
    out <- do.call(rbind, metrics)
    write.csv(out, file.path(out_dir, "confidence_metrics.csv"), row.names = FALSE)
    invisible(out)
}
