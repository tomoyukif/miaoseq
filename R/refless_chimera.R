################################################################################
# Chimera detection (vsearch uchime_denovo)
################################################################################

#' Parse derep UC file: map each read to its derep representative
#' @keywords internal
.parse_derep_uc <- function(uc_fn){
    lines <- readLines(uc_fn)
    mapping <- list()
    for(ln in lines){
        x <- strsplit(ln, "\t", fixed = TRUE)[[1]]
        if(length(x) < 10){
            next
        }
        typ <- x[1]
        if(typ == "S"){
            mapping[[x[9]]] <- x[9]
        } else if(typ == "H"){
            mapping[[x[9]]] <- x[10]
        }
    }
    if(length(mapping) == 0){
        return(data.frame(
            read_id = character(),
            rep_id = character(),
            stringsAsFactors = FALSE
        ))
    }
    data.frame(
        read_id = names(mapping),
        rep_id = unlist(mapping, use.names = FALSE),
        stringsAsFactors = FALSE
    )
}

#' Run vsearch derep + uchime_denovo
#' @keywords internal
.run_uchime_denovo <- function(fasta_fn, out_dir, abskew = 2.0, sif_path = NULL){
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    local_fa <- file.path(out_dir, "reads.fa")
    if(normalizePath(dirname(fasta_fn)) != normalizePath(out_dir)){
        file.copy(fasta_fn, local_fa, overwrite = TRUE)
    } else {
        local_fa <- fasta_fn
    }
    derep_fa <- "derep.fa"
    derep_uc <- "derep.uc"
    chim_fa <- "chimeras.fa"
    good_fa <- "nonchimeras.fa"

    derep_args <- paste(
        "--derep_fulllength reads.fa",
        "--output", derep_fa,
        "--sizeout",
        "--uc", derep_uc,
        "--minuniquesize 1"
    )
    status <- .refless_run_tool("vsearch", derep_args, work_dir = out_dir, sif_path = sif_path, prefer_host = FALSE)
    if(!is.null(status) && status != 0){
        stop("vsearch --derep_fulllength failed")
    }

    uchime_args <- paste(
        "--uchime_denovo", derep_fa,
        "--chimeras", chim_fa,
        "--nonchimeras", good_fa,
        "--abskew", abskew
    )
    status <- .refless_run_tool("vsearch", uchime_args, work_dir = out_dir, sif_path = sif_path, prefer_host = FALSE)
    if(!is.null(status) && status != 0){
        stop("vsearch --uchime_denovo failed")
    }

    list(
        derep_uc = file.path(out_dir, derep_uc),
        chimeras = file.path(out_dir, chim_fa)
    )
}

#' Detect PCR chimeras among reads using vsearch uchime_denovo
#'
#' Dereplicates reads, runs de novo chimera detection, and maps chimera labels
#' back to individual read IDs.
#'
#' @param reads_fasta Sample FASTA/FASTQ.
#' @param out_dir Output directory.
#' @param abskew Abundance skew parameter for uchime_denovo (default 2.0).
#' @param sif_path Optional Apptainer image path.
#' @param exclude_from_clusters If TRUE, return only non-chimera read IDs.
#' @return Data frame with read_id, rep_id, is_chimera.
#' @export
detectChimeras <- function(reads_fasta,
                           out_dir,
                           abskew = 2.0,
                           sif_path = NULL,
                           exclude_from_clusters = FALSE){
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    fa_in <- reads_fasta
    if(grepl("\\.(fq|fastq)(\\.gz)?$", reads_fasta, ignore.case = TRUE)){
        fa_in <- file.path(out_dir, "reads.fa")
        seqs <- .refless_read_seqs(reads_fasta)
        Biostrings::writeXStringSet(seqs, fa_in, width = 80)
    }

    uchime <- .run_uchime_denovo(fa_in, out_dir, abskew = abskew, sif_path = sif_path)
    map_df <- .parse_derep_uc(uchime$derep_uc)

    chimera_ids <- character()
    if(file.exists(uchime$chimeras)){
        chim_lines <- readLines(uchime$chimeras, warn = FALSE)
        if(length(chim_lines) > 0 && any(nzchar(chim_lines))){
            chim_seqs <- Biostrings::readDNAStringSet(uchime$chimeras)
            chimera_ids <- names(chim_seqs)
            chimera_ids <- sub(";size=.*", "", chimera_ids)
        }
    }

    map_df$is_chimera <- map_df$rep_id %in% chimera_ids
    write.csv(map_df, file.path(out_dir, "chimera_assignments.csv"), row.names = FALSE)

    summary <- data.frame(
        total_reads = nrow(map_df),
        chimera_reads = sum(map_df$is_chimera),
        chimera_fraction = if(nrow(map_df)) sum(map_df$is_chimera) / nrow(map_df) else 0,
        stringsAsFactors = FALSE
    )
    write.csv(summary, file.path(out_dir, "chimera_summary.csv"), row.names = FALSE)

    if(exclude_from_clusters){
        map_df <- map_df[!map_df$is_chimera, , drop = FALSE]
    }
    invisible(map_df)
}
