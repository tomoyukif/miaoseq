################################################################################
# Read clustering for refless pipeline
################################################################################

#' Parse MMseqs2 easy-linclust cluster TSV
#' @keywords internal
.parse_mmseqs_cluster_tsv <- function(tsv_fn){
    if(!file.exists(tsv_fn)){
        stop("MMseqs2 cluster output not found: ", tsv_fn)
    }
    x <- read.table(tsv_fn, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    names(x) <- c("cluster_rep", "read_id")
    x$cluster_rep <- sub(";.*", "", x$cluster_rep)
    x
}

#' Parse vsearch UC output into read->cluster mapping
#' @keywords internal
.parse_vsearch_uc <- function(uc_fn, read_ids){
    if(!file.exists(uc_fn)){
        stop("vsearch UC output not found: ", uc_fn)
    }
    lines <- readLines(uc_fn)
    mapping <- character(length(read_ids))
    names(mapping) <- read_ids
    for(ln in lines){
        x <- strsplit(ln, "\t", fixed = TRUE)[[1]]
        if(length(x) < 10){
            next
        }
        typ <- x[1]
        if(typ == "S"){
            mapping[x[9]] <- x[9]
        } else if(typ == "H"){
            mapping[x[9]] <- x[10]
        }
    }
    data.frame(
        read_id = names(mapping),
        cluster_id = unname(mapping),
        stringsAsFactors = FALSE
    )
}

#' Run MMseqs2 easy-linclust via Apptainer or host
#' @keywords internal
.run_mmseqs_cluster <- function(fasta_fn,
                                 out_dir,
                                 min_seq_id = 0.85,
                                 c_cov = 0.8,
                                 sif_path = NULL){
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    local_fa <- file.path(out_dir, "reads.fa")
    if(normalizePath(dirname(fasta_fn)) != normalizePath(out_dir)){
        file.copy(fasta_fn, local_fa, overwrite = TRUE)
    } else {
        local_fa <- fasta_fn
    }
    tmp <- "mmseqs_tmp"
    prefix <- "clusterResult"
    args <- paste(
        "easy-linclust reads.fa", prefix, tmp,
        "--min-seq-id", min_seq_id,
        "-c", c_cov,
        "--cov-mode", 0,
        "--alignment-mode", 3
    )
    status <- .refless_run_tool("mmseqs", args, work_dir = out_dir, sif_path = sif_path)
    if(!is.null(status) && status != 0){
        stop("mmseqs easy-linclust failed")
    }
    file.path(out_dir, paste0(prefix, "_cluster.tsv"))
}

#' Run vsearch cluster_fast (MeShClust stand-in when meshclust binary unavailable)
#' @keywords internal
.run_vsearch_cluster <- function(fasta_fn,
                                 out_dir,
                                 id = 0.90,
                                 sif_path = NULL){
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    local_fa <- file.path(out_dir, "reads.fa")
    if(normalizePath(dirname(fasta_fn)) != normalizePath(out_dir)){
        file.copy(fasta_fn, local_fa, overwrite = TRUE)
    } else {
        local_fa <- fasta_fn
    }
    clstr <- file.path(out_dir, "vsearch_clusters.uc")
    centroids <- file.path(out_dir, "vsearch_centroids.fa")
    args <- paste(
        "--cluster_fast reads.fa",
        "--id", id,
        "--centroids vsearch_centroids.fa",
        "--uc vsearch_clusters.uc",
        "--clusterout_sort"
    )
    status <- .refless_run_tool("vsearch", args, work_dir = out_dir, sif_path = sif_path)
    if(!is.null(status) && status != 0){
        stop("vsearch --cluster_fast failed")
    }
    list(uc = clstr, centroids = centroids)
}

#' Run UMAP pre-cluster assignment script on host or in container
#' @keywords internal
.run_umap_prefilter <- function(fasta_fn, out_dir, min_samples = 50, sif_path = NULL){
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    out_tsv <- file.path(out_dir, "umap_groups.tsv")
    status <- .refless_run_umap_script(
        fasta_fn = fasta_fn,
        out_tsv = out_tsv,
        min_samples = min_samples,
        sif_path = sif_path
    )
    if(!is.null(status) && status != 0){
        stop("UMAP pre-cluster script failed")
    }
    read.table(out_tsv, sep = "\t", header = FALSE,
               col.names = c("read_id", "umap_group"), stringsAsFactors = FALSE)
}

#' Cluster reads within a sample
#'
#' Supports three methods: meshclust (vsearch cluster_fast), mmseqs2, umap_meshclust.
#'
#' @param reads_fasta Path to sample FASTA/FASTQ.
#' @param out_dir Output directory for clustering results.
#' @param cluster_method One of meshclust, mmseqs2, umap_meshclust.
#' @param refine If TRUE, run minimap2-based cluster refinement (placeholder).
#' @param meshclust_id Identity threshold for meshclust/vsearch method.
#' @param mmseqs_min_seq_id MMseqs2 --min-seq-id.
#' @param umap_min_samples UMAP/HDBSCAN min cluster size.
#' @param sif_path Path to Apptainer SIF image.
#' @param detect_chimeras If TRUE, run \code{detectChimeras()} and exclude chimeras from clustering.
#' @param chimera_abskew Abundance skew for uchime_denovo.
#' @return Data frame with read_id and cluster_id.
#' @export
doCluster <- function(reads_fasta,
                      out_dir,
                      cluster_method = c("meshclust", "mmseqs2", "umap_meshclust"),
                      refine = TRUE,
                      meshclust_id = 0.90,
                      mmseqs_min_seq_id = 0.85,
                      umap_min_samples = 50,
                      sif_path = NULL,
                      detect_chimeras = FALSE,
                      chimera_abskew = 2.0){
    cluster_method <- match.arg(cluster_method)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    chimera_df <- NULL
    if(detect_chimeras){
        chim_dir <- file.path(out_dir, "chimera")
        chimera_df <- detectChimeras(
            reads_fasta = reads_fasta,
            out_dir = chim_dir,
            abskew = chimera_abskew,
            sif_path = sif_path
        )
    }

    fa_in <- reads_fasta
    if(grepl("\\.(fq|fastq)(\\.gz)?$", reads_fasta, ignore.case = TRUE)){
        fa_in <- file.path(out_dir, "reads.fa")
        seqs <- .refless_read_seqs(reads_fasta)
        Biostrings::writeXStringSet(seqs, fa_in, width = 80)
    }

    if(cluster_method == "mmseqs2"){
        tsv <- .run_mmseqs_cluster(
            fa_in, out_dir,
            min_seq_id = mmseqs_min_seq_id,
            sif_path = sif_path
        )
        cl <- .parse_mmseqs_cluster_tsv(tsv)
        cl$cluster_id <- cl$cluster_rep
        out <- cl[, c("read_id", "cluster_id")]

    } else if(cluster_method == "meshclust"){
        vs <- .run_vsearch_cluster(fa_in, out_dir, id = meshclust_id, sif_path = sif_path)
        seqs <- Biostrings::readDNAStringSet(fa_in)
        read_ids <- names(seqs)
        if(is.null(read_ids)){
            read_ids <- paste0("read_", seq_along(seqs))
        }
        out <- .parse_vsearch_uc(vs$uc, read_ids)

    } else if(cluster_method == "umap_meshclust"){
        groups <- .run_umap_prefilter(fa_in, out_dir, min_samples = umap_min_samples, sif_path = sif_path)
        seqs <- Biostrings::readDNAStringSet(fa_in)
        read_ids <- names(seqs)
        if(is.null(read_ids)){
            read_ids <- paste0("read_", seq_along(seqs))
        }
        pieces <- split(read_ids, groups$umap_group[match(read_ids, groups$read_id)])
        out_list <- list()
        for(g in names(pieces)){
            g_dir <- file.path(out_dir, paste0("umap_group_", g))
            g_fa <- file.path(g_dir, "reads.fa")
            dir.create(g_dir, recursive = TRUE, showWarnings = FALSE)
            Biostrings::writeXStringSet(seqs[pieces[[g]]], g_fa)
            tsv <- .run_mmseqs_cluster(g_fa, g_dir, min_seq_id = mmseqs_min_seq_id, sif_path = sif_path)
            cl <- .parse_mmseqs_cluster_tsv(tsv)
            cl$cluster_id <- paste0("ug", g, "_", cl$cluster_rep)
            out_list[[g]] <- cl[, c("read_id", "cluster_id")]
        }
        out <- do.call(rbind, out_list)
    }

    if(refine){
        out <- .refine_clusters_minimap(out, fa_in, out_dir, sif_path = sif_path)
    }

    if(detect_chimeras && !is.null(chimera_df)){
        out <- out[!out$read_id %in% chimera_df$read_id[chimera_df$is_chimera], , drop = FALSE]
    }

    write.csv(out, file.path(out_dir, "cluster_assignments.csv"), row.names = FALSE)
    invisible(out)
}

#' Refine clusters by re-grouping reads with high minimap2 identity to cluster medoid
#' @keywords internal
.refine_clusters_minimap <- function(cluster_df, fasta_fn, out_dir, sif_path = NULL){
    seqs <- Biostrings::readDNAStringSet(fasta_fn)
    clusters <- split(cluster_df$read_id, cluster_df$cluster_id)
    refined <- list()
    for(cid in names(clusters)){
        ids <- clusters[[cid]]
        if(length(ids) < 2){
            refined[[cid]] <- data.frame(read_id = ids, cluster_id = cid, stringsAsFactors = FALSE)
            next
        }
        sub <- seqs[ids]
        lens <- Biostrings::width(sub)
        medoid <- names(sub)[which.max(lens)]
        medoid_fa <- file.path(out_dir, paste0("medoid_", cid, ".fa"))
        Biostrings::writeXStringSet(sub[medoid], medoid_fa)
        query_fa <- file.path(out_dir, paste0("cluster_", cid, ".fa"))
        Biostrings::writeXStringSet(sub, query_fa)
        paf <- file.path(out_dir, paste0("cluster_", cid, ".paf"))
        cmd <- sprintf(
            "minimap2 -x map-ont %s %s > %s",
            shQuote(basename(query_fa)),
            shQuote(basename(medoid_fa)),
            shQuote(basename(paf))
        )
        status <- .refless_run_shell(cmd, work_dir = out_dir, sif_path = sif_path)
        if(is.null(status) || status != 0 || !file.exists(paf) || !nzchar(readLines(paf, n = 1, warn = FALSE))){
            refined[[cid]] <- data.frame(read_id = ids, cluster_id = cid, stringsAsFactors = FALSE)
            next
        }
        paf_lines <- readLines(paf)
        map_id <- vapply(strsplit(paf_lines, "\t"), function(x) x[1], character(1))
        map_match <- vapply(strsplit(paf_lines, "\t"), function(x) as.integer(x[10]), integer(1))
        keep <- map_id[map_match >= 80]
        if(length(keep) == 0){
            keep <- ids
        }
        refined[[cid]] <- data.frame(read_id = keep, cluster_id = cid, stringsAsFactors = FALSE)
    }
    do.call(rbind, refined)
}
