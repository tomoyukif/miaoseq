################################################################################
# miaoPipeline: core refless orchestration
################################################################################

#' Run the refless core pipeline (demux -> cluster -> consensus -> confidence)
#'
#' @param fastq_fn Path to multiplexed FASTQ.
#' @param index_list Index CSV (5 columns).
#' @param primer_list Primer CSV (_F / _R).
#' @param out_dir Output root directory.
#' @param cluster_method Clustering method: meshclust, mmseqs2, or umap_meshclust.
#' @param demultiplex_params Named list passed to doDemultiplex2.
#' @param cluster_params Named list passed to doCluster.
#' @param sif_path Path to miaoseq-refless Apptainer image.
#' @param samples Optional character vector of sample IDs to process (default all assigned).
#' @param resume Skip steps when outputs exist.
#' @param detect_chimeras If TRUE, run chimera detection per sample (passed to doCluster).
#' @param write_report If TRUE, generate quality report after processing.
#' @param report_params Named list passed to writeReflessReport.
#' @return List with demultiplex, per-sample cluster/consensus/confidence results.
#' @export
miaoPipeline <- function(fastq_fn,
                         index_list,
                         primer_list,
                         out_dir,
                         cluster_method = "mmseqs2",
                         demultiplex_params = list(),
                         cluster_params = list(),
                         sif_path = NULL,
                         samples = NULL,
                         resume = FALSE,
                         detect_chimeras = FALSE,
                         write_report = TRUE,
                         report_params = list()){
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    demult_dir <- file.path(out_dir, "demultiplex")
    demult_fn <- file.path(demult_dir, "demultiplex_assignments.csv")

    if(resume && file.exists(demult_fn)){
        demult_out <- read.csv(demult_fn, stringsAsFactors = FALSE)
    } else {
        demult_args <- c(
            list(
                fastq_fn = fastq_fn,
                index_list = index_list,
                primer_list = primer_list,
                out_dir = demult_dir,
                write_per_sample = TRUE
            ),
            demultiplex_params
        )
        demult_out <- do.call(doDemultiplex2, demult_args)
    }

    sample_dir <- file.path(demult_dir, "per_sample")
    if(is.null(samples)){
        samples <- unique(demult_out$sample_id[demult_out$class %in% c("High", "Low")])
        samples <- samples[!is.na(samples)]
    }

    results <- list(demultiplex = demult_out, samples = list())

    for(sid in samples){
        sample_fq <- file.path(sample_dir, paste0(sid, ".fastq"))
        if(!file.exists(sample_fq)){
            next
        }
        sample_out <- file.path(out_dir, "samples", sid)
        cluster_dir <- file.path(sample_out, "cluster")
        consensus_dir <- file.path(sample_out, "consensus")
        confidence_dir <- file.path(sample_out, "confidence")

        cluster_fn <- file.path(cluster_dir, "cluster_assignments.csv")
        if(resume && file.exists(cluster_fn)){
            cluster_df <- read.csv(cluster_fn, stringsAsFactors = FALSE)
        } else {
            cluster_args <- c(
                list(
                    reads_fasta = sample_fq,
                    out_dir = cluster_dir,
                    cluster_method = cluster_method,
                    sif_path = sif_path,
                    detect_chimeras = detect_chimeras
                ),
                cluster_params
            )
            cluster_df <- do.call(doCluster, cluster_args)
        }

        chimera_fn <- file.path(cluster_dir, "chimera", "chimera_assignments.csv")
        chimera_df <- if(file.exists(chimera_fn)){
            read.csv(chimera_fn, stringsAsFactors = FALSE)
        } else {
            NULL
        }

        consensus_fn <- file.path(consensus_dir, "consensus_summary.csv")
        if(resume && file.exists(consensus_fn)){
            consensus_df <- read.csv(consensus_fn, stringsAsFactors = FALSE)
        } else {
            consensus_df <- doConsensus(sample_fq, cluster_df, consensus_dir)
        }

        confidence_fn <- file.path(confidence_dir, "confidence_metrics.csv")
        if(resume && file.exists(confidence_fn)){
            confidence_df <- read.csv(confidence_fn, stringsAsFactors = FALSE)
        } else {
            confidence_df <- evalConfidence(sample_fq, cluster_df, consensus_df, confidence_dir)
        }

        results$samples[[sid]] <- list(
            cluster = cluster_df,
            consensus = consensus_df,
            confidence = confidence_df,
            chimera = chimera_df
        )
    }

    saveRDS(results, file.path(out_dir, "miaoPipeline_results.rds"))

    if(write_report){
        report_args <- c(
            list(
                pipeline_results = results,
                pipeline_out = out_dir
            ),
            report_params
        )
        results$report <- do.call(writeReflessReport, report_args)
    }

    invisible(results)
}
