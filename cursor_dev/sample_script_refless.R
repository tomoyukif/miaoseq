#!/usr/bin/env Rscript
# Sample script for refless core pipeline (demux -> cluster -> consensus)
#
# Prerequisites:
#   1. Build Apptainer image (optional, for vsearch/umap when host tools missing):
#        bash cursor_dev/apptainer/build.sh
#   2. export MIAOSEQ_REFLESS_SIF=cursor_dev/apptainer/images/miaoseq-refless.sif
#   3. Optional host overrides:
#        export MIAOSEQ_MMSEQS_PATH=/path/to/mmseqs
#        export MIAOSEQ_MINIMAP2_PATH=/path/to/minimap2

library(miaoseq)

fastq_fn <- Sys.getenv("REFLESS_FASTQ", "path/to/multiplexed.fastq.gz")
out_dir <- Sys.getenv("REFLESS_OUT", "refless_output")
index_list <- system.file("extdata", "index_list.csv", package = "miaoseq")
primer_list <- system.file("extdata", "amplicon_primers.csv", package = "miaoseq")
sif_default <- "cursor_dev/apptainer/images/miaoseq-refless.sif"
sif_path <- Sys.getenv("MIAOSEQ_REFLESS_SIF", sif_default)

results <- miaoPipeline(
    fastq_fn = fastq_fn,
    index_list = index_list,
    primer_list = primer_list,
    out_dir = out_dir,
    cluster_method = "mmseqs2",
    sif_path = if(file.exists(sif_path)) sif_path else NULL,
    detect_chimeras = FALSE,
    write_report = TRUE,
    demultiplex_params = list(
        end_window = 200,
        max_reads = NULL
    ),
    cluster_params = list(
        refine = TRUE,
        mmseqs_min_seq_id = 0.85
    ),
    resume = FALSE
)

print(table(results$demultiplex$class))
cat("Samples processed:", length(results$samples), "\n")
if(!is.null(results$report)){
    cat("Quality report:", results$report$html_report, "\n")
}
