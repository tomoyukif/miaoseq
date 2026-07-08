################################################################################
# Step 1: Basecalling
################################################################################

#' Create minimap2 index for Dorado reference-guided basecalling
#'
#' @param genome_fn Path to the reference genome FASTA file.
#' @param mmi_fn Output path for the minimap2 index (.mmi).
#' @param minimap2_path Path to the minimap2 executable.
#' @export
makeMMI <- function(genome_fn, mmi_fn, minimap2_path){
    system2(command = minimap2_path,
            args = paste("-x map-ont -d", mmi_fn, genome_fn))
}

#' Basecalling using Dorado
#'
#' This function performs basecalling on raw sequencing data using the Dorado basecaller.
#' It processes the raw data in a specified input directory, applies quality filtering,
#' and size selection, and outputs the basecalled reads in FASTA format.
#' @param in_dir Path to the input directory containing raw sequencing data in pod5 format.
#' @param basecall_dir Path to the output directory where basecalling results will be saved.
#' @param size_sel A numeric vector of length 2 specifying the minimum and maximum read lengths (in bp) to retain after size selection.
#' @param dorado_path Path to the Dorado executable for basecalling.
#' @param samtools_path Path to the Samtools executable for BAM file processing.
#' @param mmi_fn Optional path to a minimap2 index for reference-guided basecalling.
#' @param bed_fn Optional BED file for reference-guided basecalling (requires `mmi_fn`).
#' @param n_core Number of CPU cores to use for parallel processing.
#' @return A character vector containing paths to the size-selected FASTA files of basecalled reads.
#' @import Biostrings
#' @importFrom BiocGenerics width
#' @importFrom Biostrings readDNAStringSet writeXStringSet DNAStringSet
#' @importFrom utils write.table
#' @importFrom stats setNames
#' @importFrom methods as
#' @export
#'
doBasecall <- function(in_dir,
                       basecall_dir,
                       size_sel = c(0, Inf),
                       dorado_path,
                       samtools_path,
                       mmi_fn = NULL,
                       bed_fn = NULL,
                       n_core){
    bam_fn <- file.path(basecall_dir, "basecall.bam")
    args <- paste("duplex sup", in_dir,
                  "--threads", n_core,
                  "--min-qscore 10")
    if(!is.null(mmi_fn)){
        args <- paste(args, "--reference", mmi_fn)
        if(!is.null(bed_fn)){
            args <- paste(args, "--bed-file", bed_fn)
        }
    }

    system2(command = dorado_path,
            args = paste(args, ">", bam_fn))

    summary_out <- file.path(basecall_dir, "basecalls_summary.tsv")
    system2(command = dorado_path,
            args = paste("summary", bam_fn,
                         ">", summary_out
            ))

    filt_bam <- file.path(basecall_dir, "basecall_filt.bam")
    system2(command = samtools_path,
            args = paste("view -h", bam_fn,
                         "|",
                         "awk '!/\tdx:i:-1\b/' |",
                         samtools_path,
                         " view -bS - >",
                         filt_bam
            ))

    fq_fn <- file.path(basecall_dir, "basecall_filt.fq")
    system2(command = samtools_path,
            args = paste("fastq", filt_bam,
                         ">", fq_fn
            ))

    reads <- readDNAStringSet(fq_fn, format = "fastq")
    reads <- reads[width(reads) >= size_sel[1] & width(reads) <= size_sel[2]]
    start <- seq(1, length(reads), by = 1e4)
    end <- c(tail(start, -1) - 1, length(reads))
    n_digit <- nchar(length(start))
    fq_fn <- NULL
    for(i in seq_along(start)){
        i_fq_fn <- file.path(basecall_dir,
                             paste0("basecall_filt_sizeselected_reads_",
                                    sprintf(paste0("%0", n_digit, "d"), i),
                                    ".fa"))
        writeXStringSet(reads[start[i]:end[i]], i_fq_fn)
        fq_fn <- c(fq_fn, i_fq_fn)
    }
    write("", file = file.path(basecall_dir, "basecall.finish"))
    return(fq_fn)
}
