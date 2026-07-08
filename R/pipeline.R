################################################################################
# Pipeline orchestration
################################################################################
#' Main pipeline wrapper for MiaoEditcall
#'
#' This function serves as the main pipeline wrapper for MiaoEditcall,
#' orchestrating the steps of basecalling, demultiplexing, read alignment,
#' and edit-calling. It takes in raw sequencing data and processes it
#' through each step, producing a summary of edit-calling results.
#'
#' @param in_dir Path to the input directory containing raw sequencing data in pod5 format.
#' @param out_dir Path to the output directory where results will be saved.
#' @param dorado_path Path to the dorado executable for basecalling.
#' @param samtools_path Path to the samtools executable for BAM file processing.
#' @param blast_path Path to a directory containing BLAST executables (makeblastdb and blastn).
#' @param primer_list Path to a CSV file containing primer sequences.
#' @param pam_list Path to a CSV file containing PAM site information.
#' @param index_list Path to a CSV file containing index sequences for demultiplexing.
#' @param genome_fn Path to the reference genome sequence in FASTA format.
#' @param amplicon_fn Path to the amplicon database FASTA file.
#' @param size_sel A numeric vector of length 2 specifying the minimum and maximum read lengths (in bp) to retain for edit-calling.
#' @param check_window An integer specifying the window size (in bp) around the expected cut site to search for edits.
#' @param n_core Number of CPU cores to use for parallel processing.
#' @param resume Logical indicating whether to resume from a previous run if output files already exist.
#' @param sample_list Optional path to a CSV file mapping index pairs to sample information.
#'   If provided, adds sample names and total read counts to the summary.
#'
#' @export
#'
#' @import Biostrings
#' @import dplyr
#' @importFrom BiocGenerics width
#' @importFrom Biostrings readDNAStringSet writeXStringSet DNAStringSet
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom stats setNames
#' @importFrom utils read.csv write.csv
#' @importFrom tools file_ext
#' @importFrom parallel mclapply
#' @importFrom methods as
#'
#' @examples
#' # in_dir <- "/path/to/pod5/directory" # MinION outputs raw sequence data files in a pod5 directory
#' # out_dir <- "/path/to/your/working/directory/output_directory_name"
#' # genome_fn <- '/reference/genome/sequence.fa' # Reference genome sequence in fasta format
#' # pam_list <- "inst/extdata/agr8_pam_list.csv"
#' # index_list <- "inst/extdata/index_list.csv"
#' # primer_list <- "inst/extdata/amplicon_primers.csv"
#' # dorado_path <- "/path/to/dorado" # Path to dorado executable
#' # samtools_path <- "/path/to/samtools" # Path to samtools executable
#' # blast_path <- "/path/to/blast/bin" # Path to a directory containing blast executables (makeblastdb and blastn)
#' # n_core <- 30 # Number of CPU cores to use
#' # amplicon_fn <- "/path/to/your/working/directory/output_directory_name/ref/amplicon.fa"
#' # editcall_out <- miaoEditcall(in_dir = in_dir,
#' #                              out_dir = out_dir,
#' #                              dorado_path = dorado_path,
#' #                              samtools_path = samtools_path,
#' #                              blast_path = blast_path,
#' #                              primer_list = primer_list,
#' #                              pam_list = pam_list,
#' #                              index_list = index_list,
#' #                              genome_fn = genome_fn,
#' #                              amplicon_fn = amplicon_fn,
#' #                              size_sel = c(300, 450), # Set the valid range of read length for edit-calling in bp. Adjust based on your amplicon size.
#' #                              check_window = 10, # Set the window size (in bp) around the expected cut site to search for edits. Adjust based on your experimental design.
#' #                              n_core = n_core,
#' #                              resume = FALSE) # Set resume = TRUE to resume from previous run
#' @return A data frame summarizing the edit-calling results.
#'
miaoEditcall <- function(in_dir,
                         out_dir,
                         dorado_path = "dorado",
                         samtools_path = "samtools",
                         blast_path = "blastn",
                         primer_list,
                         pam_list,
                         index_list,
                         genome_fn,
                         mmi_fn = NULL,
                         bed_fn = NULL,
                         amplicon_fn,
                         size_sel,
                         check_window = 10,
                         strict = FALSE,
                         n_core = 1,
                         resume = FALSE,
                         sample_list = NULL){
    basecall_dir <- file.path(out_dir, "basecall")
    demult_dir <- file.path(out_dir, "demultiplex")
    align_dir <- file.path(out_dir, "align")
    editcall_dir <- file.path(out_dir, "editcall")

    if(!resume){
        if(dir.exists(basecall_dir)){
            unlink(basecall_dir, recursive = TRUE, force = TRUE)
        }
        if(dir.exists(demult_dir)){
            unlink(demult_dir, recursive = TRUE, force = TRUE)
        }
        if(dir.exists(align_dir)){
            unlink(align_dir, recursive = TRUE, force = TRUE)
        }
        if(dir.exists(editcall_dir)){
            unlink(editcall_dir, recursive = TRUE, force = TRUE)
        }
    }

    dir.create(basecall_dir, recursive = TRUE, showWarnings = FALSE)
    basecall_fn <- file.path(basecall_dir, "basecall.finish")
    if(resume && file.exists(basecall_fn)){
        message("Basecalling has already completed")
        basecall_fn <- list.files(path = basecall_dir,
                                  pattern = "basecall_filt_sizeselected_reads_.+.fa$",
                                  full.names = TRUE)
    } else {
        basecall_fn <- doBasecall(in_dir = in_dir,
                                  basecall_dir = basecall_dir,
                                  size_sel = size_sel,
                                  dorado_path = dorado_path,
                                  samtools_path = samtools_path,
                                  mmi_fn = mmi_fn,
                                  bed_fn = bed_fn,
                                  n_core = n_core)

    }

    dir.create(demult_dir, recursive = TRUE, showWarnings = FALSE)
    demult_fn <- file.path(demult_dir, "assignments.tsv")
    if(resume && file.exists(demult_fn)){
        message("Demultiplexing has already completed")
        demult_out <- utils::read.delim(demult_fn, stringsAsFactors = FALSE)

    } else {
        # New edlib demux expects user FASTQ; legacy basecall FASTA chunks are accepted for transition.
        demult_res <- doDemultiplex(fastq = basecall_fn,
                                    demult_dir = demult_dir,
                                    index_list = index_list,
                                    sample_list = sample_list,
                                    n_core = n_core,
                                    split_reads = FALSE)
        demult_out <- demult_res$assignments
    }

    dir.create(align_dir, recursive = TRUE, showWarnings = FALSE)
    align_fn <- file.path(align_dir, "alignment_list.csv")
    if(resume && file.exists(align_fn)){
        message("Alignment has already completed")
        align_out <- read.csv(align_fn)
        intact_seq <- readDNAStringSet(file.path(align_dir, "intact_seq.fa"))
        attributes(align_out) <- c(attributes(align_out), list(intact_seq = intact_seq))

    } else {
        align_out <- doAlign(blast_path = blast_path,
                             basecall_fn = basecall_fn,
                             amplicon_fn = amplicon_fn,
                             align_dir = align_dir,
                             primer_list = primer_list,
                             pam_list = pam_list,
                             genome_fn = genome_fn,
                             check_window = check_window,
                             n_core = n_core)
    }

    dir.create(editcall_dir, recursive = TRUE, showWarnings = FALSE)
    editcall_fn <- file.path(editcall_dir, "editcall_summary.csv")
    if(resume && file.exists(editcall_fn)){
        message("Editcalling has already completed")
        editcall_out <- read.csv(editcall_fn)

    } else {
        editcall_out <- doEditcall(demult_out = demult_out,
                                   align_out = align_out,
                                   editcall_dir = editcall_dir,
                                   sample_list = sample_list,
                                   strict = strict)
    }
    return(editcall_out)
}
