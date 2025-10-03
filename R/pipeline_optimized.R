################################################################################
# Optimized Main Pipeline Wrapper
################################################################################

#' Optimized main pipeline wrapper for MiaoEditcall
#'
#' This function serves as the optimized main pipeline wrapper for MiaoEditcall,
#' using Rcpp-optimized functions for improved performance.
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
#' @param use_optimized Logical indicating whether to use optimized functions (default: TRUE).
#'
#' @export
#'
#' @import Biostrings
#' @import dplyr
#' @import Rcpp
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
#' # Load optimized functions
#' source("R/optimized_functions.R")
#' 
#' # Run optimized pipeline
#' editcall_out <- miaoEditcall_optimized(in_dir = in_dir,
#'                                      out_dir = out_dir,
#'                                      dorado_path = dorado_path,
#'                                      samtools_path = samtools_path,
#'                                      blast_path = blast_path,
#'                                      primer_list = primer_list,
#'                                      pam_list = pam_list,
#'                                      index_list = index_list,
#'                                      genome_fn = genome_fn,
#'                                      amplicon_fn = amplicon_fn,
#'                                      size_sel = c(300, 450),
#'                                      check_window = 10,
#'                                      n_core = n_core,
#'                                      resume = FALSE,
#'                                      use_optimized = TRUE)
#' @return A data frame summarizing the edit-calling results.
#'
miaoEditcall_optimized <- function(in_dir,
                                 out_dir,
                                 dorado_path = "dorado",
                                 samtools_path = "samtools",
                                 blast_path = "blastn",
                                 primer_list,
                                 pam_list,
                                 index_list,
                                 genome_fn,
                                 amplicon_fn,
                                 size_sel,
                                 check_window = 10,
                                 n_core = 1,
                                 resume = FALSE,
                                 use_optimized = TRUE) {
    
    # Load optimized functions if requested
    if (use_optimized) {
        source("R/optimized_functions.R")
    }
    
    # Basecalling step (unchanged - external tool)
    basecall_dir <- file.path(out_dir, "basecall")
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
                                  n_core = n_core)
    }
    
    # Demultiplexing step (unchanged - external tool)
    demult_dir <- file.path(out_dir, "demultiplex")
    dir.create(demult_dir, recursive = TRUE, showWarnings = FALSE)
    demult_fn <- file.path(demult_dir, "demultiplex_list.csv")
    if(resume && file.exists(demult_fn)){
        message("Demultiplexing has already completed")
        demult_out <- read.csv(demult_fn)
    } else {
        demult_out <- doDemultiplex(blast_path = blast_path,
                                    basecall_fn = basecall_fn,
                                    demult_dir = demult_dir,
                                    index_list = index_list)
    }
    
    # Alignment step - use optimized version if available
    align_dir <- file.path(out_dir, "align")
    dir.create(align_dir, recursive = TRUE, showWarnings = FALSE)
    align_fn <- file.path(align_dir, "alignment_list.csv")
    if(resume && file.exists(align_fn)){
        message("Alignment has already completed")
        align_out <- read.csv(align_fn)
    } else {
        if (use_optimized && exists("doAlign_optimized")) {
            message("Using optimized alignment function")
            align_out <- doAlign_optimized(blast_path = blast_path,
                                         basecall_fn = basecall_fn,
                                         demult_out = demult_out,
                                         amplicon_fn = amplicon_fn,
                                         align_dir = align_dir,
                                         primer_list = primer_list,
                                         pam_list = pam_list,
                                         genome_fn = genome_fn,
                                         check_window = check_window)
        } else {
            message("Using standard alignment function")
            align_out <- doAlign(blast_path = blast_path,
                                 basecall_fn = basecall_fn,
                                 demult_out = demult_out,
                                 amplicon_fn = amplicon_fn,
                                 align_dir = align_dir,
                                 primer_list = primer_list,
                                 pam_list = pam_list,
                                 genome_fn = genome_fn,
                                 check_window = check_window)
        }
    }
    
    # Edit-calling step - use optimized version if available
    editcall_dir <- file.path(out_dir, "editcall")
    dir.create(editcall_dir, recursive = TRUE, showWarnings = FALSE)
    editcall_fn <- file.path(editcall_dir, "editcall_summary.csv")
    if(resume && file.exists(editcall_fn)){
        message("Editcalling has already completed")
        editcall_out <- read.csv(editcall_fn)
    } else {
        if (use_optimized && exists("doEditcall_optimized")) {
            message("Using optimized edit-calling function")
            editcall_out <- doEditcall_optimized(demult_out = demult_out,
                                               align_out = align_out,
                                               editcall_dir = editcall_dir)
        } else {
            message("Using standard edit-calling function")
            editcall_out <- doEditcall(demult_out = demult_out,
                                       align_out = align_out,
                                       editcall_dir = editcall_dir)
        }
    }
    
    return(editcall_out)
}

################################################################################
# Performance Benchmarking Function
################################################################################

#' Benchmark performance comparison between standard and optimized functions
#'
#' @param test_data List containing test data for benchmarking
#' @param iterations Number of iterations to run for benchmarking
#' @return List with timing results
#'
#' @export
benchmark_miaoseq <- function(test_data, iterations = 10) {
    
    # Load both standard and optimized functions
    source("R/functions.R")
    source("R/optimized_functions.R")
    
    results <- list()
    
    # Benchmark string operations
    if ("sequences" %in% names(test_data)) {
        message("Benchmarking string operations...")
        
        # Pattern counting
        standard_time <- system.time({
            for (i in 1:iterations) {
                sapply(gregexpr("-", test_data$sequences), length)
            }
        })
        
        optimized_time <- system.time({
            for (i in 1:iterations) {
                fast_pattern_count_optimized(test_data$sequences, "-")
            }
        })
        
        results$pattern_counting <- list(
            standard = standard_time,
            optimized = optimized_time,
            speedup = standard_time[3] / optimized_time[3]
        )
        
        # Sequence length calculation
        standard_time <- system.time({
            for (i in 1:iterations) {
                nchar(gsub("-", "", test_data$sequences))
            }
        })
        
        optimized_time <- system.time({
            for (i in 1:iterations) {
                fast_seq_length_no_gaps_optimized(test_data$sequences)
            }
        })
        
        results$length_calculation <- list(
            standard = standard_time,
            optimized = optimized_time,
            speedup = standard_time[3] / optimized_time[3]
        )
    }
    
    # Benchmark sequence table operations
    if ("sequence_counts" %in% names(test_data)) {
        message("Benchmarking sequence table operations...")
        
        standard_time <- system.time({
            for (i in 1:iterations) {
                table(test_data$sequence_counts)
            }
        })
        
        optimized_time <- system.time({
            for (i in 1:iterations) {
                fast_sequence_table_optimized(test_data$sequence_counts, 
                                             rep(1, length(test_data$sequence_counts)))
            }
        })
        
        results$sequence_table <- list(
            standard = standard_time,
            optimized = optimized_time,
            speedup = standard_time[3] / optimized_time[3]
        )
    }
    
    return(results)
}

################################################################################
# Memory Usage Optimization
################################################################################

#' Optimize memory usage for large datasets
#'
#' @param data_size Estimated size of dataset
#' @param available_memory Available system memory in GB
#' @return List with recommended parameters
#'
#' @export
optimize_memory_usage <- function(data_size, available_memory = 16) {
    
    recommendations <- list()
    
    # Calculate optimal chunk size
    if (data_size > available_memory * 0.8) {
        recommendations$chunk_size <- floor(available_memory * 0.6 * 1024^3 / data_size)
        recommendations$use_chunking <- TRUE
    } else {
        recommendations$chunk_size <- data_size
        recommendations$use_chunking <- FALSE
    }
    
    # Calculate optimal number of cores
    if (available_memory < 8) {
        recommendations$max_cores <- 4
    } else if (available_memory < 16) {
        recommendations$max_cores <- 8
    } else if (available_memory < 32) {
        recommendations$max_cores <- 16
    } else {
        recommendations$max_cores <- 32
    }
    
    # Memory optimization flags
    recommendations$use_optimized <- TRUE
    recommendations$clear_intermediate <- TRUE
    recommendations$use_compression <- TRUE
    
    return(recommendations)
}
