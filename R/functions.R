################################################################################
# Prepare indexed amplicon database
################################################################################
#' Prepare amplicon database for read alignment
#'
#' This function prepares an amplicon database for read alignment
#' by extracting amplicon sequences from a reference genome based
#' on provided primer sequences. It uses BLAST to identify the
#' locations of the primers in the reference genome and
#' extracts the corresponding amplicon sequences.
#'
#' @param blast_path Path to a directory containing BLAST executables (makeblastdb and blastn).
#' @param primer_list Path to a CSV file containing primer sequences.
#' The file should have two columns: primer ID and primer sequence.
#' Forward primers should have IDs ending with "_F" and reverse primers
#' with IDs ending with "_R".
#' @param genome_fn Path to the reference genome sequence in FASTA format.
#' @param out_dir Output directory where the amplicon database will be created.
#' @param n_core Number of CPU cores to use for BLAST.
#'
#' @export
#'
#' @import Biostrings
#' @import dplyr
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom BiocGenerics width
#' @importFrom Biostrings readDNAStringSet writeXStringSet
#' @importFrom stats setNames
#' @importFrom utils read.csv write.csv
#' @importFrom tools file_ext
#' @importFrom parallel mclapply
#' @importFrom methods as
#'
#' @examples
#' # blast_path <- "/path/to/blast/bin" # Path to a directory containing blast executables (makeblastdb and blastn)
#' # primer_list <- "inst/extdata/amplicon_primers.csv"
#' # genome_fn <- '/reference/genome/sequence.fa' # Reference genome sequence in fasta format
#' # out_dir <- "/path/to/your/working/directory/output_directory_name"
#' # n_core <- 30 # Number of CPU cores to use
#' # amplicon_fn <- prepAmpliconDB(blast_path = blast_path,
#' #                               primer_list = primer_list,
#' #                               genome_fn = genome_fn,
#' #                               out_dir = out_dir,
#' #                               n_core = n_core)
#'
#' @return Path to the generated amplicon FASTA file.
#'
prepAmpliconDB <- function(blast_path,
                           primer_list,
                           genome_fn,
                           out_dir,
                           n_core){
    primers <- read.csv(primer_list, header = FALSE)
    f_primers <- DNAStringSet(primers$V2[grep("_F$", primers$V1)])
    names(f_primers) <- primers$V1[grep("_F$", primers$V1)]
    r_primers <- DNAStringSet(primers$V2[grep("_R$", primers$V1)])
    names(r_primers) <- primers$V1[grep("_R$", primers$V1)]

    ref_dir <- file.path(out_dir, "ref")
    dir.create(ref_dir, recursive = TRUE, showWarnings = FALSE)
    f_primer_fn <- file.path(ref_dir, "f_primer.fa")
    writeXStringSet(f_primers, f_primer_fn)
    r_primer_fn <- file.path(ref_dir, "r_primer.fa")
    writeXStringSet(r_primers, r_primer_fn)
    db_path <- file.path(ref_dir, "ref_genome.fa")
    file.copy(genome_fn, db_path)
    .makeblastdb(blast_path = blast_path, db_path = db_path)

    blastout_fn <- sub("\\.fa", ".blastout", f_primer_fn)
    f_out <- .run_blastn(blast_path = blast_path,
                         query_fn = f_primer_fn,
                         db_path = db_path,
                         blastout_fn = blastout_fn,
                         task = "blastn-short",
                         outfmt = "short",
                         word_size = 4,
                         n_core = n_core)
    f_out <- .parse_blastout(blastout_fn = f_out, outfmt = "short")
    f_out <- tapply(seq_along(f_out$qseqid), f_out$qseqid, function(i){
        out_i <- f_out[i, ]
        hit <- out_i$length == width(f_primers[names(f_primers) %in% out_i$qseqid]) & out_i$pident >= 100
        out_i <- out_i[hit, ]
        return(out_i)
    })
    f_out <- do.call("rbind", f_out)

    blastout_fn <- sub("\\.fa", ".blastout", r_primer_fn)
    r_out <- .run_blastn(blast_path = blast_path,
                         query_fn = r_primer_fn,
                         db_path = db_path,
                         blastout_fn = blastout_fn,
                         task = "blastn-short",
                         outfmt = "short",
                         word_size = 4,
                         n_core = n_core)
    r_out <- .parse_blastout(blastout_fn = r_out, outfmt = "short")
    r_out <- tapply(seq_along(r_out$qseqid), r_out$qseqid, function(i){
        out_i <- r_out[i, ]
        hit <- out_i$length == width(r_primers[names(r_primers) %in% out_i$qseqid]) & out_i$pident >= 100
        out_i <- out_i[hit, ]
    })
    r_out <- do.call("rbind", r_out)

    f_out$qseqid <- sub("_F", "", f_out$qseqid)
    r_out$qseqid <- sub("_R", "", r_out$qseqid)
    amplicon_df <- full_join(f_out, r_out, c("qseqid", "sseqid"))
    amplicon_df <- subset(amplicon_df,
                          select = c(qseqid,
                                     sseqid,
                                     sstart.x,
                                     send.x,
                                     sstart.y,
                                     send.y))
    amplicon_df$start <- apply(amplicon_df[, -(1:2)], 1, min)
    amplicon_df$end <- apply(amplicon_df[, -(1:2)], 1, max)
    amplicon_gr <- GRanges(seqnames = amplicon_df$sseqid,
                           range = IRanges(start = amplicon_df$start,
                                           end = amplicon_df$end),
                           id = amplicon_df$qseqid)

    ref_genome <- readDNAStringSet(db_path)
    amplicon_seq <- ref_genome[amplicon_gr]
    names(amplicon_seq) <- amplicon_gr$id
    amplicon_fn <- file.path(ref_dir, "amplicon.fa")
    writeXStringSet(amplicon_seq, amplicon_fn)
    return(amplicon_fn)
}

.makeblastdb <- function(blast_path, db_path){
    blast_args <- paste(paste("-in", db_path),
                        paste("-dbtype nucl"))
    check <- try({
        system2(command = file.path(blast_path, "makeblastdb"),
                args = blast_args,
                stdout = FALSE)
    }, silent = TRUE)
    if(inherits(check, "try-error")){
        stop("Error in makeblastdb of BLAST.\n", check)
    }
}

.run_blastn <- function(blast_path,
                        query_fn,
                        db_path,
                        blastout_fn,
                        task = "blastn-short",
                        outfmt = "short",
                        word_size = 4,
                        n_core){
    if(outfmt == "short"){
        outfmt <- "'6 qseqid sseqid qstart qend sstart send length pident'"

    } else {
        outfmt <- "'6 qseqid sseqid qstart qend qlen sstart send slen mismatch gapopen sseq qseq'"
    }
    blast_args <- paste(paste("-query", query_fn),
                        paste("-db", db_path),
                        paste("-task", task),
                        "-max_target_seqs 1000000",
                        "-max_hsps 1",
                        paste("-outfmt", outfmt),
                        paste("-num_threads", n_core),
                        paste("-out", blastout_fn))

    if(!is.null(word_size)){
        blast_args <- paste(blast_args,
                            paste("-word_size", word_size))
    }

    system2(command = file.path(blast_path, "blastn"),
            args = blast_args,
            stdout = FALSE)
    return(blastout_fn)
}

.parse_blastout <- function(blastout_fn, outfmt){
    if(outfmt == "short"){
        out_names <- c("qseqid", "sseqid", "qstart", "qend",
                       "sstart", "send", "length", "pident")

    } else {
        out_names <- c("qseqid", "sseqid", "qstart", "qend",
                       "qlen", "sstart", "send", "slen", "mismatch",
                       "gapopen", "sseq", "qseq")
    }
    out <- read.table(blastout_fn, sep = "\t")
    names(out) <- out_names
    out <- subset(out, subset = !is.na(qseqid) & !is.na(sseqid))
    return(out)
}

################################################################################
# Main pipeline wrapper
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
#'
#' @export
#'
#' @import Biostrings
#' @import dplyr
#' @importFrom BiocGenerics width
#' @importFrom Biostrings readDNAStringSet writeXStringSet
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
                         amplicon_fn,
                         size_sel,
                         check_window = 10,
                         n_core = 1,
                         resume = FALSE){
    basecall_dir <- file.path(out_dir, "basecall")
    dir.create(basecall_dir, recursive = TRUE, showWarnings = FALSE)
    basecall_fn <- file.path(basecall_dir, "basecall.finish")
    if(resume & file.exists(basecall_fn)){
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

    demult_dir <- file.path(out_dir, "demultiplex")
    dir.create(demult_dir, recursive = TRUE, showWarnings = FALSE)
    demult_fn <- file.path(demult_dir, "demultiplex_list.csv")
    if(resume & file.exists(demult_fn)){
        message("Demultiplexing has already completed")
        demult_out <- read.csv(demult_fn)

    } else {
        demult_out <- doDemultiplex(blast_path = blast_path,
                                    basecall_fn = basecall_fn,
                                    demult_dir = demult_dir,
                                    index_list = index_list)
    }

    align_dir <- file.path(out_dir, "align")
    dir.create(align_dir, recursive = TRUE, showWarnings = FALSE)
    align_fn <- file.path(align_dir, "alignment_list.csv")
    if(resume & file.exists(align_fn)){
        message("Alignment has already completed")
        align_out <- read.csv(align_fn)

    } else {
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

    editcall_dir <- file.path(out_dir, "editcall")
    dir.create(editcall_dir, recursive = TRUE, showWarnings = FALSE)
    editcall_fn <- file.path(editcall_dir, "editcall_summary.csv")
    if(resume & file.exists(editcall_fn)){
        message("Editcalling has already completed")
        editcall_out <- read.csv(editcall_fn)

    } else {
        editcall_out <- doEditcall(demult_out = demult_out,
                                   align_out = align_out,
                                   editcall_dir = editcall_dir)
    }
    return(editcall_out)
}

.check_files <- function(in_dir,
                         out_dir,
                         dorado_path,
                         samtools_path,
                         primer_list,
                         pam_list,
                         index_list,
                         amplicon_fn){

}

################################################################################
# Basecalling
################################################################################
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
#' @param n_core Number of CPU cores to use for parallel processing.
#' @return A character vector containing paths to the size-selected FASTA files of basecalled reads.
#' @import Biostrings
#' @importFrom BiocGenerics width
#' @importFrom Biostrings readDNAStringSet writeXStringSet
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
                       n_core){
    bam_fn <- file.path(basecall_dir, "basecall.bam")
    system2(command = dorado_path,
            args = paste("duplex sup", in_dir,
                         "--threads", n_core,
                         "--min-qscore 10 > ", bam_fn
            ))

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

################################################################################
# Demultiplexing
################################################################################
#' Demultiplexing using BLAST
#'
#' This function performs demultiplexing of basecalled reads using BLAST.
#' It matches the reads against provided index sequences to assign them to their respective samples.
#' @param blast_path Path to a directory containing BLAST executables (makeblastdb and blastn).
#' @param basecall_fn A character vector containing paths to the FASTA files of basecalled reads.
#' @param demult_dir Path to the output directory where demultiplexing results will be saved.
#' @param index_list Path to a CSV file containing index sequences for demultiplexing.
#' The file should have five columns: index ID, forward index ID, forward index sequence, reverse index ID, and reverse index sequence.
#' @return A data frame summarizing the demultiplexing results.
#' @import Biostrings
#' @import dplyr
#' @importFrom BiocGenerics width
#' @importFrom Biostrings readDNAStringSet writeXStringSet
#' @importFrom utils read.csv write.csv
#' @importFrom stats setNames
#' @importFrom methods as
#' @export
#'
doDemultiplex <- function(blast_path, basecall_fn, demult_dir, index_list){
    index_df <- read.csv(index_list, header = FALSE)
    f_index <- index_df$V3
    f_index_id <- index_df$V2[!duplicated(index_df$V2)]
    f_index_uniq <- f_index[!duplicated(index_df$V2)]
    f_index_uniq <- DNAStringSet(f_index_uniq)
    names(f_index_uniq) <- f_index_id
    f_index_fn <- file.path(demult_dir, "f_index.fa")
    writeXStringSet(f_index_uniq, f_index_fn)
    r_index <- index_df$V5
    r_index_id <- index_df$V4[!duplicated(index_df$V4)]
    r_index_uniq <- r_index[!duplicated(index_df$V4)]
    r_index_uniq <- DNAStringSet(r_index_uniq)
    names(r_index_uniq) <- r_index_id
    r_index_fn <- file.path(demult_dir, "r_index.fa")
    writeXStringSet(r_index_uniq, r_index_fn)

    demult_out_fn <- NULL
    for(i in seq_along(basecall_fn)){
        i_demult_out_fn <- .loopBLAST(blast_path = blast_path,
                                      db_path = basecall_fn[i],
                                      f_index_fn = f_index_fn,
                                      r_index_fn = r_index_fn)
        demult_out_fn <- c(demult_out_fn, list(i_demult_out_fn))
    }

    demult_out <- NULL
    demult_rest <- NULL
    for(i in seq_along(demult_out_fn)){
        i_demult_out <- .loopParse(demult_out_fn = demult_out_fn[[i]])
        demult_out <- rbind(demult_out, i_demult_out$demult_out)
        demult_rest <- rbind(demult_rest, i_demult_out$fr_rest)
    }

    demult_out$index_pair <- paste(demult_out$qseqid.f, demult_out$qseqid.r, sep = "_")
    index_pair <- paste(index_df$V2, index_df$V4, sep = "_")
    index_pair <- gsub("\\s", "_", index_pair)
    index_pair_hit <- match(demult_out$index_pair, index_pair)
    demult_out$index_pair_id <- index_df$V1[index_pair_hit]

    demult_fn <- file.path(demult_dir, "demultiplex_list.csv")
    write.csv(demult_out, demult_fn, row.names = FALSE)
    demult_fn <- file.path(demult_dir, "undemultiplex_list.csv")
    write.csv(demult_rest, demult_fn, row.names = FALSE)
    return(demult_out)
}

.loopBLAST <- function(blast_path, db_path, f_index_fn, r_index_fn){
    .makeblastdb(blast_path = blast_path, db_path = db_path)

    blastout_prefix <- sub("\\/basecall\\/", "/demultiplex/", sub("\\.fq", "", db_path))
    blastout_fn <- paste(blastout_prefix,
                         sub("\\.fa", ".blastout", basename(f_index_fn)),
                         sep = "_")
    f_hit <- .run_blastn(blast_path = blast_path,
                         query_fn = f_index_fn,
                         db_path = db_path,
                         blastout_fn = blastout_fn,
                         task = "blastn-short",
                         outfmt = "long",
                         word_size = 4,
                         n_core = n_core)

    blastout_fn <- paste(blastout_prefix,
                         sub("\\.fa", ".blastout", basename(r_index_fn)),
                         sep = "_")
    r_hit <- .run_blastn(blast_path = blast_path,
                         query_fn = r_index_fn,
                         db_path = db_path,
                         blastout_fn = blastout_fn,
                         task = "blastn-short",
                         outfmt = "long",
                         word_size = 4,
                         n_core = n_core)

    out <- list(f_hit = f_hit, r_hit = r_hit)
    return(out)
}

.loopParse <- function(demult_out_fn){
    f_hit <- .parse_blastout(blastout_fn = demult_out_fn$f_hit, outfmt = "long")
    r_hit <- .parse_blastout(blastout_fn = demult_out_fn$r_hit, outfmt = "long")
    f_hit$posdiff <- f_hit$slen - f_hit$sstart
    front <- f_hit$slen - f_hit$posdiff <= f_hit$slen / 2
    rear <- f_hit$posdiff < f_hit$slen / 2
    f_hit$barpos<- NA
    f_hit$barpos[front] <- "front"
    f_hit$barpos[rear] <- "rear"

    r_hit$posdiff <- r_hit$slen - r_hit$sstart
    front <- r_hit$slen - r_hit$posdiff <= r_hit$slen / 2
    rear <- r_hit$posdiff < r_hit$slen / 2
    r_hit$barpos<- NA
    r_hit$barpos[front] <- "front"
    r_hit$barpos[rear] <- "rear"
    fr_best <- full_join(f_hit, r_hit, by = "sseqid", suffix = c(".f", ".r"),
                         relationship = "many-to-many")
    fr_best <- subset(fr_best, subset = barpos.f != barpos.r)

    f_no_variant <- (fr_best$mismatch.f + fr_best$gapopen.f) == 0
    f_full_cov <- (fr_best$qend.f - fr_best$qstart.f + 1) == fr_best$qlen.f
    r_no_variant <- (fr_best$mismatch.r + fr_best$gapopen.r) == 0
    r_full_cov <- (fr_best$qend.r - fr_best$qstart.r + 1) == fr_best$qlen.r
    fr_compl_match <- f_no_variant & f_full_cov & r_no_variant & r_full_cov
    if(sum(fr_compl_match) > 0){
        fr_compl_match <- fr_best[fr_compl_match, ]
        fr_compl_match$class <- "complete_match"
        fr_best_rest <- subset(fr_best, subset = !sseqid %in% fr_compl_match$sseqid)

    } else {
        fr_compl_match <- NULL
    }


    fr_best_subset <- subset(fr_best_rest, select = c(qseqid.f, sseqid, qseqid.r))
    fr_best_subset <- unique(fr_best_subset)
    fr_single_match <- table(fr_best_subset$sseqid)
    fr_single_match <- fr_best_rest$sseqid %in% names(fr_single_match[fr_single_match == 1])
    if(sum(fr_single_match) > 0){
        fr_single_match <- fr_best_rest[fr_single_match, ]
        fr_single_match$class <- "single_match"
        fr_best_rest <- subset(fr_best_rest, subset = !sseqid %in% fr_single_match$sseqid)

    } else {
        fr_single_match <- NULL
    }

    fr_partial_match_variants <- (fr_best_rest$qlen.f - fr_best_rest$qend.f + fr_best_rest$qstart.f) +
        fr_best_rest$mismatch.f +
        fr_best_rest$gapopen.f +
        (fr_best_rest$qlen.r - fr_best_rest$qend.r + fr_best_rest$qstart.r) +
        fr_best_rest$mismatch.r +
        fr_best_rest$gapopen.r
    valid_i <- tapply(seq_along(fr_partial_match_variants), fr_best_rest$sseqid, function(i){
        x <- fr_partial_match_variants[i]
        min_x <- min(x)
        check <- sum(x == min_x)
        if(check == 1){
            return(i[x == min_x])

        } else {
            NULL
        }
    })
    valid_i <- unlist(valid_i)
    if(length(valid_i) > 0){
        fr_partial_match <- fr_best_rest[valid_i, ]
        fr_partial_match$class <- "partial_match"
        fr_rest <- subset(fr_best_rest, subset = !sseqid %in% fr_partial_match$sseqid)

    } else {
        fr_partial_match <- NULL
        fr_rest <- fr_best_rest
    }

    fr_rest <- tapply(seq_along(fr_rest$sseqid), fr_rest$sseqid, function(i){
        i_f_min <- min(fr_rest$qstart.f[i])
        i_r_min <- min(fr_rest$qstart.r[i])
        sel <- fr_rest$qstart.f[i] == i_f_min & fr_rest$qstart.r[i] == i_r_min
        sel <- i[which(sel)][1]
        return(fr_rest[sel, ])
    })
    fr_rest <- do.call("rbind", fr_rest)

    demult_out <- rbind(fr_single_match, fr_compl_match, fr_partial_match)
    return(list(demult_out = demult_out, fr_rest = fr_rest))
}

################################################################################
# align to amplicon
################################################################################
#' Align reads to amplicon sequences using BLAST
#' This function aligns basecalled reads to a provided amplicon database using BLAST.
#' It processes the basecalled reads, performs the alignment, and extracts relevant alignment information.
#' @param blast_path Path to a directory containing BLAST executables (makeblastdb and blastn).
#' @param basecall_fn A character vector containing paths to the FASTA files of basecalled reads.
#' @param demult_out A data frame summarizing the demultiplexing results.
#' @param amplicon_fn Path to the amplicon database FASTA file.
#' @param align_dir Path to the output directory where alignment results will be saved.
#' @param primer_list Path to a CSV file containing primer sequences.
#' @param pam_list Path to a CSV file containing PAM site information.
#' @param genome_fn Path to the reference genome sequence in FASTA format.
#' @param check_window An integer specifying the window size (in bp) around the expected cut site to search for edits.
#' @return A data frame summarizing the alignment results.
#' @import Biostrings
#' @import dplyr
#' @importFrom pwalign pairwiseAlignment aligned
#' @importFrom BiocGenerics width
#' @importFrom Biostrings readDNAStringSet writeXStringSet
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom utils read.csv write.csv
#' @importFrom stats setNames
#' @importFrom methods as
#' @export
#'
doAlign <- function(blast_path,
                    basecall_fn,
                    demult_out,
                    amplicon_fn,
                    align_dir,
                    primer_list,
                    pam_list,
                    genome_fn,
                    check_window = 10){

    blastout_suffix <- sub("\\.fa", ".blastout", basename(amplicon_fn))
    ampl_hit_fn <- NULL
    for(i in seq_along(basecall_fn)){
        blastout_fn <- paste(sub("\\/basecall\\/", "/align/",
                                 sub("\\.fq", "", basecall_fn[i])),
                             blastout_suffix,
                             sep = "_")
        i_ampl_hit_fn <- .run_blastn(blast_path = blast_path,
                                     query_fn = amplicon_fn,
                                     db_path = basecall_fn[i],
                                     blastout_fn = blastout_fn,
                                     task = "blastn",
                                     outfmt = "long",
                                     word_size = NULL,
                                     n_core = n_core)
        ampl_hit_fn <- c(ampl_hit_fn, i_ampl_hit_fn)
    }

    ampl_hit <- NULL
    for(i in seq_along(ampl_hit_fn)){
        i_ampl_hit <- .parse_blastout(blastout_fn = ampl_hit_fn[i],
                                      outfmt = "long")
        ampl_hit <- rbind(ampl_hit, i_ampl_hit)
    }

    n_id_hit <- table(ampl_hit$sseqid)
    single_hit <- names(n_id_hit)[n_id_hit == 1]
    multi_hit <- ampl_hit[!ampl_hit$sseqid %in% single_hit, ]
    single_hit <- ampl_hit[ampl_hit$sseqid %in% single_hit, ]
    multi_hit_check <- tapply(seq_along(multi_hit$sseqid), multi_hit$sseqid, function(i){
        q_start <- multi_hit$qstat[i]
        q_end <- multi_hit$qend[i]
        inside <- q_start >= q_start[1] & q_end <= tail(q_end, 1)
        return(all(inside))
    })
    multi_hit_check <- unlist(multi_hit_check)
    if(any(!multi_hit_check)){
        stop("Multiple hit reads")
    }

    ampl_hit <- ampl_hit[!duplicated(ampl_hit$sseqid), ]

    amplicon_seq <- readDNAStringSet(amplicon_fn)
    genome <- readDNAStringSet(genome_fn)
    edit_site <- read.csv(pam_list, header = FALSE)
    edit_site_gr <- GRanges(seqnames = paste0("chr", sprintf("%02d", edit_site$V2)),
                            ranges = IRanges(start = edit_site$V3 - check_window,
                                             end = edit_site$V3 + check_window),
                            target = edit_site$V1)
    genome_seq <- genome[edit_site_gr]
    names(genome_seq) <- edit_site_gr$target
    genome_seq <- genome_seq[order(names(genome_seq))]
    amplicon_seq <- amplicon_seq[order(names(amplicon_seq))]

    aln <- Map(f = function(x, y){
        pairwiseAlignment(x, y)
    }, genome_seq, amplicon_seq)

    aln <- lapply(aln, aligned)
    aln <- lapply(aln, as.character)
    aln <- lapply(aln, gregexpr, pattern = "-")
    aln <- lapply(aln, unlist)
    aln_diff <- lapply(aln, diff)
    aln_diff <- lapply(aln_diff, ">", 1)
    aln_diff <- lapply(aln_diff, which)
    aln_pos <- Map(f = function(x, y){
        c(y[x], y[x + 1])
    }, aln_diff, aln)
    aln_pos <- do.call("rbind", aln_pos)
    aln_pos <- as.data.frame(aln_pos)
    names(aln_pos) <- c("start", "end")
    aln_pos$start <- aln_pos$start + 1
    aln_pos$end <- aln_pos$end - 1

    primers <- read.csv(primer_list, header = FALSE)
    f_primers <- DNAStringSet(primers$V2[grep("_F$", primers$V1)])
    names(f_primers) <- primers$V1[grep("_F$", primers$V1)]
    r_primers <- DNAStringSet(primers$V2[grep("_R$", primers$V1)])
    names(r_primers) <- primers$V1[grep("_R$", primers$V1)]
    f_primers <- f_primers[order(names(f_primers))]
    r_primers <- r_primers[order(names(r_primers))]
    aln_pos$f_primer_len <- width(f_primers)
    aln_pos$r_primer_len <- width(r_primers)
    aln_pos$dist_to_end <- width(amplicon_seq) - aln_pos$end
    aln_pos$target_len <- aln_pos$end - aln_pos$start

    target <- names(amplicon_seq)
    aln_pos$target <- target
    align_out <- NULL
    for(i in seq_along(target)){
        target_gene <- target[i]
        i_hit <- subset(ampl_hit, subset = qseqid == target_gene)
        i_aln_pos <- subset(aln_pos, subset = target == target_gene)
        s_cover_edit_site <- i_hit$qstart <= i_aln_pos$start
        e_cover_edit_site <- i_hit$qend >= i_aln_pos$start
        cover_edit_site <- s_cover_edit_site & e_cover_edit_site
        i_hit <- i_hit[cover_edit_site, ]
        i_target <- lapply(seq_along(i_hit$qseq), function(j){
            x <- i_hit$qseq[j]
            ins <- unlist(gregexec("-", substr(x, 1, i_aln_pos$start - 1)))
            if(ins[1] < 0){
                n_ins <- 0

            } else {
                n_ins <- length(ins)
            }
            target_start <- i_aln_pos$start + n_ins
            target_rest <- substr(x, target_start, nchar(x))
            target_end <- i_aln_pos$target_len
            detected_ins <- 0
            while(TRUE){
                ins <- unlist(gregexec("-", substr(target_rest, 1, target_end)))
                if(ins[1] < 0){
                    n_ins <- 0

                } else {
                    n_ins <- length(ins)
                }
                if(n_ins > detected_ins){
                    added_ins <- n_ins - detected_ins
                    target_end <- target_end + added_ins
                    detected_ins <- n_ins
                } else {
                    break
                }
            }
            target_end <- target_start + target_end
            j_target <- substr(i_hit$sseq[j], target_start, target_end)
            j_query <- substr(i_hit$qseq[j], target_start, target_end)
            return(c(j_target, j_query))
        })
        i_target <- do.call("rbind", i_target)
        i_out <- data.frame(target_gene = i_hit$qseqid,
                            read_name = i_hit$sseqid,
                            read_seq = i_target[, 1],
                            ref_seq = i_target[, 2])
        align_out <- rbind(align_out, i_out)
    }

    intact_seq <- as.character(genome_seq)
    align_out$intact <- align_out$read_seq %in% intact_seq

    full_amplicon <- data.frame(read_name = ampl_hit$sseqid,
                                full_read_seq = ampl_hit$sseq,
                                aln_start = ampl_hit$sstart,
                                aln_end = ampl_hit$send)
    full_amplicon <- unique(full_amplicon)
    align_out <- left_join(align_out, full_amplicon, "read_name")

    align_fn <- file.path(align_dir, "alignment_list.csv")
    write.csv(align_out, align_fn, row.names = FALSE)
    attributes(align_out) <- c(attributes(align_out), list(intact_seq = intact_seq))
    return(align_out)
}

################################################################################
# Evaluate alignments
################################################################################
#' Evaluate alignments and perform edit-calling
#' This function evaluates the alignments of reads to amplicon sequences and performs edit-calling.
#' It summarizes the edit-calling results and saves them to output files.
#' @param demult_out A data frame summarizing the demultiplexing results.
#' @param align_out A data frame summarizing the alignment results.
#' @param editcall_dir Path to the output directory where edit-calling results will be saved.
#' @return A data frame summarizing the edit-calling results.
#' @import Biostrings
#' @import dplyr
#' @importFrom BiocGenerics width
#' @importFrom Biostrings readDNAStringSet writeXStringSet
#' @importFrom utils read.csv write.csv
#' @importFrom stats setNames
#' @importFrom methods as
#' @export
#'
doEditcall <- function(demult_out, align_out, editcall_dir){
    demult_df <- data.frame(read_name = demult_out$sseqid,
                            i7_index = demult_out$qseqid.f,
                            i5_index = demult_out$qseqid.r,
                            index_pair_id = demult_out$index_pair_id)
    demult_df <- unique(demult_df)
    align_out <- left_join(align_out, demult_df, "read_name")

    intact_seq <- attributes(align_out)$intact_seq
    edit_df <- tapply(seq_len(nrow(align_out)), align_out$index_pair_id, function(i){
        i_align_out <- align_out[i, ]
        i_out <- tapply(seq_len(nrow(i_align_out)), i_align_out$target_gene, function(j){
            ij_align_out <- i_align_out[j, ]
            intact_read <- table(ij_align_out$intact)
            read_seq_tbl <- table(ij_align_out$read_seq)
            ij_out <- data.frame(target_gene = ij_align_out$target_gene[1],
                                 read_seq = names(read_seq_tbl),
                                 count = as.numeric(read_seq_tbl))
            ij_out$intact <- ij_out$read_seq %in% intact_seq[names(intact_seq) %in% ij_align_out$target_gene[1]]
            ij_out <- ij_out[order(ij_out$count, decreasing = TRUE), ]
            return(ij_out)
        })
        i_out <- do.call("rbind", i_out)
        i_out <- cbind(index_pair_id = i_align_out$index_pair_id[1],
                       i_out)
        return(i_out)
    })
    edit_df <- do.call("rbind", edit_df)
    write.csv(edit_df, file.path(editcall_dir, "editcall_all.csv"), row.names = FALSE)

    edit_df_filtered <- tapply(seq_len(nrow(edit_df)), edit_df$index_pair_id, function(i){
        i_df <- edit_df[i, ]
        i_out <- tapply(seq_len(nrow(i_df)), i_df$target_gene, function(j){
            ij_df <- i_df[j, ]
            top_count <- max(ij_df$count)
            top_df <- ij_df[ij_df$count > top_count / 2, ]
            top_df$vs_intact_ratio <- 0
            top_df$intact_count <- 0
            if(any(top_df$intact)){
                top_df$intact_count <- top_df$count[top_df$intact]
            }
            top_df$vs_intact_ratio <- top_df$count / (top_df$count + top_df$intact_count)
            top_df$vs_intact_ratio[top_df$intact] <- NA
            return(top_df)
        })
        i_out <- do.call("rbind", i_out)
        return(i_out)
    })
    edit_df_filtered <- do.call("rbind", edit_df_filtered)
    write.csv(edit_df_filtered, file.path(editcall_dir, "editcall_filtered.csv"), row.names = FALSE)

    editcall_out <- tapply(seq_len(nrow(edit_df_filtered)), edit_df_filtered$index_pair_id, function(i){
        i_df <- edit_df_filtered[i, ]
        out <- data.frame(target_gene = names(intact_seq))
        i_out <- tapply(seq_len(nrow(i_df)), i_df$target_gene, function(j){
            ij_df <- i_df[j, ]
            ij_target_gene <- ij_df$target_gene[1]
            genotype <- rep("ref", nrow(ij_df))
            ij_intact_seq <- intact_seq[names(intact_seq) %in% ij_target_gene]
            ij_del <- grepl("-", ij_df$read_seq)
            ij_ins <- nchar(gsub("-", "", ij_df$read_seq)) > nchar(ij_intact_seq)
            genotype[!ij_df$intact] <- "sub"
            genotype[ij_del & ij_ins] <- "indel"
            genotype[ij_del & !ij_ins] <- "del"
            genotype[!ij_del & ij_ins] <- "ins"
            genotype_order <- factor(genotype, levels = c("ref", "sub", "ins", "del", "indel"))
            genotype_order <- order(as.numeric(genotype_order))
            if(any(ij_del & !ij_ins)){
                n_ij_del <- sapply(gregexec("-", ij_df$read_seq[ij_del & !ij_ins]), length)
                genotype[ij_del & !ij_ins] <- paste0(genotype[ij_del & !ij_ins], n_ij_del)
            }
            if(any(!ij_del & ij_ins)){
                n_ij_ins <- nchar(gsub("-", "", ij_df$read_seq[!ij_del & ij_ins])) - nchar(ij_intact_seq)
                genotype[!ij_del & ij_ins] <- paste0(genotype[!ij_del & ij_ins], n_ij_ins)
            }
            if(any(ij_del & ij_ins)){
                n_ij_del <- sapply(gregexec("-", ij_df$read_seq[ij_del & ij_ins]), length)
                n_ij_ins <- nchar(gsub("-", "", ij_df$read_seq[ij_del & ij_ins])) - nchar(ij_intact_seq)
                genotype[ij_del & ij_ins] <- paste0(genotype[ij_del & ij_ins], n_ij_ins, "-", n_ij_ins)
            }
            genotype <- genotype[genotype_order]
            which_ref <- genotype == "ref"
            genotype <- paste(genotype, collapse = "/")
            alt_patterns <- ij_df$read_seq[genotype_order]
            alt_patterns <- paste(alt_patterns[!which_ref], collapse = "/")
            count <- paste0(sum(ij_df$count),
                            "(",
                            paste(ij_df$count[genotype_order], collapse = ","),
                            ")")
            ij_out <- data.frame(target_gene = ij_target_gene,
                                 genotype = genotype,
                                 alt_patterns = alt_patterns,
                                 count = count)
            return(ij_out)
        })
        i_out <- do.call("rbind", i_out)
        out <- left_join(out, i_out, "target_gene")
        out <- t(out)
        colnames(out) <- out[1, ]
        out <- out[-1, ]
        out <- cbind(index_pair_id = i_df$index_pair_id[1],
                     sample_id = i_df$sample_id[1],
                     data_type = c("genotype", "seq", "count"),
                     out)
        return(out)
    })
    editcall_out <- do.call("rbind", editcall_out)
    editcall_fn <- file.path(editcall_dir, "editcall_summary.csv")
    write.csv(editcall_out, editcall_fn, row.names = FALSE)
    return(editcall_out)
}

################################################################################
# Evaluate MIAO results
################################################################################
#' @title Evaluate MIAO results
#' @description This function evaluates the results of the MIAO pipeline and generates a summary report.
#' @param out_dir The output directory where MIAO results are stored.
#' @param output_reads Logical. If TRUE, outputs the sequences of non-indexed aligned reads to a FASTA file.
#' @return A summary report saved in the output directory.
#'
#' @importFrom Biostrings readDNAStringSet writeXStringSet
#' @importFrom dplyr left_join
#'
#' @export
#'
evalMiao <- function(out_dir, output_reads){
    basecall_dir <- file.path(out_dir, "basecall")
    basecall_tsv <- file.path(basecall_dir, "basecalls_summary.tsv")
    n_raw_read <- length(count.fields(file = basecall_tsv))

    basecall_fn <- list.files(path = basecall_dir,
                              pattern = "basecall_filt_sizeselected_reads_.+.fa$",
                              full.names = TRUE)
    for(i in seq_along(basecall_fn)){
        if(i == 1){
            basecall_out <- readDNAStringSet(basecall_fn[i])
        } else {
            basecall_out <- c(basecall_out,
                              readDNAStringSet(basecall_fn[i]))
        }
    }
    n_filt_read <- length(basecall_out)

    demult_dir <- file.path(out_dir, "demultiplex")
    demult_fn <- file.path(demult_dir, "demultiplex_list.csv")
    demult_out <- read.csv(file = demult_fn)
    n_demult_reads <- length(unique(demult_out$sseqid))
    n_dup_index_reads <- sum(duplicated(demult_out$sseqid))

    demult_dir <- file.path(out_dir, "demultiplex")
    undemult_fn <- file.path(demult_dir, "undemultiplex_list.csv")
    undemult_out <- read.csv(file = undemult_fn)

    n_undemult <- nrow(undemult_out)
    n_indexed_both_side <- sum(undemult_out$qstart.f <= 15 & undemult_out$qstart.r <= 15, na.rm = TRUE)
    n_indexed_one_side <- sum(undemult_out$qstart.f <= 15 | undemult_out$qstart.r <= 15, na.rm = TRUE)
    n_not_indexed <- sum(undemult_out$qstart.f > 15 & undemult_out$qstart.r > 15, na.rm = TRUE)

    n_demult_reads <- length(unique(demult_out$sseqid))
    n_dup_index_reads <- sum(duplicated(demult_out$sseqid))

    n_indexed_reads_per_index_f <- table(demult_out$qseqid.f)
    prop_indexed_reads_per_index_f <- n_indexed_reads_per_index_f / sum(n_indexed_reads_per_index_f)
    n_indexed_reads_per_index_r <- table(demult_out$qseqid.r)
    prop_indexed_reads_per_index_r <- n_indexed_reads_per_index_r / sum(n_indexed_reads_per_index_r)

    align_dir <- file.path(out_dir, "align")
    align_fn <- file.path(align_dir, "alignment_list.csv")
    align_out <- read.csv(file = align_fn)
    n_align_reads <- length(unique(align_out$read_name))
    n_dup_align_reads <- sum(duplicated(align_out$read_name))
    n_align_reads_per_gene <- table(align_out$target_gene)
    prop_align_reads_per_gene <- signif(n_align_reads_per_gene / sum(n_align_reads_per_gene), 3)

    editcall_dir <- file.path(out_dir, "editcall")
    editcall_fn <- file.path(editcall_dir, "editcall_filtered.csv")
    editcall_out <- read.csv(file = editcall_fn)
    n_edicall_reads <- sum(editcall_out$count)
    n_edicall_reads_per_gene <- tapply(editcall_out$count, editcall_out$target_gene, sum)
    prop_edicall_reads_per_gene <- n_edicall_reads_per_gene / n_align_reads_per_gene

    non_indexed_aligned_reads <- undemult_out$sseqid[undemult_out$qstart.f > 15 & undemult_out$qstart.r > 15]
    non_indexed_aligned_reads <- align_out$read_name %in% non_indexed_aligned_reads
    n_non_indexed_aligned_reads_per_gene <- table(align_out$target_gene[non_indexed_aligned_reads])
    prop_non_indexed_aligned_reads_per_gene <- n_non_indexed_aligned_reads_per_gene / n_align_reads_per_gene

    if(output_reads){
        non_indexed_aligned_read_names <- align_out$read_name[non_indexed_aligned_reads]
        non_indexed_aligned_read_seq <- basecall_out[names(basecall_out) %in% non_indexed_aligned_read_names]
        writeXStringSet(non_indexed_aligned_read_seq, file.path(out_dir, "non_indexed_aligned_read_seq.fa"))
    }

    summary_dir <- file.path(out_dir, "miao_summary")
    dir.create(summary_dir, showWarnings = FALSE, recursive = TRUE)

    out1 <- rbind(c("Raw reads: ", n_raw_read, ""),
                  c("Reads after filtering: ",
                    n_filt_read,
                    signif(n_filt_read / n_raw_read, 3) * 100),
                  c("Demultiplexed reads: ",
                    n_demult_reads,
                    signif(n_demult_reads / n_raw_read, 3) * 100),
                  c("Aligned reads: ",
                    n_align_reads,
                    signif(n_align_reads / n_raw_read, 3) * 100),
                  c("Editcalled reads: ",
                    n_edicall_reads,
                    signif(n_edicall_reads / n_raw_read, 3) * 100),
                  c("Undemultiplexed reads: ",
                    n_undemult,
                    signif(n_undemult / n_raw_read, 3) * 100),
                  c("Undemultiplexed reads (umbiguously indexed): ",
                    n_indexed_both_side,
                    signif(n_indexed_both_side / n_raw_read, 3) * 100),
                  c("Undemultiplexed reads (single index): ",
                    n_indexed_one_side,
                    signif(n_indexed_one_side / n_raw_read, 3) * 100),
                  c("Undemultiplexed reads (no index): ",
                    n_not_indexed,
                    signif(n_not_indexed / n_raw_read, 3) * 100))
    write.table(x = out1, file = file.path(summary_dir, "read_stats.tsv"),
                row.names = FALSE, col.names = FALSE, sep = "\t")

    out2 <- rbind(names(n_indexed_reads_per_index_f),
                  n_indexed_reads_per_index_f,
                  signif(prop_indexed_reads_per_index_f, 3) * 100,
                  names(n_indexed_reads_per_index_f),
                  n_indexed_reads_per_index_r,
                  signif(prop_indexed_reads_per_index_r, 3) * 100)
    out2 <- cbind(c("",
                    "Indexed reads per forward index",
                    "Proportion of indexed reads per forward index",
                    "",
                    "Indexed reads per reverse index",
                    "Proportion of indexed reads per reverse index"),
                  out2)
    write.table(x = out2, file = file.path(summary_dir, "indexed_reads_per_gene.tsv"),
                row.names = FALSE, col.names = FALSE, sep = "\t")

    out3 <- rbind(names(n_align_reads_per_gene),
                  n_align_reads_per_gene,
                  signif(prop_align_reads_per_gene, 3) * 100,
                  n_edicall_reads_per_gene,
                  signif(prop_edicall_reads_per_gene, 3) * 100,
                  n_non_indexed_aligned_reads_per_gene,
                  signif(prop_non_indexed_aligned_reads_per_gene, 3) * 100)
    out3 <- cbind(c("",
                    "Aligned reads per gene",
                    "Proportion of aligned reads per gene",
                    "Editcalled reads per gene",
                    "Proportion of editcalled reads per gene",
                    "Non-indexed reads per gene",
                    "Proportion of non-indexed reads per gene"),
                  out3)
    write.table(x = out3, file = file.path(summary_dir, "aligned_reads_per_gene.tsv"),
                row.names = FALSE, col.names = FALSE, sep = "\t")
}
