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
#' @importFrom Biostrings readDNAStringSet writeXStringSet DNAStringSet
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
    amplicon_fn <- file.path(ref_dir, "amplicon.csv")
    write.csv(amplicon_df, amplicon_fn)
    amplicon_fn <- file.path(ref_dir, "amplicon.fa")
    writeXStringSet(amplicon_seq, amplicon_fn)
    return(amplicon_fn)
}
