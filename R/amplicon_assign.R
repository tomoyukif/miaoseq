################################################################################
# Step 3: Amplicon assignment
################################################################################
#' Align reads to amplicon sequences using BLAST
#'
#' This function assigns basecalled reads to amplicon targets using BLAST.
#' It processes the basecalled reads, performs the alignment, and extracts relevant alignment information.
#' @param blast_path Path to a directory containing BLAST executables (makeblastdb and blastn).
#' @param basecall_fn A character vector containing paths to the FASTA files of basecalled reads.
#' @param amplicon_fn Path to the amplicon database FASTA file.
#' @param align_dir Path to the output directory where alignment results will be saved.
#' @param primer_list Path to a CSV file containing primer sequences.
#' @param pam_list Path to a CSV file containing PAM site information.
#' @param genome_fn Path to the reference genome sequence in FASTA format.
#' @param check_window An integer specifying the window size (in bp) around the expected cut site to search for edits.
#' @param n_core Number of CPU cores to use for BLAST.
#' @return A data frame summarizing the alignment results.
#' @import Biostrings
#' @import dplyr
#' @importFrom pwalign pairwiseAlignment aligned
#' @importFrom BiocGenerics width
#' @importFrom Biostrings readDNAStringSet writeXStringSet DNAStringSet
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom utils read.csv write.csv
#' @importFrom stats setNames
#' @importFrom methods as
#' @export
#'
doAlign <- function(blast_path,
                    basecall_fn,
                    amplicon_fn,
                    align_dir,
                    primer_list,
                    pam_list,
                    genome_fn,
                    check_window = 10,
                    n_core = 1){

    blastout_suffix <- sub("\\.fa", ".blastout", basename(amplicon_fn))
    ampl_hit_fn <- NULL
    for(i in seq_along(basecall_fn)){
        blastout_fn <- paste(sub("\\/basecall\\/", "/align/",
                                 sub("\\.fq", "", basecall_fn[i])),
                             blastout_suffix,
                             sep = "_")
        if(file.exists(blastout_fn)){
            i_ampl_hit_fn <- blastout_fn

        } else {
            i_ampl_hit_fn <- .run_blastn(blast_path = blast_path,
                                         query_fn = amplicon_fn,
                                         db_path = basecall_fn[i],
                                         blastout_fn = blastout_fn,
                                         task = "blastn",
                                         outfmt = "long",
                                         word_size = NULL,
                                         n_core = n_core,
                                         not_run = FALSE)
        }
        ampl_hit_fn <- c(ampl_hit_fn, i_ampl_hit_fn)
    }

    ampl_hit <- NULL
    for(i in seq_along(ampl_hit_fn)){
        i_ampl_hit <- .parse_blastout(blastout_fn = ampl_hit_fn[i],
                                      outfmt = "long")
        ampl_hit <- rbind(ampl_hit, i_ampl_hit)
    }
    ampl_hit$index <- seq_len(nrow(ampl_hit))
    n_id_hit <- table(ampl_hit$sseqid)
    single_hit <- names(n_id_hit)[n_id_hit == 1]
    multi_hit <- subset(ampl_hit,
                        subset = !sseqid %in% single_hit,
                        select = c(sseqid, sstart, send, index))
    flip <- multi_hit$sstart > multi_hit$send
    tmp <- multi_hit$sstart[flip]
    multi_hit$sstart[flip] <- multi_hit$send[flip]
    multi_hit$send[flip] <- tmp

    invalid_multi_hit_reads <- tapply(seq_along(multi_hit$sseqid),
                                      multi_hit$sseqid,
                                      function(i){
                                          inside <- sapply(i, function(j){
                                              k <- setdiff(i, j)
                                              check_start <- multi_hit$sstart[j] >= multi_hit$sstart[k]
                                              check_end <- multi_hit$send[j] <= multi_hit$send[k]
                                              return(any(check_start & check_end))
                                          })
                                          if(sum(!inside) > 1){
                                              return(multi_hit$index[i])

                                          } else {
                                              len <- multi_hit$send[i] - multi_hit$sstart[i]
                                              max_len <- max(len)
                                              which_max <- len == max_len
                                              if(sum(which_max) > 1){
                                                  return(multi_hit$index[i])

                                              } else {
                                                  return(multi_hit$index[i[!which_max]])
                                              }
                                          }
                                      })

    ampl_hit <- subset(ampl_hit, subset = !index %in% unlist(invalid_multi_hit_reads))

    amplicon_seq <- readDNAStringSet(amplicon_fn)
    genome <- readDNAStringSet(genome_fn)
    edit_site <- read.csv(pam_list, header = FALSE)
    edit_site_gr <- GRanges(seqnames = paste0("chr", sprintf("%02d", edit_site$V2)),
                            ranges = IRanges(start = edit_site$V3 - check_window,
                                             end = edit_site$V3 + check_window),
                            target = edit_site$V1)
    intact_seq <- genome[edit_site_gr]
    names(intact_seq) <- edit_site_gr$target
    intact_seq <- intact_seq[order(names(intact_seq))]
    amplicon_seq <- amplicon_seq[order(names(amplicon_seq))]

    aln <- Map(f = function(x, y){
        pairwiseAlignment(x, y)
    }, intact_seq, amplicon_seq)

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
    if(!is.null(edit_site$V4)){
        aln_pos$guide_id <- edit_site$V4
    }
    align_out <- NULL
    for(i in seq_along(aln_pos$target)){
        target_gene <- aln_pos$target[i]
        i_hit <- subset(ampl_hit, subset = qseqid == target_gene)
        i_aln_pos <- aln_pos[i, ]
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

        if(!is.null(aln_pos$guide_id)){
            if(all(!aln_pos$guide_id[i] %in% c("", NA))){
                target_gene <- paste(target_gene, aln_pos$guide_id[i], sep = "_")
                names(intact_seq)[i] <- target_gene
            }
        }

        i_out <- data.frame(target_gene = target_gene,
                            read_name = i_hit$sseqid,
                            read_seq = i_target[, 1],
                            ref_seq = i_target[, 2])
        align_out <- rbind(align_out, i_out)
    }

    align_out$intact <- align_out$read_seq %in% intact_seq

    full_amplicon <- data.frame(read_name = ampl_hit$sseqid,
                                full_read_seq = ampl_hit$sseq,
                                aln_start = ampl_hit$sstart,
                                aln_end = ampl_hit$send,
                                query_start = ampl_hit$qstart,
                                query_end = ampl_hit$qend,
                                query_length = ampl_hit$qlen)
    full_amplicon$query_align <- full_amplicon$query_end - full_amplicon$query_start + 1
    full_amplicon <- unique(full_amplicon)
    align_out <- left_join(align_out, full_amplicon, "read_name")

    align_fn <- file.path(align_dir, "alignment_list.csv")
    write.csv(align_out, align_fn, row.names = FALSE)
    writeXStringSet(intact_seq, file.path(align_dir, "intact_seq.fa"))
    attributes(align_out) <- c(attributes(align_out), list(intact_seq = intact_seq))
    return(align_out)
}
