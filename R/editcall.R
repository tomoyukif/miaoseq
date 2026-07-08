################################################################################
# Step 4: Edit-calling
################################################################################
#' Evaluate alignments and perform edit-calling
#'
#' This function evaluates the alignments of reads to amplicon sequences and performs edit-calling.
#' It summarizes the edit-calling results and saves them to output files.
#' @param demult_out A data frame summarizing the demultiplexing results.
#' @param align_out A data frame summarizing the alignment results.
#' @param editcall_dir Path to the output directory where edit-calling results will be saved.
#' @param sample_list Optional path to a CSV file mapping index pairs to sample information.
#'   If provided, adds sample names and total read counts to the summary.
#' @param strict Logical. If TRUE, apply strict filtering on barcode and alignment positions.
#' @return A data frame summarizing the edit-calling results.
#' @import Biostrings
#' @import dplyr
#' @importFrom BiocGenerics width
#' @importFrom Biostrings readDNAStringSet writeXStringSet DNAStringSet
#' @importFrom utils read.csv write.csv
#' @importFrom stats setNames
#' @importFrom methods as
#' @export
#'
doEditcall <- function(demult_out, align_out, editcall_dir, sample_list = NULL, strict = FALSE){
    demult_out$read_name <- demult_out$sseqid
    align_out <- left_join(align_out,
                           subset(demult_out,
                                  select = c(read_name, qstart.f:qlen.f,
                                             qstart.r:qlen.r, send.f, send.r,
                                             qseqid.f, qseqid.r, index_pair_id)),
                           "read_name")
    align_out <- align_out[!is.na(align_out$target_gene), ]
    align_out <- align_out[!is.na(align_out$index_pair_id), ]

    if(strict){
        check1 <- align_out$qstart.f == 1 & align_out$qend.f == align_out$qlen.f
        check2 <- align_out$qstart.r == 1 & align_out$qend.r == align_out$qlen.r
        check3 <- align_out$query_align == align_out$query_length

        check4 <- align_out$send.f + 1 == align_out$aln_start
        check5 <- align_out$send.f - 1 == align_out$aln_start
        check6 <- align_out$send.r - 1 == align_out$aln_start
        check7 <- align_out$send.r + 1 == align_out$aln_start
        check8 <- check4 | check5 | check6 | check7

        check9 <- align_out$send.f + 1 == align_out$aln_end
        check10 <- align_out$send.f - 1 == align_out$aln_end
        check11 <- align_out$send.r - 1 == align_out$aln_end
        check12 <- align_out$send.r + 1 == align_out$aln_end
        check13 <- check9| check10 | check11 | check12
        valid <- check1 & check2 & check3 & check8 & check13
        align_out <- align_out[valid, ]
    }

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

    edit_df <- subset(edit_df, subset = count > 4)
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
                     data_type = c("genotype", "seq", "count"),
                     out)
        return(out)
    })
    editcall_out <- do.call("rbind", editcall_out)
    editcall_out <- as.data.frame(editcall_out)

    total_reads_per_sample <- tapply(edit_df_filtered$count, edit_df_filtered$index_pair_id, sum)
    total_reads_df <- data.frame(index_pair_id = names(total_reads_per_sample),
                                 total_reads = as.numeric(total_reads_per_sample))

    if(!is.null(sample_list)){
        sample_info <- read.csv(sample_list, header = FALSE)
        names(sample_info) <- c("index_pair_id", "sample_name", "plate_id", "row_id", "col_id")

        editcall_out <- left_join(editcall_out, sample_info, by = "index_pair_id")
        editcall_out <- left_join(editcall_out, total_reads_df, by = "index_pair_id")

        gene_cols <- setdiff(names(editcall_out), c("index_pair_id", "data_type",
                                                    "sample_name", "plate_id", "row_id", "col_id", "total_reads"))
        editcall_out <- editcall_out[, c("index_pair_id", "data_type", gene_cols,
                                         "sample_name", "plate_id", "row_id", "col_id", "total_reads")]
    } else {
        editcall_out <- left_join(editcall_out, total_reads_df, by = "index_pair_id")
    }

    editcall_fn <- file.path(editcall_dir, "editcall_summary.csv")
    write.csv(editcall_out, editcall_fn, row.names = FALSE)
    return(editcall_out)
}
