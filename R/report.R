################################################################################
# Reporting and visualization
################################################################################

#' @title Evaluate MIAO results
#' @description This function evaluates the results of the MIAO pipeline and generates a summary report.
#' @param out_dir The output directory where MIAO results are stored.
#' @param output_reads Logical. If TRUE, outputs the sequences of non-indexed aligned reads to a FASTA file.
#' @return A summary report saved in the output directory.
#'
#' @importFrom Biostrings readDNAStringSet writeXStringSet DNAStringSet
#' @importFrom dplyr left_join full_join
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
    n_edicall_reads_per_gene <- tapply(editcall_out$count,
                                       editcall_out$target_gene,
                                       sum)
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

#' Visualize edit-calling results per plate as PDF heatmaps
#'
#' This function takes the edit-calling summary and produces per-plate PDF
#' heatmaps summarizing genotype categories per target for each sample well.
#'
#' @param out_dir Output directory of a completed run. Must contain
#'   `editcall/editcall_summary.csv` with the same structure as `edit_result`.
#'   PDFs will be saved into `file.path(out_dir, "editviewer")`.
#' @param sample_list Path to a CSV file mapping index pairs to sample and plate
#'   layout. This file is expected to have NO header and exactly five columns in
#'   this order:
#'   1) `index_pair_id`, 2) `sample_name`, 3) `plate_id`, 4) `row_id`,
#'   5) `col_id`. Example rows:
#'
#'   miaoBC0001,Sample_A,1,A,1
#'   miaoBC0002,Sample_B,1,A,2
#'
#'   The values are used to annotate plots (sample name) and to facet by plate
#'   and well coordinates (row, col).
#' @param onefile If TRUE, draw plots in one PDF file, otherwise separate PDF files.
#' @param fill_plate If TRUE, fill missing wells on the plate layout.
#' @return Invisibly, a character vector of generated PDF file paths.
#' @export
editViewer <- function(out_dir, sample_list, onefile = FALSE, fill_plate = TRUE){
    if(!dir.exists(out_dir)){
        stop("out_dir does not exist: ", out_dir)
    }

    plot_dir <- file.path(out_dir, "editviewer")
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

    csv_fn <- file.path(out_dir, "editcall", "editcall_summary.csv")
    if(!file.exists(csv_fn)){
        stop("Cannot find editcall summary: ", csv_fn)
    }
    edit_result <- read.csv(csv_fn)

    edit_result <- subset(edit_result,
                          select = -data_type,
                          subset = data_type == "genotype")

    sample_list <- read.csv(sample_list, header = FALSE)
    hit <- match(edit_result$index_pair_id, sample_list$V1)
    edit_result$name <- sample_list$V2[hit]
    edit_result$plate <- sample_list$V3[hit]
    edit_result$row <- sample_list$V4[hit]
    edit_result$col <- sample_list$V5[hit]

    if("total_reads" %in% names(edit_result)){
        edit_result$total_reads_display <- paste0("n=", edit_result$total_reads)
    } else {
        edit_result$total_reads_display <- ""
    }
    edit_result_pat <- apply(subset(edit_result,
                                    select = -c(index_pair_id, name:col)),
                             1,
                             paste,
                             collapse = "_")
    dup_list <- tapply(seq_along(edit_result_pat), sub("_.+", "", edit_result$name), function(i){
        dup <- !duplicated(edit_result_pat[i])
        out <- data.frame(i = i, dup = dup)
        return(out)
    })
    dup_list <- do.call("rbind", dup_list)
    dup_list <- dup_list[order(dup_list$i), ]
    edit_result$uniq <- "Dup"
    edit_result$uniq[dup_list$dup] <- "Uniq"
    edit_result$uniq[edit_result$name == "" | is.na(edit_result$name)] <- ""
    long_edit_result <- tidyr::pivot_longer(edit_result,
                                            cols = -c(index_pair_id, sample_name:uniq),
                                            names_to = "gene",
                                            values_to = "edit")

    long_edit_result$edit_eval <- sapply(long_edit_result$edit, function(x){
        x <- unlist(strsplit(x, "/"))
        if(all(is.na(x))){
            return(NA)
        }
        x1 <- gsub("[0-9]", "", x)
        x2 <- as.numeric(gsub("[a-zA-Z]", "", x))
        is_ref <- x1 %in% "ref"
        is_inframe <- x2 %% 3 %in% 0
        if(all(is_ref)){
            return("ref")
        } else if(all(!is_ref)){
            if(all(!is_inframe)){
                return("alt")
            } else if(all(is_inframe)){
                return("alt_inframe_homo")
            } else {
                return("alt_inframe_het")
            }
        } else {
            if(all(!is_inframe)){
                return("het")
            } else {
                return("het_inframe")
            }
        }
    })

    eval_levels <- c("ref",
                     "alt", "alt_inframe_het", "alt_inframe_homo",
                     "het", "het_inframe")
    long_edit_result$edit_eval <- factor(long_edit_result$edit_eval, eval_levels)
    n_gene <- length(unique(long_edit_result$gene))

    if(fill_plate){
        wells <- paste(long_edit_result$plate,
                       long_edit_result$row,
                       long_edit_result$col, sep = "_")
        all_wells <- expand.grid(plate = unique(long_edit_result$plate),
                                 row = LETTERS[1:8],
                                 col = as.character(1:12))
        all_wells_id <- apply(all_wells, 1, paste, collapse = "_")
        missing_wells <- !all_wells_id %in% wells
        add_result <- long_edit_result[rep(1, sum(missing_wells)), ]
        add_result$name <- add_result$total_reads_display <- add_result$edit_eval <- add_result$uniq <- NA
        add_result$plate <- all_wells$plate[missing_wells]
        add_result$row <- all_wells$row[missing_wells]
        add_result$col <- all_wells$col[missing_wells]
        add_result <- lapply(unique(long_edit_result$gene), function(x){
            add_result$gene <- x
            return(add_result)
        })
        add_result <- do.call(rbind, add_result)
        long_edit_result <- rbind(long_edit_result, add_result)
    }

    long_edit_result$col <- factor(long_edit_result$col, 1:12)
    out_files <- character(0)
    if(onefile){
        pdf_fn <- file.path(plot_dir, "edit_viewer_plate.pdf")
        pdf(pdf_fn, width = 11.69, height = 8.27)
    }
    for(i in unique(long_edit_result$plate)){
        p <- ggplot2::ggplot(subset(long_edit_result, plate == i)) +
            ggplot2::geom_tile(ggplot2::aes(x = gene, y = 0, fill = edit_eval)) +
            ggplot2::geom_text(ggplot2::aes(x = (n_gene + 1) / 2, y = 2, label = name), vjust = 1, hjust = 0.5, size = 4) +
            ggplot2::geom_text(ggplot2::aes(x = (n_gene + 1) / 2, y = 1.4, label = uniq), vjust = 1, hjust = 0.5, size = 3) +
            ggplot2::geom_text(ggplot2::aes(x = (n_gene + 1) / 2, y = 0.9, label = total_reads_display), vjust = 1, hjust = 0.5, size = 3) +
            ggplot2::facet_grid(rows = ggplot2::vars(row), cols = ggplot2::vars(col), switch = "y", drop = FALSE) +
            ggplot2::scale_fill_manual(values = c("yellow", "darkblue", "blue",
                                                  "lightblue", "darkgreen", "green"),
                                       breaks = eval_levels,
                                       name = NULL) +
            ggplot2::labs(title = paste0("Plate ", i)) +
            ggplot2::theme(axis.title = ggplot2::element_blank(),
                           axis.text.y = ggplot2::element_blank(),
                           axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, size = 7),
                           axis.ticks.y = ggplot2::element_blank(),
                           panel.grid = ggplot2::element_blank())

        if(!onefile){
            pdf_fn <- file.path(plot_dir, paste0("edit_viewer_plate", i, ".pdf"))
            pdf(pdf_fn, width = 11.69, height = 8.27)
            print(p)
            dev.off()
        } else {
            print(p)
        }
        out_files <- c(out_files, pdf_fn)
    }

    if(onefile){
        dev.off()
    }
    invisible(out_files)
}
