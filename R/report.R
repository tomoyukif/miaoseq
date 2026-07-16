################################################################################
# Reporting and visualization
################################################################################

#' @title Evaluate MIAO results
#' @description Summarize demultiplex, amplicon resolve, and edit-calling outputs
#'   under a single run directory.
#' @param out_dir Output directory containing `demultiplex/`, and optionally
#'   `amplicon_assign/`, `amplicon/`, and `editcall/` subdirectories.
#' @param output_reads If `TRUE`, write read IDs that received no gene assignment
#'   to `miao_summary/amplicon_unassigned_read_ids.tsv`.
#' @return Invisibly, the path to `miao_summary/`.
#' @importFrom utils read.delim write.table
#' @export
evalMiao <- function(out_dir, output_reads = FALSE) {
    demult_dir <- file.path(out_dir, "demultiplex")
    amplicon_assign_dir <- file.path(out_dir, "amplicon_assign")
    amplicon_dir <- file.path(out_dir, "amplicon")
    editcall_dir <- file.path(out_dir, "editcall")
    summary_dir <- file.path(out_dir, "miao_summary")
    dir.create(summary_dir, showWarnings = FALSE, recursive = TRUE)

    assignments_fn <- file.path(demult_dir, "assignments.tsv")
    if (!file.exists(assignments_fn)) {
        stop("Cannot find demultiplex assignments: ", assignments_fn)
    }
    assignments <- utils::read.delim(assignments_fn, stringsAsFactors = FALSE)

    unassigned_fn <- file.path(demult_dir, "unassigned.tsv")
    unassigned <- if (file.exists(unassigned_fn)) {
        utils::read.delim(unassigned_fn, stringsAsFactors = FALSE)
    } else {
        data.frame(read_id = character(), reason = character(), stringsAsFactors = FALSE)
    }

    n_demult_reads <- nrow(assignments)
    n_undemult <- nrow(unassigned)
    n_input_reads <- n_demult_reads + n_undemult

    n_unassign_by_reason <- table(unassigned$reason)
    n_indexed_reads_per_index_f <- table(assignments$f_index_id)
    n_indexed_reads_per_index_r <- table(assignments$r_index_id)
    prop_indexed_reads_per_index_f <- n_indexed_reads_per_index_f / sum(n_indexed_reads_per_index_f)
    prop_indexed_reads_per_index_r <- n_indexed_reads_per_index_r / sum(n_indexed_reads_per_index_r)

    gene_assignments <- NULL
    gene_assignments_fn <- file.path(amplicon_assign_dir, "gene_assignments.tsv")
    if (!file.exists(gene_assignments_fn)) {
        gene_assignments_fn <- file.path(amplicon_dir, "gene_assignments.tsv")
    }
    if (file.exists(gene_assignments_fn)) {
        gene_assignments <- utils::read.delim(gene_assignments_fn, stringsAsFactors = FALSE)
    }

    n_gene_assigned_reads <- NA_integer_
    n_gene_unassigned_reads <- NA_integer_
    n_gene_reads_per_gene <- NULL
    prop_gene_reads_per_gene <- NULL
    if (!is.null(gene_assignments) && nrow(gene_assignments) > 0L) {
        is_assigned <- !is.na(gene_assignments$gene_id) &
            gene_assignments$gene_id != "" &
            gene_assignments$gene_id != "NA"
        n_gene_assigned_reads <- sum(is_assigned)
        n_gene_unassigned_reads <- sum(!is_assigned)
        n_gene_reads_per_gene <- table(gene_assignments$gene_id[is_assigned])
        prop_gene_reads_per_gene <- signif(
            n_gene_reads_per_gene / sum(n_gene_reads_per_gene),
            3
        )
    }

    n_edicall_reads <- NA_integer_
    n_edicall_reads_per_gene <- NULL
    prop_edicall_reads_per_gene <- NULL
    editcall_fn <- file.path(editcall_dir, "editcall_filtered.csv")
    if (file.exists(editcall_fn)) {
        editcall_out <- utils::read.csv(editcall_fn, stringsAsFactors = FALSE)
        n_edicall_reads <- sum(editcall_out$count)
        n_edicall_reads_per_gene <- tapply(
            editcall_out$count,
            editcall_out$target_gene,
            sum
        )
        if (!is.null(n_gene_reads_per_gene)) {
            common_genes <- intersect(
                names(n_gene_reads_per_gene),
                names(n_edicall_reads_per_gene)
            )
            prop_edicall_reads_per_gene <- n_edicall_reads_per_gene[common_genes] /
                n_gene_reads_per_gene[common_genes]
        }
    }

    if (isTRUE(output_reads) && !is.null(gene_assignments)) {
        unassigned_gene_reads <- gene_assignments$read_id[
            is.na(gene_assignments$gene_id) |
                gene_assignments$gene_id == "" |
                gene_assignments$gene_id == "NA"
        ]
        utils::write.table(
            data.frame(read_id = unassigned_gene_reads),
            file.path(summary_dir, "amplicon_unassigned_read_ids.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE
        )
    }

  pct <- function(x, denom) {
    if (is.na(x) || is.na(denom) || denom == 0L) {
      return(NA)
    }
    signif(x / denom, 3) * 100
  }

    out1 <- rbind(
        c("Input reads: ", n_input_reads, ""),
        c("Demultiplexed reads: ", n_demult_reads, pct(n_demult_reads, n_input_reads)),
        c("Gene-assigned reads: ", n_gene_assigned_reads, pct(n_gene_assigned_reads, n_input_reads)),
        c("Amplicon unassigned reads: ", n_gene_unassigned_reads, pct(n_gene_unassigned_reads, n_input_reads)),
        c("Editcalled reads: ", n_edicall_reads, pct(n_edicall_reads, n_input_reads)),
        c("Undemultiplexed reads: ", n_undemult, pct(n_undemult, n_input_reads))
    )
    if (length(n_unassign_by_reason) > 0L) {
        reason_rows <- cbind(
            paste0("Undemultiplexed (", names(n_unassign_by_reason), "): "),
            as.character(n_unassign_by_reason),
            signif(as.numeric(n_unassign_by_reason) / n_input_reads, 3) * 100
        )
        out1 <- rbind(out1, reason_rows)
    }
    utils::write.table(
        out1,
        file.path(summary_dir, "read_stats.tsv"),
        row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE
    )

    out2 <- rbind(
        names(n_indexed_reads_per_index_f),
        n_indexed_reads_per_index_f,
        signif(prop_indexed_reads_per_index_f, 3) * 100,
        names(n_indexed_reads_per_index_r),
        n_indexed_reads_per_index_r,
        signif(prop_indexed_reads_per_index_r, 3) * 100
    )
    out2 <- cbind(
        c(
            "",
            "Indexed reads per forward index",
            "Proportion of indexed reads per forward index",
            "",
            "Indexed reads per reverse index",
            "Proportion of indexed reads per reverse index"
        ),
        out2
    )
    utils::write.table(
        out2,
        file.path(summary_dir, "indexed_reads_per_index.tsv"),
        row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE
    )

    if (!is.null(n_gene_reads_per_gene)) {
        gene_names <- names(n_gene_reads_per_gene)
        out3 <- rbind(
            gene_names,
            as.character(n_gene_reads_per_gene),
            as.character(prop_gene_reads_per_gene)
        )
        out3_labels <- c(
            "",
            "Gene-assigned reads per gene",
            "Proportion of gene-assigned reads per gene"
        )
        if (!is.null(n_edicall_reads_per_gene)) {
            editcall_vals <- ifelse(
                gene_names %in% names(n_edicall_reads_per_gene),
                as.character(n_edicall_reads_per_gene[gene_names]),
                NA_character_
            )
            prop_vals <- ifelse(
                gene_names %in% names(prop_edicall_reads_per_gene),
                as.character(signif(prop_edicall_reads_per_gene[gene_names], 3)),
                NA_character_
            )
            out3 <- rbind(out3, editcall_vals, prop_vals)
            out3_labels <- c(
                out3_labels,
                "Editcalled reads per gene",
                "Proportion of editcalled reads per gene-assigned reads"
            )
        }
        out3 <- cbind(out3_labels, out3)
        utils::write.table(
            out3,
            file.path(summary_dir, "gene_reads_per_gene.tsv"),
            row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE
        )
    }

    invisible(summary_dir)
}

#' Visualize edit-calling results per plate as PDF heatmaps
#'
#' Reads `editcall/editcall_summary.csv` and draws per-plate genotype heatmaps.
#' Genotype labels (`ref` / `sub` / `delN` / `insN` / `indelD-I`) are parsed
#' structurally; SNP (`sub`) is never treated as frameshift. In-frame classes
#' use **net indel length mod 3** (not CDS phase).
#'
#' @param out_dir Output directory of a completed run. Must contain
#'   `editcall/editcall_summary.csv`. PDFs are written to `editviewer/`.
#' @param sample_list Path to a headerless CSV with five columns:
#'   `index_pair_id`, sample name, plate id, row id (`A`–`H`), column id (`1`–`12`).
#'   This must be a plate layout file, not `index_list.csv`. When `NULL`, layout
#'   is taken from `editcall_summary.csv` only if `row_id` / `col_id` already
#'   look like well coordinates.
#' @param onefile If `TRUE`, draw plots in one PDF file, otherwise separate PDF files.
#' @param fill_plate If `TRUE`, fill missing wells on the plate layout.
#' @return Invisibly, a character vector of generated PDF file paths.
#' @export
editViewer <- function(out_dir, sample_list = NULL, onefile = FALSE, fill_plate = TRUE) {
    if (!dir.exists(out_dir)) {
        stop("out_dir does not exist: ", out_dir)
    }

    plot_dir <- file.path(out_dir, "editviewer")
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

    csv_fn <- file.path(out_dir, "editcall", "editcall_summary.csv")
    if (!file.exists(csv_fn)) {
        stop("Cannot find editcall summary: ", csv_fn)
    }
    edit_result <- utils::read.csv(csv_fn, stringsAsFactors = FALSE)
    edit_result <- edit_result[edit_result$data_type == "genotype", , drop = FALSE]
    edit_result$data_type <- NULL

    meta_cols <- c(
        "index_pair_id", "sample_name", "plate_id", "row_id", "col_id", "total_reads"
    )
    gene_cols <- setdiff(names(edit_result), meta_cols)

    if (!is.null(sample_list)) {
        layout <- .read_plate_layout(sample_list)
        edit_result <- .join_plate_layout(edit_result, layout)
    } else if (.has_valid_plate_layout(edit_result)) {
        edit_result$name <- if ("sample_name" %in% names(edit_result)) {
            edit_result$sample_name
        } else {
            edit_result$index_pair_id
        }
        edit_result$plate <- edit_result$plate_id
        edit_result$row <- toupper(trimws(as.character(edit_result$row_id)))
        edit_result$col <- as.character(suppressWarnings(as.integer(edit_result$col_id)))
    } else {
        stop(
            "Plate layout is missing or invalid in editcall_summary.csv. ",
            "Pass sample_list with columns: index_pair_id, sample_name, plate, ",
            "row (A-H), col (1-12). Do not use index_list.csv."
        )
    }

    if (nrow(edit_result) < 1L) {
        stop("No samples with valid plate layout to plot.")
    }

    if ("total_reads" %in% names(edit_result)) {
        edit_result$total_reads_display <- paste0("n=", edit_result$total_reads)
    } else {
        edit_result$total_reads_display <- ""
    }

    edit_result$uniq <- ""
    edit_result_pat <- apply(
        edit_result[, gene_cols, drop = FALSE],
        1,
        paste,
        collapse = "_"
    )
    dup_list <- tapply(seq_along(edit_result_pat), edit_result$name, function(i) {
        data.frame(i = i, dup = !duplicated(edit_result_pat[i]))
    })
    dup_list <- do.call("rbind", dup_list)
    dup_list <- dup_list[order(dup_list$i), , drop = FALSE]
    edit_result$uniq <- "Dup"
    edit_result$uniq[dup_list$dup] <- "Uniq"
    edit_result$uniq[edit_result$name == "" | is.na(edit_result$name)] <- ""

    long_edit_result <- tidyr::pivot_longer(
        edit_result,
        cols = dplyr::all_of(gene_cols),
        names_to = "gene",
        values_to = "edit"
    )

    long_edit_result$edit_eval <- vapply(
        long_edit_result$edit,
        .classify_genotype_edit_eval,
        character(1)
    )

    eval_levels <- c(
        "ref",
        "alt", "alt_inframe_het", "alt_inframe_homo",
        "het", "het_inframe"
    )
    long_edit_result$edit_eval <- factor(long_edit_result$edit_eval, eval_levels)
    n_gene <- length(unique(long_edit_result$gene))

    long_edit_result$col <- as.character(long_edit_result$col)
    long_edit_result$row <- as.character(long_edit_result$row)
    long_edit_result$plate <- as.character(long_edit_result$plate)

    if (fill_plate) {
        wells <- paste(
            long_edit_result$plate,
            long_edit_result$row,
            long_edit_result$col,
            sep = "_"
        )
        all_wells <- expand.grid(
            plate = unique(long_edit_result$plate),
            row = LETTERS[1:8],
            col = as.character(1:12),
            stringsAsFactors = FALSE
        )
        all_wells_id <- apply(all_wells, 1, paste, collapse = "_")
        missing_wells <- !all_wells_id %in% wells
        if (any(missing_wells)) {
            add_result <- long_edit_result[rep(1L, sum(missing_wells)), , drop = FALSE]
            add_result$name <- add_result$total_reads_display <- NA_character_
            add_result$edit_eval <- add_result$uniq <- NA_character_
            add_result$plate <- all_wells$plate[missing_wells]
            add_result$row <- all_wells$row[missing_wells]
            add_result$col <- all_wells$col[missing_wells]
            add_result <- do.call(
                "rbind",
                lapply(unique(long_edit_result$gene), function(g) {
                    add_result$gene <- g
                    add_result
                })
            )
            long_edit_result <- rbind(long_edit_result, add_result)
        }
    }

    long_edit_result$col <- factor(long_edit_result$col, levels = as.character(1:12))

    out_files <- character(0)
    if (onefile) {
        pdf_fn <- file.path(plot_dir, "edit_viewer_plate.pdf")
        grDevices::pdf(pdf_fn, width = 11.69, height = 8.27)
    }
    for (i in unique(long_edit_result$plate)) {
        p <- ggplot2::ggplot(subset(long_edit_result, plate == i)) +
            ggplot2::geom_tile(ggplot2::aes(x = gene, y = 0, fill = edit_eval)) +
            ggplot2::geom_text(
                ggplot2::aes(x = (n_gene + 1) / 2, y = 2, label = name),
                vjust = 1, hjust = 0.5, size = 4
            ) +
            ggplot2::geom_text(
                ggplot2::aes(x = (n_gene + 1) / 2, y = 1.4, label = uniq),
                vjust = 1, hjust = 0.5, size = 3
            ) +
            ggplot2::geom_text(
                ggplot2::aes(x = (n_gene + 1) / 2, y = 0.9, label = total_reads_display),
                vjust = 1, hjust = 0.5, size = 3
            ) +
            ggplot2::facet_grid(
                rows = ggplot2::vars(row),
                cols = ggplot2::vars(col),
                switch = "y",
                drop = FALSE
            ) +
            ggplot2::scale_fill_manual(
                values = c("yellow", "darkblue", "blue", "lightblue", "darkgreen", "green"),
                breaks = eval_levels,
                name = NULL
            ) +
            ggplot2::labs(title = paste0("Plate ", i)) +
            ggplot2::theme(
                axis.title = ggplot2::element_blank(),
                axis.text.y = ggplot2::element_blank(),
                axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, size = 7),
                axis.ticks.y = ggplot2::element_blank(),
                panel.grid = ggplot2::element_blank()
            )

        if (!onefile) {
            pdf_fn <- file.path(plot_dir, paste0("edit_viewer_plate", i, ".pdf"))
            grDevices::pdf(pdf_fn, width = 11.69, height = 8.27)
            print(p)
            grDevices::dev.off()
        } else {
            print(p)
        }
        out_files <- c(out_files, pdf_fn)
    }

    if (onefile) {
        grDevices::dev.off()
    }
    invisible(out_files)
}

################################################################################
# Genotype label parsing for editViewer (net indel length class, not CDS frame)
################################################################################

#' Parse one allele label from editcall summary genotype strings.
#'
#' Recognized forms: `ref`, `sub`, `delN`, `insN`, `indelD-I`.
#' @keywords internal
.parse_allele_label <- function(label) {
    label <- as.character(label)[[1]]
    if (is.na(label) || !nzchar(label)) {
        return(list(
            type = NA_character_,
            n_del = NA_integer_,
            n_ins = NA_integer_,
            is_ref = NA,
            is_frameshift = NA
        ))
    }
    if (identical(label, "ref")) {
        return(list(
            type = "ref", n_del = 0L, n_ins = 0L,
            is_ref = TRUE, is_frameshift = FALSE
        ))
    }
    if (identical(label, "sub")) {
        return(list(
            type = "sub", n_del = 0L, n_ins = 0L,
            is_ref = FALSE, is_frameshift = FALSE
        ))
    }
    m_del <- regexec("^del([0-9]+)$", label)
    if (m_del[[1]][1] >= 0L) {
        n <- as.integer(regmatches(label, m_del)[[1]][2])
        return(list(
            type = "del", n_del = n, n_ins = 0L,
            is_ref = FALSE, is_frameshift = (n %% 3L) != 0L
        ))
    }
    m_ins <- regexec("^ins([0-9]+)$", label)
    if (m_ins[[1]][1] >= 0L) {
        n <- as.integer(regmatches(label, m_ins)[[1]][2])
        return(list(
            type = "ins", n_del = 0L, n_ins = n,
            is_ref = FALSE, is_frameshift = (n %% 3L) != 0L
        ))
    }
    m_indel <- regexec("^indel([0-9]+)-([0-9]+)$", label)
    if (m_indel[[1]][1] >= 0L) {
        parts <- regmatches(label, m_indel)[[1]]
        n_del <- as.integer(parts[2])
        n_ins <- as.integer(parts[3])
        net <- abs(n_ins - n_del)
        return(list(
            type = "indel", n_del = n_del, n_ins = n_ins,
            is_ref = FALSE, is_frameshift = (net %% 3L) != 0L
        ))
    }
    # Unrecognized label: treat as non-ref frameshift so it cannot look "clean".
    list(
        type = "unknown", n_del = NA_integer_, n_ins = NA_integer_,
        is_ref = FALSE, is_frameshift = TRUE
    )
}

#' Classify a genotype string (alleles joined by `/`) for plate coloring.
#'
#' Uses net indel length mod 3 for frameshift classes. `sub` is length-neutral
#' (never frameshift). This is not a CDS-phase codon model.
#' @keywords internal
.classify_genotype_edit_eval <- function(genotype) {
    if (length(genotype) != 1L || is.na(genotype) || !nzchar(as.character(genotype))) {
        return(NA_character_)
    }
    alleles <- unlist(strsplit(as.character(genotype), "/", fixed = TRUE), use.names = FALSE)
    alleles <- alleles[nzchar(alleles) & !is.na(alleles)]
    if (length(alleles) < 1L) {
        return(NA_character_)
    }
    parsed <- lapply(alleles, .parse_allele_label)
    is_ref <- vapply(parsed, function(p) isTRUE(p$is_ref), logical(1))
    is_fs <- vapply(parsed, function(p) isTRUE(p$is_frameshift), logical(1))

    if (all(is_ref)) {
        return("ref")
    }

    alt_fs <- is_fs[!is_ref]
    if (length(alt_fs) < 1L) {
        return(NA_character_)
    }
    alt_inframe <- !alt_fs

    if (all(!is_ref)) {
        if (all(alt_fs)) {
            "alt"
        } else if (all(alt_inframe)) {
            "alt_inframe_homo"
        } else {
            "alt_inframe_het"
        }
    } else if (all(alt_fs)) {
        "het"
    } else {
        "het_inframe"
    }
}
