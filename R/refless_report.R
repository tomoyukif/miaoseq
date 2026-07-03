################################################################################
# Quality report generation for refless pipeline
################################################################################

#' Identify likely contaminant clusters (small fraction, low support)
#' @keywords internal
.flag_contaminant_clusters <- function(consensus_df, confidence_df, min_fraction = 0.05){
    if(is.null(consensus_df) || nrow(consensus_df) == 0){
        return(data.frame())
    }
    out <- consensus_df
    if(!is.null(confidence_df) && nrow(confidence_df)){
        out <- merge(out, confidence_df[, c("cluster_id", "mean_identity", "min_base_support")],
                     by = "cluster_id", all.x = TRUE)
    }
    out$is_contaminant <- out$fraction < min_fraction
    out
}

#' Build per-sample quality metrics from pipeline results
#' @keywords internal
.collect_sample_qc <- function(sample_id, sample_res, demult_out){
    demux <- demult_out[demult_out$sample_id == sample_id, , drop = FALSE]
    n_assigned <- nrow(demux)
    class_tab <- if(n_assigned) table(demux$class) else table(character())

    cluster_df <- sample_res$cluster
    consensus_df <- sample_res$consensus
    confidence_df <- sample_res$confidence
    chimera_df <- sample_res$chimera

    n_reads <- if(!is.null(cluster_df)) nrow(cluster_df) else 0
    n_clusters <- if(!is.null(consensus_df)) nrow(consensus_df) else 0
    primary_fraction <- if(!is.null(consensus_df) && nrow(consensus_df)){
        max(consensus_df$fraction, na.rm = TRUE)
    } else {
        NA_real_
    }
    n_contaminant <- if(!is.null(consensus_df)){
        sum(consensus_df$fraction < 0.05, na.rm = TRUE)
    } else {
        0L
    }
    n_chimera <- if(!is.null(chimera_df)){
        sum(chimera_df$is_chimera, na.rm = TRUE)
    } else {
        NA_integer_
    }
    chimera_fraction <- if(!is.null(chimera_df) && nrow(chimera_df)){
        mean(chimera_df$is_chimera)
    } else {
        NA_real_
    }

    data.frame(
        sample_id = sample_id,
        assigned_reads = n_assigned,
        high_confidence = as.integer(if("High" %in% names(class_tab)) class_tab[["High"]] else 0),
        low_confidence = as.integer(if("Low" %in% names(class_tab)) class_tab[["Low"]] else 0),
        ambiguous = as.integer(if("Ambiguous" %in% names(class_tab)) class_tab[["Ambiguous"]] else 0),
        clustered_reads = n_reads,
        n_clusters = n_clusters,
        primary_cluster_fraction = primary_fraction,
        contaminant_clusters = n_contaminant,
        chimera_reads = n_chimera,
        chimera_fraction = chimera_fraction,
        stringsAsFactors = FALSE
    )
}

#' Write simple HTML quality report
#' @keywords internal
.write_refless_html_report <- function(out_dir, run_summary, sample_qc, demult_summary){
    html_fn <- file.path(out_dir, "quality_report.html")
    rows <- apply(sample_qc, 1, function(r){
        sprintf(
            "<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>",
            r["sample_id"], r["assigned_reads"], r["n_clusters"],
            ifelse(is.na(r["primary_cluster_fraction"]), "-",
                   sprintf("%.1f%%", 100 * as.numeric(r["primary_cluster_fraction"]))),
            r["contaminant_clusters"],
            ifelse(is.na(r["chimera_reads"]), "-", r["chimera_reads"]),
            ifelse(is.na(r["chimera_fraction"]), "-",
                   sprintf("%.1f%%", 100 * as.numeric(r["chimera_fraction"])))
        )
    })
    html <- paste0(
        "<!DOCTYPE html><html><head><meta charset='utf-8'>",
        "<title>miaoseq refless quality report</title>",
        "<style>body{font-family:sans-serif;margin:2em}",
        "table{border-collapse:collapse}td,th{border:1px solid #ccc;padding:4px 8px}",
        "th{background:#f0f0f0}</style></head><body>",
        "<h1>miaoseq refless quality report</h1>",
        "<h2>Run summary</h2><ul>",
        sprintf("<li>Total reads processed: %s</li>", run_summary$total_reads),
        sprintf("<li>Assigned (High+Low): %s (%.1f%%)</li>",
                run_summary$assigned_reads,
                100 * run_summary$assigned_fraction),
        sprintf("<li>Unclassified: %s (%.1f%%)</li>",
                run_summary$unclassified_reads,
                100 * run_summary$unclassified_fraction),
        sprintf("<li>Ambiguous: %s (%.1f%%)</li>",
                run_summary$ambiguous_reads,
                100 * run_summary$ambiguous_fraction),
        sprintf("<li>Samples processed: %s</li>", run_summary$n_samples),
        "</ul>",
        "<h2>Demultiplex classification</h2><ul>",
        paste(sprintf("<li>%s: %s</li>", names(demult_summary), demult_summary), collapse = ""),
        "</ul>",
        "<h2>Per-sample metrics</h2>",
        "<table><tr><th>Sample</th><th>Assigned reads</th><th>Clusters</th>",
        "<th>Primary fraction</th><th>Contaminant clusters</th>",
        "<th>Chimera reads</th><th>Chimera fraction</th></tr>",
        paste(rows, collapse = ""),
        "</table></body></html>"
    )
    writeLines(html, html_fn)
    html_fn
}

#' Generate quality report for refless pipeline output
#'
#' Writes TSV summaries and an HTML report with demultiplex rates, cluster counts,
#' contaminant candidates, and chimera statistics.
#'
#' @param pipeline_results Output from \code{miaoPipeline}.
#' @param out_dir Report output directory (default: \code{file.path(pipeline_out, "report")}).
#' @param pipeline_out Root pipeline output directory (used when \code{out_dir} is NULL).
#' @param edit_results Optional output from \code{doEditCalling}.
#' @param contaminant_fraction Threshold below which clusters are flagged as contaminants.
#' @return List with paths to generated report files.
#' @export
writeReflessReport <- function(pipeline_results,
                               out_dir = NULL,
                               pipeline_out = NULL,
                               edit_results = NULL,
                               contaminant_fraction = 0.05){
    if(is.null(out_dir)){
        if(is.null(pipeline_out)){
            stop("Provide out_dir or pipeline_out")
        }
        out_dir <- file.path(pipeline_out, "report")
    }
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    demult_out <- pipeline_results$demultiplex
    total_reads <- nrow(demult_out)
    class_tab <- table(demult_out$class)
    assigned <- sum(demult_out$class %in% c("High", "Low"), na.rm = TRUE)
    unclassified <- as.integer(if("Unclassified" %in% names(class_tab)) class_tab[["Unclassified"]] else 0)
    ambiguous <- as.integer(if("Ambiguous" %in% names(class_tab)) class_tab[["Ambiguous"]] else 0)

    run_summary <- data.frame(
        total_reads = total_reads,
        assigned_reads = assigned,
        assigned_fraction = if(total_reads) assigned / total_reads else 0,
        unclassified_reads = unclassified,
        unclassified_fraction = if(total_reads) unclassified / total_reads else 0,
        ambiguous_reads = ambiguous,
        ambiguous_fraction = if(total_reads) ambiguous / total_reads else 0,
        n_samples = length(pipeline_results$samples),
        stringsAsFactors = FALSE
    )
    write.csv(run_summary, file.path(out_dir, "run_summary.csv"), row.names = FALSE)
    write.csv(as.data.frame(class_tab), file.path(out_dir, "demultiplex_class_counts.csv"),
              row.names = FALSE)

    sample_qc <- do.call(rbind, lapply(names(pipeline_results$samples), function(sid){
        .collect_sample_qc(sid, pipeline_results$samples[[sid]], demult_out)
    }))
    if(is.null(sample_qc)){
        sample_qc <- data.frame()
    }
    write.csv(sample_qc, file.path(out_dir, "sample_quality.csv"), row.names = FALSE)

    contaminant_rows <- list()
    for(sid in names(pipeline_results$samples)){
        s <- pipeline_results$samples[[sid]]
        if(is.null(s$consensus)){
            next
        }
        flagged <- .flag_contaminant_clusters(s$consensus, s$confidence, contaminant_fraction)
        if(nrow(flagged)){
            flagged$sample_id <- sid
            contaminant_rows[[sid]] <- flagged
        }
    }
    contaminant_df <- if(length(contaminant_rows)) do.call(rbind, contaminant_rows) else data.frame()
    write.csv(contaminant_df, file.path(out_dir, "contaminant_clusters.csv"), row.names = FALSE)

    edit_summary <- NULL
    if(!is.null(edit_results) && !is.null(edit_results$primary)){
        edit_summary <- edit_results$primary
        write.csv(edit_summary, file.path(out_dir, "edit_primary_summary.csv"), row.names = FALSE)
    }

    html_fn <- .write_refless_html_report(out_dir, run_summary, sample_qc, class_tab)

    invisible(list(
        run_summary = file.path(out_dir, "run_summary.csv"),
        sample_quality = file.path(out_dir, "sample_quality.csv"),
        contaminant_clusters = file.path(out_dir, "contaminant_clusters.csv"),
        html_report = html_fn,
        edit_primary = if(!is.null(edit_summary)) file.path(out_dir, "edit_primary_summary.csv") else NULL
    ))
}
