#!/usr/bin/env Rscript
# 20k demux feature ablation (baseline + each idea solo)
suppressPackageStartupMessages({
  if (!requireNamespace("devtools", quietly = TRUE))
    stop("devtools required")
})
devtools::load_all(".", quiet = TRUE)

BASE <- "cursor_dev/tool_test/miao20250812"
fastq <- file.path(BASE, "basecall_filt_20k.fq")
index_list <- "/home/ftom/01_wd/invitro_domestication/genomeEdit/edit_check/nanopore/index_list.csv"
blast_csv <- "/home/ftom/01_wd/invitro_domestication/genomeEdit/edit_check/nanopore/20250918_pipeline_test/demultiplex/demultiplex_list.csv"
out_root <- file.path(BASE, "ablation_20k")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(fastq), file.exists(index_list), file.exists(blast_csv))

blast <- utils::read.csv(blast_csv, stringsAsFactors = FALSE)
blast <- blast[!is.na(blast$index_pair_id) & nzchar(as.character(blast$index_pair_id)), ]
blast <- blast[!duplicated(blast$sseqid), ]

# All read ids from FASTQ (via short header skim: names from doDemultiplex outputs)
n_core <- 16L
chunk_size <- 20000L

runs <- list(
  baseline = list(),
  pair_constrained = list(pair_constrained = TRUE),
  sandwich_rescue = list(sandwich_rescue = TRUE),
  pass_b_full_index = list(pass_b_full_index = TRUE),
  one_end_constrained = list(one_end_constrained = TRUE),
  confidence_gate = list(confidence_gate = TRUE)
)

results <- list()
baseline_assign <- NULL
all_ids <- NULL

for (nm in names(runs)) {
  demult_dir <- file.path(out_root, nm)
  unlink(demult_dir, recursive = TRUE)
  dir.create(demult_dir, recursive = TRUE, showWarnings = FALSE)
  args <- c(
    list(
      fastq = fastq,
      demult_dir = demult_dir,
      index_list = index_list,
      n_core = n_core,
      split_reads = FALSE,
      end_window = 120L,
      chunk_size = chunk_size
    ),
    runs[[nm]]
  )
  cat("\n=== run:", nm, "===\n")
  t0 <- Sys.time()
  res <- do.call(doDemultiplex, args)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  all_read_ids <- unique(c(res$assignments$read_id, res$unassigned$read_id))
  if (is.null(all_ids)) all_ids <<- all_read_ids

  met <- evalDemultiplexOracle(res$assignments, blast, all_read_ids)
  row <- data.frame(
    run = nm,
    elapsed_sec = round(elapsed, 2),
    n_edlib = res$n_edlib,
    n_edlib_per_read = round(res$n_edlib / met$n_total, 2),
    n_total = met$n_total,
    n_assigned = met$n_assigned,
    assign_rate = round(met$assign_rate, 4),
    n_blast = met$n_blast,
    BLAST_both = met$n_both,
    n_agree = met$n_agree,
    pair_agree = round(met$pair_agree, 4),
    false_vs_BLAST = met$n_both - met$n_agree,
    new_vs_baseline_agree = NA_real_,
    n_new_vs_baseline = NA_integer_,
    n_high = NA_integer_,
    assign_rate_high = NA_real_,
    BLAST_both_high = NA_integer_,
    pair_agree_high = NA_real_,
    stringsAsFactors = FALSE
  )

  if (nm == "baseline") {
    baseline_assign <<- res$assignments
  } else if (!is.null(baseline_assign)) {
    base_ids <- unique(as.character(baseline_assign$read_id))
    new_ids <- setdiff(unique(as.character(res$assignments$read_id)), base_ids)
    row$n_new_vs_baseline <- length(new_ids)
    if (length(new_ids)) {
      bp <- stats::setNames(as.character(blast$index_pair_id), as.character(blast$sseqid))
      ap <- stats::setNames(as.character(res$assignments$index_pair_id),
                            as.character(res$assignments$read_id))
      both_new <- intersect(new_ids, names(bp))
      if (length(both_new)) {
        row$new_vs_baseline_agree <- mean(ap[both_new] == bp[both_new], na.rm = TRUE)
      }
    }
  }

  if (nm == "confidence_gate" && nrow(res$assignments)) {
    high <- res$assignments[res$assignments$assign_tier == "high", , drop = FALSE]
    m_hi <- evalDemultiplexOracle(high, blast, all_read_ids)
    row$n_high <- m_hi$n_assigned
    row$assign_rate_high <- round(m_hi$assign_rate, 4)
    row$BLAST_both_high <- m_hi$n_both
    row$pair_agree_high <- round(m_hi$pair_agree, 4)
  }

  cat("elapsed:", row$elapsed_sec, "s  assign:", row$n_assigned,
      sprintf("(%.1f%%)", 100 * row$assign_rate),
      " pair_agree:", row$pair_agree, "\n")
  results[[nm]] <- row
  write.table(res$assignments, file.path(demult_dir, "assignments.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

tab <- do.call(rbind, results)
rownames(tab) <- NULL
tsv <- file.path(BASE, "ablation_metrics.tsv")
utils::write.table(tab, tsv, sep = "\t", quote = FALSE, row.names = FALSE)

md <- file.path(BASE, "ablation_results.md")
lines <- c(
  "# Demux ablation results (20k)",
  "",
  paste0("Date: ", format(Sys.time(), "%Y-%m-%d %H:%M")),
  paste0("FASTQ: `", fastq, "`"),
  paste0("n_core=", n_core, ", end_window=120"),
  "",
  "| run | assign_rate | BLAST_both | pair_agree | new_agree | n_new | elapsed | edlib/read |",
  "|-----|-------------|------------|------------|-----------|-------|---------|------------|"
)
for (i in seq_len(nrow(tab))) {
  r <- tab[i, ]
  na_g <- if (is.na(r$new_vs_baseline_agree)) "—" else sprintf("%.3f", r$new_vs_baseline_agree)
  na_n <- if (is.na(r$n_new_vs_baseline)) "—" else as.character(r$n_new_vs_baseline)
  lines <- c(lines, sprintf(
    "| %s | %.4f | %d | %.4f | %s | %s | %.2f | %.1f |",
    r$run, r$assign_rate, r$BLAST_both, r$pair_agree, na_g, na_n,
    r$elapsed_sec, r$n_edlib_per_read
  ))
}
if (!is.na(tab$pair_agree_high[tab$run == "confidence_gate"][1])) {
  rh <- tab[tab$run == "confidence_gate", ][1, ]
  lines <- c(
    lines, "",
    "### confidence_gate high-tier only",
    "",
    sprintf("- n_high: %d (assign_rate_high=%.4f)", rh$n_high, rh$assign_rate_high),
    sprintf("- BLAST_both_high: %d  pair_agree_high: %.4f",
            rh$BLAST_both_high, rh$pair_agree_high)
  )
}
lines <- c(lines, "", paste0("Raw metrics: `", tsv, "`"))
writeLines(lines, md)
cat("\nWrote", tsv, "and", md, "\n")
print(tab[, c("run", "assign_rate", "BLAST_both", "pair_agree",
              "new_vs_baseline_agree", "elapsed_sec")])
