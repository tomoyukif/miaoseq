suppressPackageStartupMessages(devtools::load_all("."))
BASE <- "cursor_dev/tool_test/miao20250812"
fastq <- "/home/ftom/01_wd/invitro_domestication/genomeEdit/edit_check/nanopore/20250918_pipeline_test/basecall/basecall_filt.fq"
index_list <- "/home/ftom/01_wd/invitro_domestication/genomeEdit/edit_check/nanopore/index_list.csv"
demult_dir <- file.path(BASE, "demultiplex_full")
# Fresh outputs (keep prior diagnostics overwritable)
unlink(file.path(demult_dir, c(
  "assignments.tsv", "unassigned.tsv", "summary_by_sample.tsv", "run_stats.txt"
)))
dir.create(demult_dir, recursive = TRUE, showWarnings = FALSE)

t0 <- Sys.time()
res <- doDemultiplex(
  fastq = fastq,
  demult_dir = demult_dir,
  index_list = index_list,
  n_core = 32,
  split_reads = FALSE,
  end_window = 120,
  # Large chunks: dict is rebuilt per demux_reads_cpp() call
  chunk_size = 100000
)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
n_a <- nrow(res$assignments)
n_u <- nrow(res$unassigned)
cat("elapsed_sec:", round(elapsed, 1), "\n")
cat("n_assigned:", n_a, "\n")
cat("n_unassigned:", n_u, "\n")
cat("assign_rate:", round(n_a / (n_a + n_u), 4), "\n")
cat("n_edlib:", format(res$n_edlib, scientific = FALSE), "\n")
cat("n_edlib_per_read:", round(res$n_edlib / (n_a + n_u), 2), "\n")
cat("n_samples:", nrow(res$summary), "\n")
print(table(res$assignments$match_class))
print(table(res$assignments$anchor_status))
print(table(res$unassigned$reason))
print(head(res$summary[order(-res$summary$n_reads),
                       c("sample_id", "n_reads", "n_complete", "n_fuzzy")], 15))

old <- "/home/ftom/01_wd/invitro_domestication/genomeEdit/edit_check/nanopore/20250918_pipeline_test/demultiplex/demultiplex_list.csv"
agree <- NA_real_
both_n <- 0L
if (file.exists(old)) {
  od <- read.csv(old, stringsAsFactors = FALSE)
  od <- od[!is.na(od$index_pair_id) & nzchar(od$index_pair_id), ]
  od <- od[!duplicated(od$sseqid), ]
  m <- merge(res$assignments[, c("read_id", "index_pair_id")],
             od[, c("sseqid", "index_pair_id")],
             by.x = "read_id", by.y = "sseqid", all = TRUE,
             suffixes = c(".edlib", ".blast"))
  both <- !is.na(m$index_pair_id.edlib) & !is.na(m$index_pair_id.blast)
  both_n <- sum(both)
  agree <- if (both_n) mean(m$index_pair_id.edlib[both] == m$index_pair_id.blast[both]) else NA_real_
  cat("BLAST assigned (unique, with pair):", nrow(od), "\n")
  cat("both:", both_n, "agree:", agree, "\n")
  cat("edlib only:", sum(!is.na(m$index_pair_id.edlib) & is.na(m$index_pair_id.blast)), "\n")
  cat("blast only:", sum(is.na(m$index_pair_id.edlib) & !is.na(m$index_pair_id.blast)), "\n")
}

stats <- c(
  paste("elapsed_sec", round(elapsed, 1)),
  paste("n_assigned", n_a),
  paste("n_unassigned", n_u),
  paste("assign_rate", round(n_a / (n_a + n_u), 4)),
  paste("n_edlib", format(res$n_edlib, scientific = FALSE)),
  paste("blast_both", both_n),
  paste("blast_agree", round(agree, 4))
)
writeLines(stats, file.path(demult_dir, "run_stats.txt"))
cat("Wrote", file.path(demult_dir, "run_stats.txt"), "\n")
