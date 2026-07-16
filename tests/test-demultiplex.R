#!/usr/bin/env Rscript
# Minimal demultiplex tests for allow_single_end mode.

suppressPackageStartupMessages({
  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("devtools required")
  }
  devtools::load_all(".", quiet = TRUE)
})

tmpdir <- tempfile("demux_test_")
dir.create(tmpdir, recursive = TRUE)

index_path <- file.path(tmpdir, "index_list.csv")
writeLines(
  readLines(system.file("extdata", "index_list.csv", package = "miaoseq"))[1:2],
  index_path
)

idx <- utils::read.csv(index_path, header = FALSE, stringsAsFactors = FALSE)
f1 <- idx[1, 3]
r1 <- idx[1, 5]
f2 <- idx[2, 3]
r2 <- idx[2, 5]

pad <- paste(rep("A", 400), collapse = "")
single_f <- paste0(f1, pad)
dual_ok <- paste0(f1, pad, r1)
conflict <- paste0(f1, pad, r2)

fq <- file.path(tmpdir, "reads.fq")
writeLines(c(
  "@single_f", single_f, "+", paste(rep("I", nchar(single_f)), collapse = ""),
  "@dual_ok", dual_ok, "+", paste(rep("I", nchar(dual_ok)), collapse = ""),
  "@conflict", conflict, "+", paste(rep("I", nchar(conflict)), collapse = "")
), fq)

run_demux <- function(allow_single_end) {
  out <- file.path(tmpdir, if (allow_single_end) "single" else "dual")
  unlink(out, recursive = TRUE)
  doDemultiplex(
    fastq = fq,
    demult_dir = out,
    index_list = index_path,
    n_core = 1L,
    split_reads = FALSE,
    allow_single_end = allow_single_end
  )
}

res_dual <- run_demux(FALSE)
res_single <- run_demux(TRUE)

stopifnot(nrow(res_dual$assignments) == 1L)
stopifnot(res_dual$assignments$read_id == "dual_ok")
stopifnot(res_dual$assignments$assign_mode == "dual_end")
stopifnot("single_f" %in% res_dual$unassigned$read_id)
stopifnot(res_dual$unassigned$reason[res_dual$unassigned$read_id == "single_f"] == "barcode_fail")

stopifnot(nrow(res_single$assignments) == 2L)
stopifnot(all(c("dual_ok", "single_f") %in% res_single$assignments$read_id))
single_row <- res_single$assignments[res_single$assignments$read_id == "single_f", ]
stopifnot(single_row$assign_mode %in% c("single_f", "single_r"))
stopifnot(single_row$index_pair_id == "miaoBC0001")
stopifnot("conflict" %in% res_single$unassigned$read_id)
stopifnot(res_single$unassigned$reason[res_single$unassigned$read_id == "conflict"] == "ambiguous_ends")

res_stats <- doDemultiplex(
  fastq = fq,
  demult_dir = file.path(tmpdir, "stats"),
  index_list = index_path,
  n_core = 1L,
  split_reads = FALSE,
  allow_single_end = TRUE,
  stats_unassign = TRUE
)
stopifnot(!is.null(res_stats$stats_unassigned))
stats <- res_stats$stats_unassigned
stopifnot(nrow(stats) == 2L)
stopifnot(all(c("scope", "sample_id", "reason", "n", "fraction") %in% names(stats)))
stopifnot(stats$scope[1] == "overall")
stopifnot(all(stats$reason %in% c("barcode_fail", "ambiguous_ends")))

# Combinatorial plate (bundled 384-plex): allow_single_end must abort early.
combo_path <- system.file("extdata", "index_list.csv", package = "miaoseq")
combo_err <- tryCatch(
  doDemultiplex(
    fastq = fq,
    demult_dir = file.path(tmpdir, "combo"),
    index_list = combo_path,
    n_core = 1L,
    split_reads = FALSE,
    allow_single_end = TRUE
  ),
  error = function(e) conditionMessage(e)
)
stopifnot(is.character(combo_err))
stopifnot(grepl("combinatorial|allow_single_end", combo_err, ignore.case = TRUE))

cat("All demultiplex allow_single_end tests passed.\n")
