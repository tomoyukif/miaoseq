# Phase K Assemble API / backend invariants
# Skip integration if vsearch/abpoa are not on PATH.

library(testthat)

test_that("doAssembleAmplicons rejects non-vsearch backends", {
  expect_error(
    doAssembleAmplicons(
      gene_assign = data.frame(),
      out_dir = tempfile(),
      cluster_backend = "mmseqs"
    ),
    "vsearch"
  )
  expect_error(
    doAssembleAmplicons(
      gene_assign = data.frame(),
      out_dir = tempfile(),
      cluster_backend = "internal"
    ),
    "vsearch"
  )
})

test_that("Assemble formals match Phase K defaults", {
  f <- formals(doAssembleAmplicons)
  expect_identical(eval(f$cluster_backend), "vsearch")
  expect_identical(eval(f$max_clusters), Inf)
  expect_equal(eval(f$min_cluster_identity), 0.95)
  expect_null(eval(f$min_cluster_purity))
  expect_false("max_cluster_edit" %in% names(f))
})

test_that("polish/racon and mmseqs cluster helpers are gone", {
  expect_false(exists(".polish_via_racon", where = asNamespace("miaoseq"),
                       inherits = FALSE))
  expect_false(exists(".cluster_ids_mmseqs", where = asNamespace("miaoseq"),
                       inherits = FALSE))
})

has_tools <- nzchar(Sys.which("vsearch")) && nzchar(Sys.which("abpoa"))

test_that("tiny vsearch+abpoa assemble run (skipped without binaries)", {
  skip_if_not(has_tools, "vsearch/abpoa not on PATH")

  # Two near-identical inserts + one distant outlier
  seq_a <- paste(rep("ACGT", 50), collapse = "")
  seq_b <- paste0(substr(seq_a, 1, nchar(seq_a) - 1), "A")
  seq_c <- paste(rep("TGCA", 50), collapse = "")
  ga <- data.frame(
    read_id = c("r1", "r2", "r3", "r4", "r5", "r6"),
    sample_id = "s1",
    gene_id = "g1",
    assign_status = "assigned",
    flipped = FALSE,
    oriented_seq = c(seq_a, seq_a, seq_b, seq_b, seq_b, seq_c),
    stringsAsFactors = FALSE
  )
  out <- tempfile("assemble_phase_k_")
  amp <- doAssembleAmplicons(
    gene_assign = ga,
    out_dir = out,
    primer_list = NULL,
    min_reads = 2L,
    min_cluster_reads = 2L,
    max_clusters = Inf,
    min_cluster_identity = 0.90,
    overwrite = TRUE,
    n_core = 1L
  )
  expect_true(nrow(amp$table) >= 1L)
  expect_true(all(amp$table$method %in% c("abpoa", "greedy")))
  expect_true("fraction_bucket" %in% names(amp$table))
  expect_false("fraction_sample" %in% names(amp$table))
  cons <- file.path(out, "s1", "consensus.fasta")
  clus <- file.path(out, "s1", "clusters.fasta")
  expect_true(file.exists(cons))
  expect_true(file.exists(clus))
  stats <- readLines(file.path(out, "run_stats.txt"))
  expect_true(any(grepl("^r_version=", stats)))
  expect_true(any(grepl("^min_cluster_identity=", stats)))
  expect_true(any(grepl("^version_vsearch=", stats)))
  expect_true(any(grepl("^version_abpoa=", stats)))
})
