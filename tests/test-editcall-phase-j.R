# Phase J: pam_list mapping, multi-guide Plan A / A′
library(testthat)
library(miaoseq)
library(Biostrings)

tmp <- tempfile("miaoseq_phase_j_")
dir.create(tmp)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

# Synthetic genome / amplicon / primers
# Amplicon = Fprim(10) + insert(200) + Rprim(10); cuts at insert 60 and 140
fprim <- "AAAAAAAAAA"
rprim <- "TTTTTTTTTT"
insert <- paste0(
  paste(rep("C", 50), collapse = ""),
  "G",  # cut1 at insert 51 when pam without strand at genome coord
  paste(rep("G", 89), collapse = ""),
  "T",  # cut2
  paste(rep("A", 59), collapse = "")
)
stopifnot(nchar(insert) == 200L)
amp <- paste0(fprim, insert, rprim)
flank_l <- paste(rep("N", 20), collapse = "")
flank_r <- paste(rep("N", 20), collapse = "")
chr_seq <- paste0(flank_l, amp, flank_r)
# cut genome coords: amp starts at 21; cut1 insert 51 → amp 61 → genome 81
# cut2 insert 141 → amp 151 → genome 171

genome_fn <- file.path(tmp, "genome.fa")
amp_fn <- file.path(tmp, "amplicon.fa")
primer_fn <- file.path(tmp, "primers.csv")
pam_fn <- file.path(tmp, "pam.csv")
pam_bad_fn <- file.path(tmp, "pam_bad.csv")
pam_multi_fn <- file.path(tmp, "pam_multi.csv")

writeXStringSet(DNAStringSet(c(chr01 = chr_seq)), genome_fn)
writeXStringSet(DNAStringSet(c(GENE1 = amp)), amp_fn)
write.table(
  data.frame(V1 = c("GENE1_F", "GENE1_R"), V2 = c(fprim, rprim), stringsAsFactors = FALSE),
  primer_fn,
  sep = ",",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

# Single PAM without strand → cut = pam_start
writeLines("GENE1,chr01,81", pam_fn)
# Wrong seqname
writeLines("GENE1,chr99,81", pam_bad_fn)
# Dual guide
writeLines(c(
  "GENE1,chr01,81,g1",
  "GENE1,chr01,171,g2"
), pam_multi_fn)

test_that("seqname mismatch stops", {
  expect_error(
    miaoseq:::.build_editcall_metadata(
      pam_list = pam_bad_fn,
      genome_fn = genome_fn,
      amplicon_fn = amp_fn,
      primer_list = primer_fn,
      check_window = 5L
    ),
    regexp = "seqname"
  )
})

test_that("single PAM metadata maps cut to insert", {
  meta <- miaoseq:::.build_editcall_metadata(
    pam_list = pam_fn,
    genome_fn = genome_fn,
    amplicon_fn = amp_fn,
    primer_list = primer_fn,
    check_window = 5L
  )
  expect_equal(meta$aln_pos$target_gene, "GENE1")
  expect_equal(as.integer(meta$aln_pos$cut_insert), 51L)
  expect_equal(as.integer(meta$aln_pos$start), 46L)
  expect_equal(as.integer(meta$aln_pos$end), 56L)
})

test_that("multi-PAM keeps both target_genes", {
  meta <- miaoseq:::.build_editcall_metadata(
    pam_list = pam_multi_fn,
    genome_fn = genome_fn,
    amplicon_fn = amp_fn,
    primer_list = primer_fn,
    check_window = 5L
  )
  expect_equal(sort(meta$aln_pos$target_gene), c("GENE1_g1", "GENE1_g2"))
  expect_equal(sort(as.integer(meta$aln_pos$cut_insert)), c(51L, 141L))
})

test_that("joint excision classifies large inter-cut deletion", {
  # Perfect alignment with query gaps spanning cut1..cut2 (90 bp expected)
  ref_insert <- insert
  # Build gapped alignment strings: delete ref bases 52..141 in query
  # Positions 1..51 kept, then gaps for 52..141, then 142..200 kept
  q <- paste0(
    substr(ref_insert, 1, 51),
    substr(ref_insert, 142, 200)
  )
  # Global align
  aln <- pwalign::pairwiseAlignment(q, ref_insert, type = "global")
  ref_aln <- as.character(pwalign::aligned(pwalign::subject(aln)))
  query_aln <- as.character(pwalign::aligned(pwalign::pattern(aln)))

  pam_rows <- data.frame(
    start = c(46L, 136L),
    end = c(56L, 146L),
    cut_insert = c(51L, 141L),
    guide_id = c("g1", "g2"),
    target_gene = c("GENE1_g1", "GENE1_g2"),
    stringsAsFactors = FALSE
  )
  local_wins <- list(
    GENE1_g1 = NULL,
    GENE1_g2 = NULL
  )
  events <- miaoseq:::.classify_joint_events(
    pam_rows = pam_rows,
    local_wins = local_wins,
    ref_aln = ref_aln,
    query_aln = query_aln,
    check_window = 5L,
    anchor_bp = 3L,
    max_expand = 40L,
    min_span_bp = 30L,
    excision_tol_bp = 20L
  )
  expect_true(length(events) >= 1L)
  expect_equal(events[[1]]$event_class, "both_cut_excision")
  expect_true(abs(events[[1]]$del_span - events[[1]]$expected_span) <= 20L)
})

test_that("min_span_bp blocks excision class", {
  ref_insert <- paste(rep("A", 80), collapse = "")
  q <- ref_insert
  aln <- pwalign::pairwiseAlignment(q, ref_insert, type = "global")
  ref_aln <- as.character(pwalign::aligned(pwalign::subject(aln)))
  query_aln <- as.character(pwalign::aligned(pwalign::pattern(aln)))
  pam_rows <- data.frame(
    start = c(20L, 35L),
    end = c(30L, 45L),
    cut_insert = c(25L, 40L),
    guide_id = c("g1", "g2"),
    target_gene = c("G_g1", "G_g2"),
    stringsAsFactors = FALSE
  )
  local_wins <- list(
    G_g1 = list(read_seq = "AAAAAAAAAAA", ref_seq = "AAAAAAAAAAA"),
    G_g2 = list(read_seq = "AAAAAAAAAAA", ref_seq = "AAAAAAAAAAA")
  )
  events <- miaoseq:::.classify_joint_events(
    pam_rows = pam_rows,
    local_wins = local_wins,
    ref_aln = ref_aln,
    query_aln = query_aln,
    check_window = 5L,
    anchor_bp = 3L,
    max_expand = 10L,
    min_span_bp = 30L,
    excision_tol_bp = 20L
  )
  expect_true(length(events) >= 1L)
  expect_false(identical(events[[1]]$event_class, "both_cut_excision"))
  expect_equal(events[[1]]$event_class, "wt")
})

cat("Phase J unit checks completed via testthat expectations.\n")
