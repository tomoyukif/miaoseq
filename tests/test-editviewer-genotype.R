#!/usr/bin/env Rscript
# Unit tests for editViewer genotype classification.

suppressPackageStartupMessages({
  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("devtools required")
  }
  devtools::load_all(".", quiet = TRUE)
})

parse <- miaoseq:::.parse_allele_label
classify <- miaoseq:::.classify_genotype_edit_eval

p <- parse("ref")
stopifnot(isTRUE(p$is_ref), !isTRUE(p$is_frameshift))

p <- parse("sub")
stopifnot(!isTRUE(p$is_ref), !isTRUE(p$is_frameshift))

p <- parse("del3")
stopifnot(identical(p$n_del, 3L), !isTRUE(p$is_frameshift))

p <- parse("del1")
stopifnot(isTRUE(p$is_frameshift))

p <- parse("ins6")
stopifnot(identical(p$n_ins, 6L), !isTRUE(p$is_frameshift))

p <- parse("indel5-3")
stopifnot(identical(p$n_del, 5L), identical(p$n_ins, 3L), isTRUE(p$is_frameshift))

p <- parse("indel6-3")
stopifnot(identical(p$n_del, 6L), identical(p$n_ins, 3L), !isTRUE(p$is_frameshift))

p <- parse("indel5-2")
stopifnot(!isTRUE(p$is_frameshift))  # |5-2|=3

stopifnot(identical(classify("ref"), "ref"))
stopifnot(identical(classify("sub"), "alt_inframe_homo"))  # SNP ≠ frameshift
stopifnot(identical(classify("del1"), "alt"))
stopifnot(identical(classify("del3"), "alt_inframe_homo"))
stopifnot(identical(classify("del1/del3"), "alt_inframe_het"))
stopifnot(identical(classify("ref/sub"), "het_inframe"))
stopifnot(identical(classify("ref/del1"), "het"))
stopifnot(identical(classify("ref/indel6-3"), "het_inframe"))
stopifnot(identical(classify("indel5-3"), "alt"))
stopifnot(identical(classify("indel6-3"), "alt_inframe_homo"))
stopifnot(is.na(classify(NA_character_)))

cat("All editViewer genotype classification tests passed.\n")
