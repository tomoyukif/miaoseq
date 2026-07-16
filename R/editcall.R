################################################################################
# Step 4: Edit-calling (per-read edit-window exact match)
################################################################################

#' Perform edit-calling from gene-assigned reads
#'
#' Consumes [doAssignGenes()] output (`gene_assign`) and per-read FASTQ
#' sequences, extracts PAM / cut-site windows, and aggregates alleles by exact
#' `read_seq` match. [doAssembleAmplicons()] is not required and is not used.
#'
#' @param gene_assign `GeneAssignResult` from [doAssignGenes()], path to
#'   `amplicon_assign/` or `gene_assignments.tsv`, or a data.frame. Required
#'   unless recoverable from `amplicon_out`.
#' @param amplicon_out Optional migration fallback: result list from
#'   [doAssembleAmplicons()] or path that contains `gene_assignments.tsv`.
#' @param pam_list Path to a CSV file with PAM / cut-site coordinates (no header).
#'   Columns: (1) gene ID matching primer / amplicon names, (2) chromosome /
#'   seqname matching `genome_fn` FASTA headers exactly, (3) PAM start (1-based),
#'   (4) optional guide ID (required when a gene has multiple rows),
#'   (5) optional strand `+` / `-` (Cas9 cut offset on genome: `+` → pam−3,
#'   `-` → pam+3; missing → cut = pam start).
#' @param genome_fn Path to the reference genome FASTA.
#' @param amplicon_fn Path to the amplicon reference FASTA from [prepAmpliconDB()].
#' @param primer_list Path to primer CSV (`primer_id`, `seq` with `_F` / `_R` suffixes).
#' @param editcall_dir Output directory for edit-calling tables.
#' @param sample_list Optional CSV mapping `index_pair_id` to plate layout (no header).
#' @param check_window Window size (bp) around the cut site for initial edit
#'   extraction (adaptive expansion starts from this window).
#' @param anchor_bp Consecutive identity (bp) required at each end of the edit
#'   window versus the reference. Ends that fail this check are expanded
#'   outward (see section 15.7 in `cursor_dev/pipeline_revise.md`).
#' @param max_expand Maximum additional bp to expand beyond `check_window` on
#'   each side while searching for anchors. Reads that cannot form anchors
#'   within this limit are discarded.
#' @param min_span_bp Minimum expected cut–cut span (bp) required before classifying
#'   a dual-guide `both_cut_excision` event (Plan A′). Default 30.
#' @param excision_tol_bp Absolute tolerance (bp) for `|del_span - expected|` when
#'   calling `both_cut_excision`. Default 20.
#' @param min_count Minimum allele count to retain (default 5).
#' @param fastq Source FASTQ path(s). Required unless `sample_fastq_dir` is set.
#' @param sample_fastq_dir Optional `by_sample/*.fq.gz` directory.
#' @param max_primer_edit Maximum edlib edit distance for primer trimming.
#' @param n_core OpenMP threads for primer trimming.
#' @return A data frame (`editcall_summary.csv` structure).
#' @import Biostrings
#' @import dplyr
#' @importFrom BiocGenerics width
#' @importFrom Biostrings readDNAStringSet writeXStringSet DNAStringSet
#'   reverseComplement matchPattern
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom pwalign pairwiseAlignment aligned pattern subject
#' @importFrom utils read.csv write.csv read.delim
#' @importFrom stats setNames aggregate
#' @importFrom methods as
#' @export
#'
doEditcall <- function(gene_assign = NULL,
                       amplicon_out = NULL,
                       pam_list,
                       genome_fn,
                       amplicon_fn,
                       primer_list,
                       editcall_dir,
                       sample_list = NULL,
                       check_window = 10L,
                       anchor_bp = 5L,
                       max_expand = 50L,
                       min_span_bp = 30L,
                       excision_tol_bp = 20L,
                       min_count = 5L,
                       fastq = NULL,
                       sample_fastq_dir = NULL,
                       max_primer_edit = 5L,
                       n_core = 1L) {
  if (is.null(pam_list) || is.null(genome_fn) || is.null(amplicon_fn) ||
      is.null(primer_list) || is.null(editcall_dir)) {
    stop("pam_list, genome_fn, amplicon_fn, primer_list, and editcall_dir are required.")
  }
  dir.create(editcall_dir, recursive = TRUE, showWarnings = FALSE)

  gene_assignments <- .resolve_gene_assignments_for_editcall(
    gene_assign = gene_assign,
    amplicon_out = amplicon_out
  )
  if (is.null(gene_assignments) || nrow(gene_assignments) < 1L) {
    stop(
      "doEditcall requires gene_assign ",
      "(GeneAssignResult / amplicon_assign / gene_assignments.tsv) ",
      "or amplicon_out containing gene_assignments."
    )
  }
  gene_assignments <- gene_assignments[
    !is.na(gene_assignments$gene_id) &
      gene_assignments$gene_id != "" &
      gene_assignments$assign_status == "assigned",
    ,
    drop = FALSE
  ]
  if (nrow(gene_assignments) < 1L) {
    warning("No gene-assigned reads found; writing empty editcall tables.")
    empty <- .empty_editcall_summary()
    write.csv(empty, file.path(editcall_dir, "editcall_summary.csv"), row.names = FALSE)
    write.csv(empty, file.path(editcall_dir, "editcall_all.csv"), row.names = FALSE)
    write.csv(empty, file.path(editcall_dir, "editcall_filtered.csv"), row.names = FALSE)
    .write_joint_editcall_outputs(
      joint_records = .empty_joint_records(),
      plan_a_read_counts = data.frame(
        index_pair_id = character(),
        gene_id = character(),
        n_plan_a_reads = integer(),
        stringsAsFactors = FALSE
      ),
      editcall_dir = editcall_dir
    )
    return(empty)
  }

  if (is.null(sample_fastq_dir) && (is.null(fastq) || !length(fastq))) {
    stop("doEditcall requires fastq or sample_fastq_dir.")
  }

  check_window <- as.integer(check_window)
  anchor_bp <- as.integer(anchor_bp)
  max_expand <- as.integer(max_expand)
  min_span_bp <- as.integer(min_span_bp)
  excision_tol_bp <- as.integer(excision_tol_bp)
  if (anchor_bp < 1L) {
    stop("anchor_bp must be >= 1.")
  }
  if (max_expand < 0L) {
    stop("max_expand must be >= 0.")
  }
  if (min_span_bp < 1L) {
    stop("min_span_bp must be >= 1.")
  }
  if (excision_tol_bp < 0L) {
    stop("excision_tol_bp must be >= 0.")
  }

  meta <- .build_editcall_metadata(
    pam_list = pam_list,
    genome_fn = genome_fn,
    amplicon_fn = amplicon_fn,
    primer_list = primer_list,
    check_window = check_window
  )
  primers <- .parse_primer_pairs_editcall(primer_list)

  read_seqs <- .load_reads_for_gene_assignments(
    gene_assignments = gene_assignments,
    fastq = fastq,
    sample_fastq_dir = sample_fastq_dir
  )
  gene_assignments$seq <- read_seqs[match(gene_assignments$read_id, names(read_seqs))]

  extracted <- .reads_to_edit_records(
    reads = gene_assignments,
    meta = meta,
    primers = primers,
    max_primer_edit = as.integer(max_primer_edit),
    n_core = as.integer(n_core),
    anchor_bp = anchor_bp,
    max_expand = max_expand,
    check_window = check_window,
    min_span_bp = min_span_bp,
    excision_tol_bp = excision_tol_bp
  )
  edit_records <- extracted$edit_records
  joint_records <- extracted$joint_records
  plan_a_read_counts <- extracted$plan_a_read_counts
  edit_records <- .filter_abnormal_edit_records(
    edit_records,
    check_window = check_window,
    max_expand = max_expand
  )

  if (nrow(edit_records) < 1L) {
    warning("No edit-site sequences could be extracted from assigned reads.")
    empty <- .empty_editcall_summary(meta$intact_seq)
    write.csv(empty, file.path(editcall_dir, "editcall_summary.csv"), row.names = FALSE)
    write.csv(empty, file.path(editcall_dir, "editcall_all.csv"), row.names = FALSE)
    write.csv(empty, file.path(editcall_dir, "editcall_filtered.csv"), row.names = FALSE)
    .write_joint_editcall_outputs(
      joint_records = joint_records,
      plan_a_read_counts = plan_a_read_counts,
      editcall_dir = editcall_dir
    )
    return(empty)
  }

  writeXStringSet(meta$intact_seq, file.path(editcall_dir, "intact_seq.fa"))
  .write_joint_editcall_outputs(
    joint_records = joint_records,
    plan_a_read_counts = plan_a_read_counts,
    editcall_dir = editcall_dir
  )
  .write_editcall_outputs(
    edit_records = edit_records,
    intact_seq = meta$intact_seq,
    editcall_dir = editcall_dir,
    sample_list = sample_list,
    min_count = as.integer(min_count)
  )
}

.resolve_gene_assignments_for_editcall <- function(gene_assign, amplicon_out) {
  if (!is.null(gene_assign)) {
    return(.as_gene_assign_df(gene_assign))
  }
  .load_gene_assignments_table(amplicon_out)
}

.load_gene_assignments_table <- function(amplicon_out) {
  if (is.null(amplicon_out)) {
    return(NULL)
  }
  if (is.list(amplicon_out) && !is.null(amplicon_out$gene_assignments) &&
      nrow(amplicon_out$gene_assignments) > 0L) {
    return(amplicon_out$gene_assignments)
  }
  out_dir <- if (is.list(amplicon_out)) amplicon_out$out_dir else amplicon_out
  ga_path <- file.path(out_dir, "gene_assignments.tsv")
  if (!file.exists(ga_path)) {
    return(NULL)
  }
  utils::read.delim(ga_path, stringsAsFactors = FALSE)
}

.load_reads_for_gene_assignments <- function(gene_assignments, fastq,
                                            sample_fastq_dir) {
  sample_ids <- sort(unique(gene_assignments$sample_id))
  id_set <- unique(gene_assignments$read_id)
  seq_map <- stats::setNames(rep(NA_character_, length(id_set)), id_set)

  if (!is.null(sample_fastq_dir)) {
    for (sid in sample_ids) {
      part <- gene_assignments[gene_assignments$sample_id == sid, , drop = FALSE]
      reads <- .read_sample_reads_if_present_editcall(sid, sample_fastq_dir)
      if (is.null(reads)) next
      hit <- match(part$read_id, reads$read_id)
      ok <- !is.na(hit)
      seq_map[part$read_id[ok]] <- reads$seq[hit[ok]]
    }
  } else {
    assn <- gene_assignments[, c("read_id", "sample_id"), drop = FALSE]
    bucketed <- bucket_fastq_assignments_cpp(
      fastq_files = as.character(fastq),
      read_ids = as.character(assn$read_id),
      sample_ids = as.character(assn$sample_id),
      target_samples = as.character(sample_ids)
    )
    for (sid in sample_ids) {
      part <- bucketed[[sid]]
      if (is.null(part)) next
      seq_map[as.character(part$read_id)] <- as.character(part$seq)
    }
  }

  seq_map
}

.read_sample_reads_if_present_editcall <- function(sample_id, sample_fastq_dir) {
  if (is.null(sample_fastq_dir)) return(NULL)
  cand <- file.path(sample_fastq_dir, paste0(sample_id, ".fq.gz"))
  if (!file.exists(cand)) {
    cand <- file.path(sample_fastq_dir, paste0(sample_id, ".fq"))
  }
  if (!file.exists(cand)) {
    return(NULL)
  }
  raw <- read_fastq_seqs_cpp(as.character(cand)[[1]])
  data.frame(
    read_id = as.character(raw$read_id),
    seq = as.character(raw$seq),
    stringsAsFactors = FALSE
  )
}

.parse_primer_pairs_editcall <- function(primer_list) {
  primers <- read.csv(primer_list, header = FALSE, stringsAsFactors = FALSE)
  f_idx <- grep("_F$", primers$V1)
  r_idx <- grep("_R$", primers$V1)
  f_primers <- primers[f_idx, , drop = FALSE]
  r_primers <- primers[r_idx, , drop = FALSE]
  f_primers$gene_id <- sub("_F$", "", f_primers$V1)
  r_primers$gene_id <- sub("_R$", "", r_primers$V1)
  pairs <- merge(
    f_primers[, c("gene_id", "V2")],
    r_primers[, c("gene_id", "V2")],
    by = "gene_id",
    suffixes = c(".f", ".r")
  )
  names(pairs) <- c("gene_id", "seq.f", "seq.r")
  pairs[order(pairs$gene_id), , drop = FALSE]
}

.reads_to_edit_records <- function(reads, meta, primers, max_primer_edit, n_core,
                                   anchor_bp = 5L, max_expand = 50L,
                                   check_window = 10L,
                                   min_span_bp = 30L, excision_tol_bp = 20L) {
  reads <- reads[!is.na(reads$seq) & nzchar(reads$seq), , drop = FALSE]
  if (nrow(reads) < 1L) {
    return(list(
      edit_records = .empty_edit_records(),
      joint_records = .empty_joint_records(),
      plan_a_read_counts = data.frame(
        index_pair_id = character(),
        gene_id = character(),
        n_plan_a_reads = integer(),
        stringsAsFactors = FALSE
      )
    ))
  }

  oriented <- reverse_complement_seqs_cpp(
    as.character(reads$seq),
    as.logical(reads$flipped)
  )
  reads$oriented_seq <- as.character(oriented)

  cache <- new.env(parent = emptyenv())
  edit_rows <- vector("list", 0L)
  joint_rows <- vector("list", 0L)
  plan_a_keys <- character()
  n_discard_expand <- 0L
  edit_i <- 0L
  joint_i <- 0L

  for (i in seq_len(nrow(reads))) {
    gid <- reads$gene_id[[i]]
    pos_idx <- which(meta$aln_pos$target == gid)
    if (length(pos_idx) < 1L) {
      next
    }
    pam_rows <- meta$aln_pos[pos_idx, , drop = FALSE]
    insert_seq <- reads$oriented_seq[[i]]

    if (!is.null(primers) && nrow(primers) > 0L && gid %in% primers$gene_id) {
      fp <- primers$seq.f[primers$gene_id == gid]
      rp <- primers$seq.r[primers$gene_id == gid]
      trimmed <- trim_amplicon_insert_cpp(
        insert_seq,
        f_primer = fp,
        r_primer = rp,
        max_edit = as.integer(max_primer_edit),
        n_core = as.integer(max(1L, n_core))
      )
      tseq <- as.character(trimmed$seq)[[1]]
      if (nzchar(tseq)) {
        insert_seq <- tseq
      }
    }

    aln_pack <- .align_insert_to_ref(
      insert_seq = insert_seq,
      gid = gid,
      meta = meta,
      cache = cache
    )
    if (is.null(aln_pack)) {
      next
    }

    local_wins <- vector("list", nrow(pam_rows))
    names(local_wins) <- pam_rows$target_gene
    for (j in seq_len(nrow(pam_rows))) {
      tg <- pam_rows$target_gene[[j]]
      local_wins[[tg]] <- .edit_seq_from_alignment(
        aln_pack = aln_pack,
        pos_row = pam_rows[j, , drop = FALSE],
        target_gene = tg,
        insert_seq = insert_seq,
        cache = cache,
        anchor_bp = as.integer(anchor_bp),
        max_expand = as.integer(max_expand)
      )
    }

    joint_events <- .classify_joint_events(
      pam_rows = pam_rows,
      local_wins = local_wins,
      ref_aln = aln_pack$ref_aln,
      query_aln = aln_pack$query_aln,
      check_window = as.integer(check_window),
      anchor_bp = as.integer(anchor_bp),
      max_expand = as.integer(max_expand),
      min_span_bp = as.integer(min_span_bp),
      excision_tol_bp = as.integer(excision_tol_bp)
    )

    excision_guide <- character()
    excision_allele <- list()
    if (length(joint_events) > 0L) {
      for (ev in joint_events) {
        joint_i <- joint_i + 1L
        joint_rows[[joint_i]] <- data.frame(
          index_pair_id = reads$sample_id[[i]],
          read_id = reads$read_id[[i]],
          gene_id = gid,
          guide_i = ev$guide_i,
          guide_j = ev$guide_j,
          target_gene_i = ev$target_gene_i,
          target_gene_j = ev$target_gene_j,
          event_class = ev$event_class,
          del_span = ev$del_span,
          expected_span = ev$expected_span,
          stringsAsFactors = FALSE
        )
        if (identical(ev$event_class, "both_cut_excision") &&
            !is.null(ev$junction)) {
          for (tg in c(ev$target_gene_i, ev$target_gene_j)) {
            excision_guide <- c(excision_guide, tg)
            excision_allele[[tg]] <- ev$junction
          }
        }
      }
    }
    excision_guide <- unique(excision_guide)

    emitted_any <- FALSE
    for (j in seq_len(nrow(pam_rows))) {
      tg <- pam_rows$target_gene[[j]]
      win <- NULL
      allele_source <- "local"
      intact <- FALSE
      if (tg %in% excision_guide) {
        win <- excision_allele[[tg]]
        allele_source <- "excision"
        intact <- FALSE
      } else {
        win <- local_wins[[tg]]
        if (!is.null(win) && !is.na(win$read_seq) && nzchar(win$read_seq)) {
          intact <- identical(win$read_seq, win$ref_seq)
          allele_source <- "local"
        } else {
          win <- NULL
        }
      }
      if (is.null(win)) {
        if (!(tg %in% excision_guide)) {
          n_discard_expand <- n_discard_expand + 1L
        }
        next
      }
      edit_i <- edit_i + 1L
      edit_rows[[edit_i]] <- data.frame(
        index_pair_id = reads$sample_id[[i]],
        target_gene = tg,
        read_seq = win$read_seq,
        ref_seq = win$ref_seq,
        count = 1L,
        intact = intact,
        allele_source = allele_source,
        stringsAsFactors = FALSE
      )
      emitted_any <- TRUE
    }
    if (emitted_any) {
      plan_a_keys <- c(
        plan_a_keys,
        paste(reads$sample_id[[i]], gid, sep = "\t")
      )
    }
  }

  if (n_discard_expand > 0L) {
    message(
      "Discarded ", n_discard_expand,
      " guide-window(s) that failed adaptive edit-window anchor search ",
      "(max_expand=", max_expand, ")."
    )
  }

  out <- if (length(edit_rows) < 1L) {
    .empty_edit_records()
  } else {
    raw <- do.call("rbind", edit_rows)
    rownames(raw) <- NULL
    agg <- aggregate(
      count ~ index_pair_id + target_gene + read_seq + ref_seq + intact + allele_source,
      data = raw,
      FUN = sum
    )
    agg[order(agg$index_pair_id, agg$target_gene, -agg$count), , drop = FALSE]
  }

  joint_out <- if (length(joint_rows) < 1L) {
    .empty_joint_records()
  } else {
    jraw <- do.call("rbind", joint_rows)
    rownames(jraw) <- NULL
    jraw
  }

  plan_a_read_counts <- if (length(plan_a_keys) < 1L) {
    data.frame(
      index_pair_id = character(),
      gene_id = character(),
      n_plan_a_reads = integer(),
      stringsAsFactors = FALSE
    )
  } else {
    tab <- table(plan_a_keys)
    parts <- strsplit(names(tab), "\t", fixed = TRUE)
    data.frame(
      index_pair_id = vapply(parts, `[[`, character(1), 1L),
      gene_id = vapply(parts, `[[`, character(1), 2L),
      n_plan_a_reads = as.integer(tab),
      stringsAsFactors = FALSE
    )
  }

  list(
    edit_records = out,
    joint_records = joint_out,
    plan_a_read_counts = plan_a_read_counts
  )
}

.empty_edit_records <- function() {
  data.frame(
    index_pair_id = character(),
    target_gene = character(),
    read_seq = character(),
    ref_seq = character(),
    count = integer(),
    intact = logical(),
    allele_source = character(),
    stringsAsFactors = FALSE
  )
}

.empty_joint_records <- function() {
  data.frame(
    index_pair_id = character(),
    read_id = character(),
    gene_id = character(),
    guide_i = character(),
    guide_j = character(),
    target_gene_i = character(),
    target_gene_j = character(),
    event_class = character(),
    del_span = integer(),
    expected_span = integer(),
    stringsAsFactors = FALSE
  )
}

.align_insert_to_ref <- function(insert_seq, gid, meta, cache) {
  if (is.na(insert_seq) || !nzchar(insert_seq)) {
    return(NULL)
  }
  if (!(gid %in% names(meta$amplicon_seq))) {
    return(NULL)
  }
  cache_key <- paste("aln", gid, insert_seq, sep = "\t")
  if (exists(cache_key, envir = cache, inherits = FALSE)) {
    return(get(cache_key, envir = cache, inherits = FALSE))
  }

  pos0 <- meta$aln_pos[meta$aln_pos$target == gid, , drop = FALSE][1, , drop = FALSE]
  amp <- as.character(meta$amplicon_seq[[gid]])
  f_len <- as.integer(pos0$f_primer_len)
  r_len <- as.integer(pos0$r_primer_len)
  ref_insert <- substr(amp, f_len + 1L, nchar(amp) - r_len)
  if (!nzchar(ref_insert)) {
    return(NULL)
  }

  aln <- pairwiseAlignment(insert_seq, ref_insert, type = "global")
  pack <- list(
    ref_aln = as.character(aligned(pwalign::subject(aln))),
    query_aln = as.character(aligned(pwalign::pattern(aln))),
    ref_insert = ref_insert
  )
  assign(cache_key, pack, envir = cache, inherits = FALSE)
  pack
}

.edit_seq_from_alignment <- function(aln_pack, pos_row, target_gene, insert_seq,
                                     cache, anchor_bp = 5L, max_expand = 50L) {
  cache_key <- paste(
    "win", target_gene, insert_seq, anchor_bp, max_expand, sep = "\t"
  )
  if (exists(cache_key, envir = cache, inherits = FALSE)) {
    return(get(cache_key, envir = cache, inherits = FALSE))
  }

  ref_len <- nchar(gsub("-", "", aln_pack$ref_aln))
  insert_aln_pos <- pos_row
  insert_aln_pos$start <- max(1L, as.integer(pos_row$start))
  insert_aln_pos$end <- min(ref_len, as.integer(pos_row$end))
  if (insert_aln_pos$end < insert_aln_pos$start) {
    assign(cache_key, NULL, envir = cache, inherits = FALSE)
    return(NULL)
  }
  win <- .extract_adaptive_edit_window(
    ref_aln = aln_pack$ref_aln,
    query_aln = aln_pack$query_aln,
    aln_pos_row = insert_aln_pos,
    anchor_bp = as.integer(anchor_bp),
    max_expand = as.integer(max_expand)
  )
  assign(cache_key, win, envir = cache, inherits = FALSE)
  win
}

.count_query_deletions <- function(ref_aln, query_aln, start_u, end_u) {
  if (end_u < start_u) {
    return(0L)
  }
  bounds <- .gapped_window_bounds(ref_aln, start_u, end_u)
  gstart <- bounds[["start"]]
  gend <- bounds[["end"]]
  if (gend < gstart || gstart < 1L || gend > nchar(ref_aln)) {
    return(0L)
  }
  r <- strsplit(substr(ref_aln, gstart, gend), "", fixed = TRUE)[[1]]
  q <- strsplit(substr(query_aln, gstart, gend), "", fixed = TRUE)[[1]]
  as.integer(sum(r != "-" & q == "-"))
}

.classify_joint_events <- function(pam_rows, local_wins, ref_aln, query_aln,
                                   check_window, anchor_bp, max_expand,
                                   min_span_bp, excision_tol_bp) {
  if (nrow(pam_rows) < 2L) {
    return(list())
  }
  ord <- order(pam_rows$cut_insert, pam_rows$guide_id, pam_rows$target_gene)
  pam_rows <- pam_rows[ord, , drop = FALSE]
  events <- vector("list", 0L)
  ei <- 0L

  for (k in seq_len(nrow(pam_rows) - 1L)) {
    row_i <- pam_rows[k, , drop = FALSE]
    row_j <- pam_rows[k + 1L, , drop = FALSE]
    tg_i <- row_i$target_gene[[1]]
    tg_j <- row_j$target_gene[[1]]
    guide_i <- if (!is.null(row_i$guide_id) && !is.na(row_i$guide_id) &&
                   nzchar(as.character(row_i$guide_id[[1]]))) {
      as.character(row_i$guide_id[[1]])
    } else {
      tg_i
    }
    guide_j <- if (!is.null(row_j$guide_id) && !is.na(row_j$guide_id) &&
                   nzchar(as.character(row_j$guide_id[[1]]))) {
      as.character(row_j$guide_id[[1]])
    } else {
      tg_j
    }

    c_i <- as.integer(row_i$cut_insert[[1]])
    c_j <- as.integer(row_j$cut_insert[[1]])
    expected <- as.integer(c_j - c_i)
    del_span <- .count_query_deletions(ref_aln, query_aln, c_i, c_j)

    left_u <- max(1L, c_i - as.integer(check_window))
    right_u <- c_j + as.integer(check_window)
    ref_len <- nchar(gsub("-", "", ref_aln))
    right_u <- min(ref_len, right_u)
    junction <- NULL
    if (right_u >= left_u) {
      junction <- .extract_adaptive_edit_window(
        ref_aln = ref_aln,
        query_aln = query_aln,
        aln_pos_row = data.frame(start = left_u, end = right_u),
        anchor_bp = as.integer(anchor_bp),
        max_expand = as.integer(max_expand)
      )
    }
    outer_ok <- !is.null(junction)

    win_i <- local_wins[[tg_i]]
    win_j <- local_wins[[tg_j]]
    i_ok <- !is.null(win_i) && !is.na(win_i$read_seq) && nzchar(win_i$read_seq)
    j_ok <- !is.null(win_j) && !is.na(win_j$read_seq) && nzchar(win_j$read_seq)
    i_wt <- i_ok && identical(win_i$read_seq, win_i$ref_seq)
    j_wt <- j_ok && identical(win_j$read_seq, win_j$ref_seq)

    event_class <- NA_character_
    if (outer_ok && expected >= as.integer(min_span_bp) &&
        abs(del_span - expected) <= as.integer(excision_tol_bp)) {
      event_class <- "both_cut_excision"
    } else if (!i_ok && !j_ok && !outer_ok) {
      next
    } else if (i_wt && j_wt) {
      event_class <- "wt"
    } else if (i_ok && !i_wt && (!j_ok || j_wt)) {
      event_class <- "g_i_only"
    } else if (j_ok && !j_wt && (!i_ok || i_wt)) {
      event_class <- "g_j_only"
    } else if (i_ok && j_ok && !i_wt && !j_wt) {
      event_class <- "both_local"
    } else if (outer_ok && (i_ok || j_ok)) {
      event_class <- "both_local"
    }
    if (is.na(event_class)) {
      next
    }

    ei <- ei + 1L
    events[[ei]] <- list(
      guide_i = guide_i,
      guide_j = guide_j,
      target_gene_i = tg_i,
      target_gene_j = tg_j,
      event_class = event_class,
      del_span = as.integer(del_span),
      expected_span = expected,
      junction = if (identical(event_class, "both_cut_excision")) junction else NULL
    )
  }
  events
}

.filter_abnormal_edit_records <- function(edit_records, check_window,
                                          max_expand = 50L) {
  if (is.null(edit_records) || nrow(edit_records) < 1L) {
    return(edit_records)
  }
  # Safety net after adaptive extraction (§15.7.6). Excision junction alleles
  # may exceed the local-window length budget.
  max_len <- 2L * as.integer(check_window) + 2L * as.integer(max_expand) + 20L
  is_excision <- if ("allele_source" %in% names(edit_records)) {
    edit_records$allele_source == "excision"
  } else {
    rep(FALSE, nrow(edit_records))
  }
  keep <- is_excision | (nchar(edit_records$read_seq) <= max_len)
  dropped <- sum(!keep)
  if (dropped > 0L) {
    warning(
      "Dropped ", dropped, " allele row(s) with read_seq longer than ",
      max_len, " bp (likely alignment failure)."
    )
  }
  edit_records[keep, , drop = FALSE]
}

.build_editcall_metadata <- function(pam_list,
                                    genome_fn,
                                    amplicon_fn,
                                    primer_list,
                                    check_window) {
  pam <- .parse_pam_list(pam_list)
  genome <- readDNAStringSet(genome_fn)
  amplicon_seq <- readDNAStringSet(amplicon_fn)
  if (is.null(names(amplicon_seq)) || any(!nzchar(names(amplicon_seq)))) {
    stop("amplicon_fn sequences must have gene ID names.")
  }

  primers <- read.csv(primer_list, header = FALSE, stringsAsFactors = FALSE)
  f_primers <- DNAStringSet(primers$V2[grep("_F$", primers$V1)])
  names(f_primers) <- sub("_F$", "", primers$V1[grep("_F$", primers$V1)])
  r_primers <- DNAStringSet(primers$V2[grep("_R$", primers$V1)])
  names(r_primers) <- sub("_R$", "", primers$V1[grep("_R$", primers$V1)])

  miss_seq <- setdiff(unique(pam$seqname), names(genome))
  if (length(miss_seq) > 0L) {
    stop(
      "pam_list seqname(s) not found in genome FASTA headers (must match ",
      "exactly, no chr/%02d rewriting): ",
      paste(miss_seq, collapse = ", "),
      call. = FALSE
    )
  }

  pam_genes <- unique(pam$gene)
  miss_amp <- setdiff(pam_genes, names(amplicon_seq))
  if (length(miss_amp) > 0L) {
    stop(
      "pam_list gene(s) missing from amplicon.fa: ",
      paste(miss_amp, collapse = ", "),
      call. = FALSE
    )
  }
  miss_fp <- setdiff(pam_genes, names(f_primers))
  miss_rp <- setdiff(pam_genes, names(r_primers))
  if (length(miss_fp) > 0L || length(miss_rp) > 0L) {
    stop(
      "pam_list gene(s) missing primer pair(s): ",
      paste(unique(c(miss_fp, miss_rp)), collapse = ", "),
      call. = FALSE
    )
  }

  locus_cache <- new.env(parent = emptyenv())
  rows <- vector("list", nrow(pam))
  intact_list <- vector("list", nrow(pam))
  intact_names <- character(nrow(pam))

  for (i in seq_len(nrow(pam))) {
    gene <- pam$gene[[i]]
    seqname <- pam$seqname[[i]]
    cut_g <- pam$cut_genome[[i]]
    tg <- pam$target_gene[[i]]

    if (!exists(gene, envir = locus_cache, inherits = FALSE)) {
      locus <- .map_amplicon_to_genome(
        amp = amplicon_seq[[gene]],
        chr_seq = genome[[seqname]],
        gene = gene,
        seqname = seqname
      )
      assign(gene, locus, envir = locus_cache, inherits = FALSE)
    } else {
      locus <- get(gene, envir = locus_cache, inherits = FALSE)
      if (!identical(locus$seqname, seqname)) {
        stop(
          "pam_list gene ", gene, " maps to multiple seqnames (",
          locus$seqname, " vs ", seqname, ").",
          call. = FALSE
        )
      }
    }

    amp_coord <- .genome_coord_to_amplicon(cut_g, locus)
    f_len <- as.integer(BiocGenerics::width(f_primers[gene]))
    r_len <- as.integer(BiocGenerics::width(r_primers[gene]))
    amp_len <- as.integer(BiocGenerics::width(amplicon_seq[gene]))
    cut_insert <- as.integer(amp_coord - f_len)
    insert_len <- amp_len - f_len - r_len
    if (is.na(cut_insert) || cut_insert < 1L || cut_insert > insert_len) {
      stop(
        "Cut for ", tg, " maps outside ref_insert (insert coord=",
        cut_insert, ", insert_len=", insert_len, ").",
        call. = FALSE
      )
    }

    win_start <- cut_insert - as.integer(check_window)
    win_end <- cut_insert + as.integer(check_window)
    rows[[i]] <- data.frame(
      start = win_start,
      end = win_end,
      cut_insert = cut_insert,
      cut_genome = cut_g,
      pam_start = pam$pam_start[[i]],
      strand = pam$strand[[i]],
      f_primer_len = f_len,
      r_primer_len = r_len,
      dist_to_end = insert_len - win_end,
      target_len = win_end - win_start,
      target = gene,
      guide_id = pam$guide_id[[i]],
      target_gene = tg,
      stringsAsFactors = FALSE
    )

    intact_start_g <- max(1L, cut_g - as.integer(check_window))
    intact_end_g <- min(as.integer(BiocGenerics::width(genome[seqname])),
                        cut_g + as.integer(check_window))
    intact_dna <- genome[[seqname]][intact_start_g:intact_end_g]
    if (identical(locus$strand, "-")) {
      intact_dna <- reverseComplement(intact_dna)
    }
    intact_list[[i]] <- intact_dna
    intact_names[[i]] <- tg
  }

  aln_pos <- do.call("rbind", rows)
  rownames(aln_pos) <- aln_pos$target_gene

  intact_seq <- Biostrings::DNAStringSet(intact_list)
  names(intact_seq) <- intact_names
  intact_seq <- intact_seq[order(names(intact_seq))]

  amp_keep <- amplicon_seq[pam_genes]
  list(
    intact_seq = intact_seq,
    amplicon_seq = amp_keep,
    aln_pos = aln_pos
  )
}

.parse_pam_list <- function(pam_list) {
  edit_site <- read.csv(pam_list, header = FALSE, stringsAsFactors = FALSE)
  if (ncol(edit_site) < 3L) {
    stop("pam_list must have at least 3 columns: gene, seqname, pam_start.")
  }
  gene <- as.character(edit_site$V1)
  seqname <- as.character(edit_site$V2)
  pam_start <- as.integer(edit_site$V3)
  if (any(is.na(pam_start))) {
    stop("pam_list column 3 (pam_start) must be integer.")
  }

  guide_id <- if (ncol(edit_site) >= 4L) {
    as.character(edit_site$V4)
  } else {
    rep(NA_character_, nrow(edit_site))
  }
  guide_id[is.na(guide_id) | !nzchar(guide_id)] <- NA_character_

  strand <- if (ncol(edit_site) >= 5L) {
    as.character(edit_site$V5)
  } else {
    rep(NA_character_, nrow(edit_site))
  }
  strand[is.na(strand) | !nzchar(strand)] <- NA_character_

  bad_strand <- !is.na(strand) & !(strand %in% c("+", "-"))
  if (any(bad_strand)) {
    stop(
      "pam_list column 5 must be '+', '-', or empty; got: ",
      paste(unique(strand[bad_strand]), collapse = ", "),
      call. = FALSE
    )
  }

  for (g in unique(gene)) {
    idx <- which(gene == g)
    if (length(idx) > 1L) {
      gids <- guide_id[idx]
      if (any(is.na(gids))) {
        stop(
          "pam_list gene ", g, " has multiple rows; guide ID (column 4) ",
          "is required for each row.",
          call. = FALSE
        )
      }
      if (any(duplicated(gids))) {
        stop(
          "pam_list gene ", g, " has duplicated guide ID(s): ",
          paste(unique(gids[duplicated(gids)]), collapse = ", "),
          call. = FALSE
        )
      }
    }
  }

  target_gene <- gene
  has_guide <- !is.na(guide_id)
  target_gene[has_guide] <- paste(gene[has_guide], guide_id[has_guide], sep = "_")
  if (any(duplicated(target_gene))) {
    stop(
      "Duplicate target_gene after pam_list parse: ",
      paste(unique(target_gene[duplicated(target_gene)]), collapse = ", "),
      call. = FALSE
    )
  }

  cut_genome <- pam_start
  cut_genome[!is.na(strand) & strand == "+"] <-
    pam_start[!is.na(strand) & strand == "+"] - 3L
  cut_genome[!is.na(strand) & strand == "-"] <-
    pam_start[!is.na(strand) & strand == "-"] + 3L

  data.frame(
    gene = gene,
    seqname = seqname,
    pam_start = pam_start,
    guide_id = guide_id,
    strand = strand,
    cut_genome = as.integer(cut_genome),
    target_gene = target_gene,
    stringsAsFactors = FALSE
  )
}

.map_amplicon_to_genome <- function(amp, chr_seq, gene, seqname) {
  amp_dna <- as(amp, "DNAString")
  hits_fwd <- Biostrings::matchPattern(amp_dna, chr_seq, max.mismatch = 0L)
  if (length(hits_fwd) > 1L) {
    stop(
      "Amplicon for gene ", gene, " has multiple exact matches on ",
      seqname, ".",
      call. = FALSE
    )
  }
  if (length(hits_fwd) == 1L) {
    return(list(
      seqname = seqname,
      start = as.integer(BiocGenerics::start(hits_fwd)),
      end = as.integer(BiocGenerics::end(hits_fwd)),
      strand = "+"
    ))
  }

  hits_rc <- Biostrings::matchPattern(
    Biostrings::reverseComplement(amp_dna),
    chr_seq,
    max.mismatch = 0L
  )
  if (length(hits_rc) > 1L) {
    stop(
      "Amplicon for gene ", gene, " has multiple exact RC matches on ",
      seqname, ".",
      call. = FALSE
    )
  }
  if (length(hits_rc) == 1L) {
    return(list(
      seqname = seqname,
      start = as.integer(BiocGenerics::start(hits_rc)),
      end = as.integer(BiocGenerics::end(hits_rc)),
      strand = "-"
    ))
  }

  stop(
    "Amplicon for gene ", gene, " not found exactly on genome seqname ",
    seqname, ".",
    call. = FALSE
  )
}

.genome_coord_to_amplicon <- function(cut_g, locus) {
  if (cut_g < locus$start || cut_g > locus$end) {
    stop(
      "Cut genome coordinate ", cut_g, " is outside amplicon locus ",
      locus$seqname, ":", locus$start, "-", locus$end, ".",
      call. = FALSE
    )
  }
  if (identical(locus$strand, "+")) {
    as.integer(cut_g - locus$start + 1L)
  } else {
    # Stored amplicon is RC of genome[locus$start:locus$end].
    as.integer(locus$end - cut_g + 1L)
  }
}

.write_joint_editcall_outputs <- function(joint_records, plan_a_read_counts,
                                          editcall_dir) {
  if (is.null(joint_records)) {
    joint_records <- .empty_joint_records()
  }
  write.csv(
    joint_records,
    file.path(editcall_dir, "editcall_joint.csv"),
    row.names = FALSE
  )

  summary_df <- .build_joint_editcall_summary(joint_records, plan_a_read_counts)
  write.csv(
    summary_df,
    file.path(editcall_dir, "editcall_joint_summary.csv"),
    row.names = FALSE
  )
  invisible(summary_df)
}

.build_joint_editcall_summary <- function(joint_records, plan_a_read_counts) {
  empty <- data.frame(
    index_pair_id = character(),
    gene_id = character(),
    guide_i = character(),
    guide_j = character(),
    n_events = integer(),
    n_wt = integer(),
    n_g_i_only = integer(),
    n_g_j_only = integer(),
    n_both_local = integer(),
    n_both_cut_excision = integer(),
    n_plan_a_reads = integer(),
    excision_rate = numeric(),
    stringsAsFactors = FALSE
  )
  if (is.null(joint_records) || nrow(joint_records) < 1L) {
    return(empty)
  }

  if (is.null(plan_a_read_counts)) {
    plan_a_read_counts <- data.frame(
      index_pair_id = character(),
      gene_id = character(),
      n_plan_a_reads = integer(),
      stringsAsFactors = FALSE
    )
  }

  keys <- unique(joint_records[, c("index_pair_id", "gene_id", "guide_i", "guide_j")])
  out <- lapply(seq_len(nrow(keys)), function(i) {
    key <- keys[i, , drop = FALSE]
    sub <- joint_records[
      joint_records$index_pair_id == key$index_pair_id &
        joint_records$gene_id == key$gene_id &
        joint_records$guide_i == key$guide_i &
        joint_records$guide_j == key$guide_j,
      ,
      drop = FALSE
    ]
    n_pa <- plan_a_read_counts$n_plan_a_reads[
      plan_a_read_counts$index_pair_id == key$index_pair_id &
        plan_a_read_counts$gene_id == key$gene_id
    ]
    if (length(n_pa) < 1L || is.na(n_pa[[1]])) {
      n_pa <- length(unique(sub$read_id))
    } else {
      n_pa <- as.integer(n_pa[[1]])
    }
    n_ex <- sum(sub$event_class == "both_cut_excision")
    data.frame(
      index_pair_id = key$index_pair_id,
      gene_id = key$gene_id,
      guide_i = key$guide_i,
      guide_j = key$guide_j,
      n_events = nrow(sub),
      n_wt = sum(sub$event_class == "wt"),
      n_g_i_only = sum(sub$event_class == "g_i_only"),
      n_g_j_only = sum(sub$event_class == "g_j_only"),
      n_both_local = sum(sub$event_class == "both_local"),
      n_both_cut_excision = n_ex,
      n_plan_a_reads = n_pa,
      excision_rate = if (n_pa > 0L) n_ex / n_pa else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  do.call("rbind", out)
}

#' Map inclusive ungapped ref coordinates to gapped alignment [start, end].
.gapped_window_bounds <- function(ref_aln, start_ungapped, end_ungapped) {
  x <- ref_aln
  target_len <- as.integer(end_ungapped) - as.integer(start_ungapped)
  ins <- unlist(gregexpr("-", substr(x, 1L, start_ungapped - 1L)))
  if (ins[1] < 0L) {
    n_ins <- 0L
  } else {
    n_ins <- length(ins)
  }
  target_start <- start_ungapped + n_ins
  target_rest <- substr(x, target_start, nchar(x))
  target_end <- target_len
  detected_ins <- 0L
  while (TRUE) {
    ins <- unlist(gregexpr("-", substr(target_rest, 1L, target_end)))
    if (ins[1] < 0L) {
      n_ins <- 0L
    } else {
      n_ins <- length(ins)
    }
    if (n_ins > detected_ins) {
      added_ins <- n_ins - detected_ins
      target_end <- target_end + added_ins
      detected_ins <- n_ins
    } else {
      break
    }
  }
  target_end <- target_start + target_end
  c(start = as.integer(target_start), end = as.integer(target_end))
}

.alignment_end_has_anchor <- function(ref_aln, query_aln, gstart, gend,
                                      anchor_bp, side = c("left", "right")) {
  side <- match.arg(side)
  if (gend < gstart || gstart < 1L || gend > nchar(ref_aln)) {
    return(FALSE)
  }
  ref_chars <- strsplit(substr(ref_aln, gstart, gend), "", fixed = TRUE)[[1]]
  qry_chars <- strsplit(substr(query_aln, gstart, gend), "", fixed = TRUE)[[1]]
  n <- length(ref_chars)
  if (n < anchor_bp) {
    return(FALSE)
  }
  idx <- if (side == "left") {
    seq_len(anchor_bp)
  } else {
    seq.int(n - anchor_bp + 1L, n)
  }
  for (i in idx) {
    rc <- ref_chars[[i]]
    qc <- qry_chars[[i]]
    if (rc == "-" || qc == "-" || rc != qc) {
      return(FALSE)
    }
  }
  TRUE
}

.extract_adaptive_edit_window <- function(ref_aln, query_aln, aln_pos_row,
                                          anchor_bp = 5L, max_expand = 50L) {
  left_u <- as.integer(aln_pos_row$start)
  right_u <- as.integer(aln_pos_row$end)
  if (is.na(left_u) || is.na(right_u) || right_u < left_u) {
    return(NULL)
  }

  ref_ungapped_len <- nchar(gsub("-", "", ref_aln))
  min_left <- max(1L, left_u - as.integer(max_expand))
  max_right <- min(ref_ungapped_len, right_u + as.integer(max_expand))

  repeat {
    bounds <- .gapped_window_bounds(ref_aln, left_u, right_u)
    gstart <- bounds[["start"]]
    gend <- bounds[["end"]]
    if (gend < gstart || gend > nchar(query_aln) || gstart < 1L) {
      return(NULL)
    }

    left_ok <- .alignment_end_has_anchor(
      ref_aln, query_aln, gstart, gend, anchor_bp, side = "left"
    )
    right_ok <- .alignment_end_has_anchor(
      ref_aln, query_aln, gstart, gend, anchor_bp, side = "right"
    )
    if (left_ok && right_ok) {
      read_seq <- substr(query_aln, gstart, gend)
      ref_seq <- substr(ref_aln, gstart, gend)
      if (!nzchar(gsub("-", "", read_seq))) {
        return(NULL)
      }
      return(list(read_seq = read_seq, ref_seq = ref_seq))
    }

    expanded <- FALSE
    if (!left_ok) {
      if (left_u <= min_left) {
        return(NULL)
      }
      left_u <- left_u - 1L
      expanded <- TRUE
    }
    if (!right_ok) {
      if (right_u >= max_right) {
        return(NULL)
      }
      right_u <- right_u + 1L
      expanded <- TRUE
    }
    if (!expanded) {
      return(NULL)
    }
  }
}

.extract_edit_seq_from_gapped <- function(ref_aln, query_aln, aln_pos_row) {
  # Compatibility wrapper: fixed window (no adaptive expansion).
  bounds <- .gapped_window_bounds(
    ref_aln,
    start_ungapped = aln_pos_row$start,
    end_ungapped = aln_pos_row$end
  )
  read_seq <- substr(query_aln, bounds[["start"]], bounds[["end"]])
  if (!nzchar(gsub("-", "", read_seq))) {
    return(NA_character_)
  }
  read_seq
}

.write_editcall_outputs <- function(edit_records,
                                    intact_seq,
                                    editcall_dir,
                                    sample_list,
                                    min_count) {
  write.csv(edit_records, file.path(editcall_dir, "editcall_all.csv"), row.names = FALSE)

  edit_df <- edit_records[edit_records$count >= min_count, , drop = FALSE]
  edit_df_filtered <- do.call("rbind", tapply(seq_len(nrow(edit_df)), edit_df$index_pair_id, function(i) {
    i_df <- edit_df[i, , drop = FALSE]
    do.call("rbind", tapply(seq_len(nrow(i_df)), i_df$target_gene, function(j) {
      ij_df <- i_df[j, , drop = FALSE]
      top_count <- max(ij_df$count)
      top_df <- ij_df[ij_df$count > top_count / 2, , drop = FALSE]
      top_df$vs_intact_ratio <- 0
      top_df$intact_count <- 0
      if (any(top_df$intact)) {
        top_df$intact_count <- top_df$count[top_df$intact]
      }
      top_df$vs_intact_ratio <- top_df$count / (top_df$count + top_df$intact_count)
      top_df$vs_intact_ratio[top_df$intact] <- NA
      top_df
    }))
  }))
  if (!is.null(edit_df_filtered)) {
    rownames(edit_df_filtered) <- NULL
    write.csv(edit_df_filtered, file.path(editcall_dir, "editcall_filtered.csv"), row.names = FALSE)
  } else {
    edit_df_filtered <- edit_df[0, , drop = FALSE]
    write.csv(edit_df_filtered, file.path(editcall_dir, "editcall_filtered.csv"), row.names = FALSE)
  }

  editcall_out <- .build_editcall_summary(edit_df_filtered, intact_seq)
  total_reads_per_sample <- tapply(edit_df_filtered$count, edit_df_filtered$index_pair_id, sum)
  total_reads_df <- data.frame(
    index_pair_id = names(total_reads_per_sample),
    total_reads = as.numeric(total_reads_per_sample),
    stringsAsFactors = FALSE
  )

  if (!is.null(sample_list)) {
    sample_info <- .read_plate_layout(sample_list)
    editcall_out <- left_join(editcall_out, sample_info, by = "index_pair_id")
    editcall_out <- left_join(editcall_out, total_reads_df, by = "index_pair_id")
    gene_cols <- setdiff(
      names(editcall_out),
      c("index_pair_id", "data_type", "sample_name", "plate_id", "row_id", "col_id", "total_reads")
    )
    editcall_out <- editcall_out[, c(
      "index_pair_id", "data_type", gene_cols,
      "sample_name", "plate_id", "row_id", "col_id", "total_reads"
    )]
  } else {
    editcall_out <- left_join(editcall_out, total_reads_df, by = "index_pair_id")
  }

  write.csv(editcall_out, file.path(editcall_dir, "editcall_summary.csv"), row.names = FALSE)
  editcall_out
}

.build_editcall_summary <- function(edit_df_filtered, intact_seq) {
  if (nrow(edit_df_filtered) < 1L) {
    return(.empty_editcall_summary(intact_seq))
  }

  editcall_out <- do.call("rbind", tapply(seq_len(nrow(edit_df_filtered)), edit_df_filtered$index_pair_id, function(i) {
    i_df <- edit_df_filtered[i, , drop = FALSE]
    out <- data.frame(target_gene = names(intact_seq))
    i_out <- do.call("rbind", tapply(seq_len(nrow(i_df)), i_df$target_gene, function(j) {
      ij_df <- i_df[j, , drop = FALSE]
      ij_target_gene <- ij_df$target_gene[1]
      genotype <- rep("ref", nrow(ij_df))
      if ("ref_seq" %in% names(ij_df) && !all(is.na(ij_df$ref_seq))) {
        n_del <- vapply(ij_df$read_seq, .count_gap_bp, integer(1))
        n_ins <- vapply(ij_df$ref_seq, .count_gap_bp, integer(1))
      } else {
        # Legacy path without paired gapped ref_seq: fall back to read gaps only.
        n_del <- vapply(ij_df$read_seq, .count_gap_bp, integer(1))
        n_ins <- rep(0L, nrow(ij_df))
      }
      genotype[!ij_df$intact] <- "sub"
      genotype[n_del > 0L & n_ins > 0L] <- "indel"
      genotype[n_del > 0L & n_ins == 0L] <- "del"
      genotype[n_del == 0L & n_ins > 0L] <- "ins"
      genotype_order <- factor(genotype, levels = c("ref", "sub", "ins", "del", "indel"))
      genotype_order <- order(as.numeric(genotype_order))
      pure_del <- genotype == "del"
      if (any(pure_del)) {
        genotype[pure_del] <- paste0("del", n_del[pure_del])
      }
      pure_ins <- genotype == "ins"
      if (any(pure_ins)) {
        genotype[pure_ins] <- paste0("ins", n_ins[pure_ins])
      }
      both <- genotype == "indel"
      if (any(both)) {
        genotype[both] <- paste0("indel", n_del[both], "-", n_ins[both])
      }
      genotype <- genotype[genotype_order]
      which_ref <- genotype == "ref"
      genotype <- paste(genotype, collapse = "/")
      alt_patterns <- ij_df$read_seq[genotype_order]
      alt_patterns <- paste(alt_patterns[!which_ref], collapse = "/")
      count <- paste0(
        sum(ij_df$count),
        "(",
        paste(ij_df$count[genotype_order], collapse = ","),
        ")"
      )
      data.frame(
        target_gene = ij_target_gene,
        genotype = genotype,
        alt_patterns = alt_patterns,
        count = count,
        stringsAsFactors = FALSE
      )
    }))
    out <- left_join(out, i_out, "target_gene")
    out <- t(out)
    colnames(out) <- out[1, ]
    out <- out[-1, , drop = FALSE]
    cbind(
      index_pair_id = i_df$index_pair_id[1],
      data_type = c("genotype", "seq", "count"),
      out
    )
  }))
  as.data.frame(editcall_out, stringsAsFactors = FALSE)
}

.empty_editcall_summary <- function(intact_seq = NULL) {
  genes <- if (is.null(intact_seq)) character() else names(intact_seq)
  out <- data.frame(
    index_pair_id = character(),
    data_type = character(),
    stringsAsFactors = FALSE
  )
  for (g in genes) {
    out[[g]] <- character()
  }
  out
}

#' Count gap bases (`-`) in a gapped alignment string.
#' @keywords internal
.count_gap_bp <- function(x) {
  x <- as.character(x)[[1]]
  if (is.na(x) || !nzchar(x)) return(0L)
  as.integer(nchar(gsub("[^-]", "", x)))
}
