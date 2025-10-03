# Optimized miaoseq functions using Rcpp for performance
# This file contains optimized versions of critical functions

# Load Rcpp functions
if (requireNamespace("Rcpp", quietly = TRUE)) {
    Rcpp::sourceCpp("src/fast_string_ops.cpp")
} else {
    warning("Rcpp not available. Optimized functions will not work.")
}

# Declare C++ functions to avoid linter warnings
fast_gregexpr <- function(text, pattern) { stop("C++ function not loaded") }
fast_substr_multiple <- function(text, starts, ends) { stop("C++ function not loaded") }
fast_pattern_count <- function(text, pattern) { stop("C++ function not loaded") }
fast_gsub <- function(text, pattern, replacement) { stop("C++ function not loaded") }
fast_seq_length_no_gaps <- function(sequences) { stop("C++ function not loaded") }
fast_has_deletions <- function(sequences) { stop("C++ function not loaded") }
fast_has_insertions <- function(sequences, references) { stop("C++ function not loaded") }
fast_sequence_table <- function(sequences, counts) { stop("C++ function not loaded") }
fast_alignment_stats <- function(query_seqs, subject_seqs, starts, ends) { stop("C++ function not loaded") }

################################################################################
# Optimized string processing functions
################################################################################

#' Fast pattern matching and extraction
#' @param text Character vector of text to search
#' @param pattern Pattern to search for
#' @return Character vector with match positions
fast_gregexpr_optimized <- function(text, pattern) {
    fast_gregexpr(text, pattern)
}

#' Fast substring extraction
#' @param text Character vector
#' @param starts Start positions
#' @param ends End positions
#' @return Character vector of extracted substrings
fast_substr_optimized <- function(text, starts, ends) {
    fast_substr_multiple(text, starts, ends)
}

#' Fast pattern counting
#' @param text Character vector
#' @param pattern Pattern to count
#' @return Integer vector of counts
fast_pattern_count_optimized <- function(text, pattern) {
    fast_pattern_count(text, pattern)
}

#' Fast string replacement
#' @param text Character vector
#' @param pattern Pattern to replace
#' @param replacement Replacement string
#' @return Character vector with replacements
fast_gsub_optimized <- function(text, pattern, replacement) {
    fast_gsub(text, pattern, replacement)
}

################################################################################
# Optimized sequence processing functions
################################################################################

#' Fast sequence length calculation excluding gaps
#' @param sequences Character vector of sequences
#' @return Integer vector of lengths
fast_seq_length_no_gaps_optimized <- function(sequences) {
    fast_seq_length_no_gaps(sequences)
}

#' Fast deletion detection
#' @param sequences Character vector of sequences
#' @return Logical vector indicating deletions
fast_has_deletions_optimized <- function(sequences) {
    fast_has_deletions(sequences)
}

#' Fast insertion detection
#' @param sequences Character vector of sequences
#' @param references Character vector of reference sequences
#' @return Logical vector indicating insertions
fast_has_insertions_optimized <- function(sequences, references) {
    fast_has_insertions(sequences, references)
}

#' Fast sequence counting
#' @param sequences Character vector of sequences
#' @param counts Integer vector of counts
#' @return List with sequence names and counts
fast_sequence_table_optimized <- function(sequences, counts) {
    fast_sequence_table(sequences, counts)
}

################################################################################
# Optimized alignment processing
################################################################################

#' Fast alignment statistics calculation
#' @param query_seqs Query sequences
#' @param subject_seqs Subject sequences
#' @param starts Start positions
#' @param ends End positions
#' @return DataFrame with alignment statistics
fast_alignment_stats_optimized <- function(query_seqs, subject_seqs, starts, ends) {
    fast_alignment_stats(query_seqs, subject_seqs, starts, ends)
}

################################################################################
# Optimized version of doEditcall function
################################################################################

#' Optimized edit-calling function
#' @param demult_out Demultiplexing results
#' @param align_out Alignment results
#' @param editcall_dir Edit-calling output directory
#' @return Edit-calling results
doEditcall_optimized <- function(demult_out, align_out, editcall_dir) {
    # Create demultiplexing data frame
    demult_df <- data.frame(read_name = demult_out$sseqid,
                            i7_index = demult_out$qseqid.f,
                            i5_index = demult_out$qseqid.r,
                            index_pair_id = demult_out$index_pair_id)
    demult_df <- unique(demult_df)
    align_out <- left_join(align_out, demult_df, "read_name")
    
    intact_seq <- attributes(align_out)$intact_seq
    
    # Use optimized sequence processing
    edit_df <- tapply(seq_len(nrow(align_out)), align_out$index_pair_id, function(i) {
        i_align_out <- align_out[i, ]
        i_out <- tapply(seq_len(nrow(i_align_out)), i_align_out$target_gene, function(j) {
            ij_align_out <- i_align_out[j, ]
            
            # Use optimized sequence counting
            seq_counts <- fast_sequence_table_optimized(ij_align_out$read_seq, 
                                                       rep(1, nrow(ij_align_out)))
            
            ij_out <- data.frame(target_gene = ij_align_out$target_gene[1],
                                 read_seq = seq_counts$names,
                                 count = seq_counts$values)
            
            # Use optimized sequence comparison
            ij_out$intact <- ij_out$read_seq %in% intact_seq[names(intact_seq) %in% ij_align_out$target_gene[1]]
            ij_out <- ij_out[order(ij_out$count, decreasing = TRUE), ]
            return(ij_out)
        })
        i_out <- do.call("rbind", i_out)
        i_out <- cbind(index_pair_id = i_align_out$index_pair_id[1], i_out)
        return(i_out)
    })
    edit_df <- do.call("rbind", edit_df)
    write.csv(edit_df, file.path(editcall_dir, "editcall_all.csv"), row.names = FALSE)
    
    # Filtered results with optimized processing
    edit_df_filtered <- tapply(seq_len(nrow(edit_df)), edit_df$index_pair_id, function(i) {
        i_df <- edit_df[i, ]
        i_out <- tapply(seq_len(nrow(i_df)), i_df$target_gene, function(j) {
            ij_df <- i_df[j, ]
            top_count <- max(ij_df$count)
            top_df <- ij_df[ij_df$count > top_count / 2, ]
            top_df$vs_intact_ratio <- 0
            top_df$intact_count <- 0
            
            if (any(top_df$intact)) {
                top_df$intact_count <- top_df$count[top_df$intact]
            }
            top_df$vs_intact_ratio <- top_df$count / (top_df$count + top_df$intact_count)
            top_df$vs_intact_ratio[top_df$intact] <- NA
            return(top_df)
        })
        i_out <- do.call("rbind", i_out)
        return(i_out)
    })
    edit_df_filtered <- do.call("rbind", edit_df_filtered)
    write.csv(edit_df_filtered, file.path(editcall_dir, "editcall_filtered.csv"), row.names = FALSE)
    
    # Final edit-calling with optimized sequence analysis
    editcall_out <- tapply(seq_len(nrow(edit_df_filtered)), edit_df_filtered$index_pair_id, function(i) {
        i_df <- edit_df_filtered[i, ]
        out <- data.frame(target_gene = names(intact_seq))
        i_out <- tapply(seq_len(nrow(i_df)), i_df$target_gene, function(j) {
            ij_df <- i_df[j, ]
            ij_target_gene <- ij_df$target_gene[1]
            genotype <- rep("ref", nrow(ij_df))
            ij_intact_seq <- intact_seq[names(intact_seq) %in% ij_target_gene]
            
            # Use optimized sequence analysis
            ij_del <- fast_has_deletions_optimized(ij_df$read_seq)
            ij_ins <- fast_has_insertions_optimized(ij_df$read_seq, ij_intact_seq)
            
            genotype[!ij_df$intact] <- "sub"
            genotype[ij_del & ij_ins] <- "indel"
            genotype[ij_del & !ij_ins] <- "del"
            genotype[!ij_del & ij_ins] <- "ins"
            
            genotype_order <- factor(genotype, levels = c("ref", "sub", "ins", "del", "indel"))
            genotype_order <- order(as.numeric(genotype_order))
            
            # Use optimized pattern counting for size determination
            if (any(ij_del & !ij_ins)) {
                n_ij_del <- fast_pattern_count_optimized(ij_df$read_seq[ij_del & !ij_ins], "-")
                genotype[ij_del & !ij_ins] <- paste0(genotype[ij_del & !ij_ins], n_ij_del)
            }
            
            if (any(!ij_del & ij_ins)) {
                n_ij_ins <- fast_seq_length_no_gaps_optimized(ij_df$read_seq[!ij_del & ij_ins]) - 
                           nchar(ij_intact_seq)
                genotype[!ij_del & ij_ins] <- paste0(genotype[!ij_del & ij_ins], n_ij_ins)
            }
            
            if (any(ij_del & ij_ins)) {
                n_ij_del <- fast_pattern_count_optimized(ij_df$read_seq[ij_del & ij_ins], "-")
                n_ij_ins <- fast_seq_length_no_gaps_optimized(ij_df$read_seq[ij_del & ij_ins]) - 
                           nchar(ij_intact_seq)
                genotype[ij_del & ij_ins] <- paste0(genotype[ij_del & ij_ins], n_ij_ins, "-", n_ij_ins)
            }
            
            genotype <- genotype[genotype_order]
            which_ref <- genotype == "ref"
            genotype <- paste(genotype, collapse = "/")
            alt_patterns <- ij_df$read_seq[genotype_order]
            alt_patterns <- paste(alt_patterns[!which_ref], collapse = "/")
            count <- paste0(sum(ij_df$count),
                            "(",
                            paste(ij_df$count[genotype_order], collapse = ","),
                            ")")
            ij_out <- data.frame(target_gene = ij_target_gene,
                                 genotype = genotype,
                                 alt_patterns = alt_patterns,
                                 count = count)
            return(ij_out)
        })
        i_out <- do.call("rbind", i_out)
        out <- left_join(out, i_out, "target_gene")
        out <- t(out)
        colnames(out) <- out[1, ]
        out <- out[-1, ]
        out <- cbind(index_pair_id = i_df$index_pair_id[1],
                     sample_id = i_df$sample_id[1],
                     data_type = c("genotype", "seq", "count"),
                     out)
        return(out)
    })
    editcall_out <- do.call("rbind", editcall_out)
    editcall_fn <- file.path(editcall_dir, "editcall_summary.csv")
    write.csv(editcall_out, editcall_fn, row.names = FALSE)
    return(editcall_out)
}

################################################################################
# Optimized version of doAlign function
################################################################################

#' Optimized alignment function
#' @param blast_path Path to BLAST executables
#' @param basecall_fn Basecalling results
#' @param demult_out Demultiplexing results
#' @param amplicon_fn Amplicon database
#' @param align_dir Alignment output directory
#' @param primer_list Primer sequences
#' @param pam_list PAM site information
#' @param genome_fn Reference genome
#' @param check_window Window size for edit detection
#' @return Alignment results
doAlign_optimized <- function(blast_path,
                             basecall_fn,
                             demult_out,
                             amplicon_fn,
                             align_dir,
                             primer_list,
                             pam_list,
                             genome_fn,
                             check_window = 10) {
    
    # Run BLAST alignment (unchanged - external tool)
    blastout_suffix <- sub("\\.fa", ".blastout", basename(amplicon_fn))
    ampl_hit_fn <- NULL
    for(i in seq_along(basecall_fn)){
        blastout_fn <- paste(sub("\\/basecall\\/", "/align/",
                                 sub("\\.fq", "", basecall_fn[i])),
                             blastout_suffix,
                             sep = "_")
        i_ampl_hit_fn <- .run_blastn(blast_path = blast_path,
                                     query_fn = amplicon_fn,
                                     db_path = basecall_fn[i],
                                     blastout_fn = blastout_fn,
                                     task = "blastn",
                                     outfmt = "long",
                                     word_size = NULL,
                                     n_core = 1)
        ampl_hit_fn <- c(ampl_hit_fn, i_ampl_hit_fn)
    }
    
    # Parse BLAST results (unchanged)
    ampl_hit <- NULL
    for(i in seq_along(ampl_hit_fn)){
        i_ampl_hit <- .parse_blastout(blastout_fn = ampl_hit_fn[i],
                                      outfmt = "long")
        ampl_hit <- rbind(ampl_hit, i_ampl_hit)
    }
    
    # Process hits (unchanged)
    n_id_hit <- table(ampl_hit$sseqid)
    single_hit <- names(n_id_hit)[n_id_hit == 1]
    multi_hit <- ampl_hit[!ampl_hit$sseqid %in% single_hit, ]
    single_hit <- ampl_hit[ampl_hit$sseqid %in% single_hit, ]
    multi_hit_check <- tapply(seq_along(multi_hit$sseqid), multi_hit$sseqid, function(i){
        q_start <- multi_hit$qstat[i]
        q_end <- multi_hit$qend[i]
        inside <- q_start >= q_start[1] & q_end <= tail(q_end, 1)
        return(all(inside))
    })
    multi_hit_check <- unlist(multi_hit_check)
    if(any(!multi_hit_check)){
        stop("Multiple hit reads")
    }
    
    ampl_hit <- ampl_hit[!duplicated(ampl_hit$sseqid), ]
    
    # Load sequences and process (unchanged)
    amplicon_seq <- readDNAStringSet(amplicon_fn)
    genome <- readDNAStringSet(genome_fn)
    edit_site <- read.csv(pam_list, header = FALSE)
    edit_site_gr <- GRanges(seqnames = paste0("chr", sprintf("%02d", edit_site$V2)),
                            ranges = IRanges(start = edit_site$V3 - check_window,
                                             end = edit_site$V3 + check_window),
                            target = edit_site$V1)
    genome_seq <- genome[edit_site_gr]
    names(genome_seq) <- edit_site_gr$target
    genome_seq <- genome_seq[order(names(genome_seq))]
    amplicon_seq <- amplicon_seq[order(names(amplicon_seq))]
    
    # Use optimized alignment processing
    aln <- Map(f = function(x, y){
        pairwiseAlignment(x, y)
    }, genome_seq, amplicon_seq)
    
    aln <- lapply(aln, aligned)
    aln <- lapply(aln, as.character)
    
    # Use optimized pattern matching
    aln_gaps <- lapply(aln, function(x) {
        positions <- fast_gregexpr_optimized(x, "-")
        positions[!is.na(positions)]
    })
    
    aln_diff <- lapply(aln_gaps, diff)
    aln_diff <- lapply(aln_diff, ">", 1)
    aln_diff <- lapply(aln_diff, which)
    aln_pos <- Map(f = function(x, y){
        c(y[x], y[x + 1])
    }, aln_diff, aln_gaps)
    aln_pos <- do.call("rbind", aln_pos)
    aln_pos <- as.data.frame(aln_pos)
    names(aln_pos) <- c("start", "end")
    aln_pos$start <- aln_pos$start + 1
    aln_pos$end <- aln_pos$end - 1
    
    # Process primers (unchanged)
    primers <- read.csv(primer_list, header = FALSE)
    f_primers <- DNAStringSet(primers$V2[grep("_F$", primers$V1)])
    names(f_primers) <- primers$V1[grep("_F$", primers$V1)]
    r_primers <- DNAStringSet(primers$V2[grep("_R$", primers$V1)])
    names(r_primers) <- primers$V1[grep("_R$", primers$V1)]
    f_primers <- f_primers[order(names(f_primers))]
    r_primers <- r_primers[order(names(r_primers))]
    aln_pos$f_primer_len <- width(f_primers)
    aln_pos$r_primer_len <- width(r_primers)
    aln_pos$dist_to_end <- width(amplicon_seq) - aln_pos$end
    aln_pos$target_len <- aln_pos$end - aln_pos$start
    
    # Process alignments with optimized string operations
    target <- names(amplicon_seq)
    aln_pos$target <- target
    align_out <- NULL
    
    for(i in seq_along(target)){
        target_gene <- target[i]
        i_hit <- subset(ampl_hit, subset = ampl_hit$qseqid == target_gene)
        i_aln_pos <- subset(aln_pos, subset = target == target_gene)
        s_cover_edit_site <- i_hit$qstart <= i_aln_pos$start
        e_cover_edit_site <- i_hit$qend >= i_aln_pos$start
        cover_edit_site <- s_cover_edit_site & e_cover_edit_site
        i_hit <- i_hit[cover_edit_site, ]
        
        # Use optimized alignment statistics
        if (nrow(i_hit) > 0) {
            aln_stats <- fast_alignment_stats_optimized(
                i_hit$qseq, i_hit$sseq, 
                rep(i_aln_pos$start, nrow(i_hit)), 
                rep(i_aln_pos$end, nrow(i_hit))
            )
            
            i_out <- data.frame(target_gene = i_hit$qseqid,
                                read_name = i_hit$sseqid,
                                read_seq = aln_stats$subject,
                                ref_seq = aln_stats$query)
            align_out <- rbind(align_out, i_out)
        }
    }
    
    # Check intact sequences (unchanged)
    intact_seq <- as.character(genome_seq)
    align_out$intact <- align_out$read_seq %in% intact_seq
    
    # Add full amplicon information (unchanged)
    full_amplicon <- data.frame(read_name = ampl_hit$sseqid,
                                full_read_seq = ampl_hit$sseq,
                                aln_start = ampl_hit$sstart,
                                aln_end = ampl_hit$send)
    full_amplicon <- unique(full_amplicon)
    align_out <- left_join(align_out, full_amplicon, "read_name")
    
    # Save results
    align_fn <- file.path(align_dir, "alignment_list.csv")
    write.csv(align_out, align_fn, row.names = FALSE)
    attributes(align_out) <- c(attributes(align_out), list(intact_seq = intact_seq))
    return(align_out)
}
