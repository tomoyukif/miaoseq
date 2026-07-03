################################################################################
# Index design evaluation utilities
################################################################################

#' Pairwise edit distance between DNA strings (unit-cost)
#' @keywords internal
.edit_distance <- function(a, b){
    a <- toupper(as.character(a))
    b <- toupper(as.character(b))
    la <- nchar(a)
    lb <- nchar(b)
    if(la == 0){
        return(lb)
    }
    if(lb == 0){
        return(la)
    }
    mat <- matrix(0L, nrow = la + 1L, ncol = lb + 1L)
    mat[, 1] <- 0:la
    mat[1, ] <- 0:lb
    ca <- strsplit(a, "", fixed = TRUE)[[1]]
    cb <- strsplit(b, "", fixed = TRUE)[[1]]
    for(i in seq_len(la)){
        for(j in seq_len(lb)){
            cost <- if(ca[i] == cb[j]) 0L else 1L
            mat[i + 1L, j + 1L] <- min(
                mat[i, j + 1L] + 1L,
                mat[i + 1L, j] + 1L,
                mat[i, j] + cost
            )
        }
    }
    mat[la + 1L, lb + 1L]
}

#' Evaluate minimum edit distance among custom index sets
#'
#' Computes pairwise Levenshtein distances for forward and reverse index
#' sequences to help assess demultiplex collision risk.
#'
#' @param index_list CSV with columns: sample_id, f_index_id, f_index, r_index_id, r_index.
#' @param out_dir Optional output directory for CSV files.
#' @return List with distance matrices and summary tables.
#' @export
evaluateIndexSet <- function(index_list, out_dir = NULL){
    samples <- .read_index_list(index_list)
    n <- nrow(samples)

    f_dist <- matrix(0L, nrow = n, ncol = n,
                     dimnames = list(samples$sample_id, samples$sample_id))
    r_dist <- f_dist
    for(i in seq_len(n)){
        for(j in seq_len(n)){
            f_dist[i, j] <- .edit_distance(samples$f_index[i], samples$f_index[j])
            r_dist[i, j] <- .edit_distance(samples$r_index[i], samples$r_index[j])
        }
    }

    diag(f_dist) <- Inf
    diag(r_dist) <- Inf
    min_f <- apply(f_dist, 1, min)
    min_r <- apply(r_dist, 1, min)
    summary_tab <- data.frame(
        sample_id = samples$sample_id,
        f_index_len = nchar(samples$f_index),
        r_index_len = nchar(samples$r_index),
        min_f_edit_distance = min_f,
        min_r_edit_distance = min_r,
        min_combined_edit_distance = pmin(min_f, min_r),
        stringsAsFactors = FALSE
    )

    if(!is.null(out_dir)){
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        write.csv(f_dist, file.path(out_dir, "f_index_edit_distance.csv"))
        write.csv(r_dist, file.path(out_dir, "r_index_edit_distance.csv"))
        write.csv(summary_tab, file.path(out_dir, "index_distance_summary.csv"), row.names = FALSE)
    }

    invisible(list(
        f_index_distance = f_dist,
        r_index_distance = r_dist,
        summary = summary_tab
    ))
}
