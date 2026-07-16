.read_plate_layout <- function(sample_list) {
    if (is.null(sample_list)) {
        return(NULL)
    }
    if (!file.exists(sample_list)) {
        stop("sample_list not found: ", sample_list)
    }
    sl <- utils::read.csv(sample_list, header = FALSE, stringsAsFactors = FALSE)
    if (ncol(sl) < 5L) {
        stop(
            "sample_list needs 5 columns: ",
            "index_pair_id, sample_name, plate_id, row_id (A-H), col_id (1-12)"
        )
    }
    if (.looks_like_index_list(sl)) {
        stop(
            "sample_list looks like index_list (barcode sequences). ",
            "Provide a plate layout CSV with columns: ",
            "index_pair_id, sample_name, plate, row (A-H), col (1-12)."
        )
    }
    out <- data.frame(
        index_pair_id = as.character(sl[[1L]]),
        sample_name = as.character(sl[[2L]]),
        plate_id = as.character(sl[[3L]]),
        row_id = toupper(trimws(as.character(sl[[4L]]))),
        col_id = as.character(suppressWarnings(as.integer(sl[[5L]]))),
        stringsAsFactors = FALSE
    )
    bad_row <- !.is_plate_row(out$row_id)
    bad_col <- !.is_plate_col(out$col_id)
    if (any(bad_row)) {
        stop(
            "Invalid row_id in sample_list (expected A-H): ",
            paste(unique(out$row_id[bad_row]), collapse = ", ")
        )
    }
    if (any(bad_col)) {
        stop(
            "Invalid col_id in sample_list (expected 1-12): ",
            paste(unique(out$col_id[bad_col]), collapse = ", ")
        )
    }
    out
}

.looks_like_index_list <- function(sl) {
    if (ncol(sl) < 5L) {
        return(FALSE)
    }
    v3 <- as.character(sl[[3L]])
    v4 <- as.character(sl[[4L]])
    v5 <- as.character(sl[[5L]])
    any(nchar(v3) > 20L) &&
        any(grepl("^miao_I[57]_index", v4, ignore.case = TRUE)) &&
        any(nchar(v5) > 20L)
}

.is_plate_row <- function(x) {
    toupper(trimws(as.character(x))) %in% LETTERS[1:8]
}

.is_plate_col <- function(x) {
    xc <- suppressWarnings(as.integer(x))
    !is.na(xc) & xc >= 1L & xc <= 12L
}

.has_valid_plate_layout <- function(df) {
    if (!all(c("row_id", "col_id") %in% names(df))) {
        return(FALSE)
    }
    all(.is_plate_row(df$row_id), na.rm = TRUE) &&
        all(.is_plate_col(df$col_id), na.rm = TRUE)
}

.join_plate_layout <- function(edit_result, layout) {
    hit <- match(edit_result$index_pair_id, layout$index_pair_id)
    edit_result$name <- layout$sample_name[hit]
    edit_result$plate <- layout$plate_id[hit]
    edit_result$row <- layout$row_id[hit]
    edit_result$col <- as.character(layout$col_id[hit])
    missing <- is.na(hit)
    if (any(missing)) {
        warning(
            sum(missing), " sample(s) in editcall_summary lack plate layout: ",
            paste(head(edit_result$index_pair_id[missing], 5L), collapse = ", "),
            if (sum(missing) > 5L) " ..." else "",
            call. = FALSE
        )
    }
    edit_result[!missing, , drop = FALSE]
}
