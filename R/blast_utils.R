################################################################################
# BLAST utilities (internal)
################################################################################

#' @keywords internal
.makeblastdb <- function(blast_path, db_path){
    blast_args <- paste(paste("-in", db_path),
                        paste("-dbtype nucl"))
    check <- try({
        system2(command = file.path(blast_path, "makeblastdb"),
                args = blast_args,
                stdout = FALSE)
    }, silent = TRUE)
    if(inherits(check, "try-error")){
        stop("Error in makeblastdb of BLAST.\n", check)
    }
}

#' @keywords internal
.run_blastn <- function(blast_path,
                        query_fn,
                        db_path,
                        blastout_fn,
                        task = "blastn-short",
                        outfmt = "short",
                        word_size = 4,
                        n_core,
                        not_run = FALSE){
    if(outfmt == "short"){
        outfmt <- "'6 qseqid sseqid qstart qend sstart send length pident'"

    } else {
        outfmt <- "'6 qseqid sseqid qstart qend qlen sstart send slen mismatch gapopen sseq qseq'"
    }
    if(!not_run){
        blast_args <- paste(paste("-query", query_fn),
                            paste("-db", db_path),
                            paste("-task", task),
                            "-max_target_seqs 1000000",
                            "-max_hsps 1",
                            paste("-outfmt", outfmt),
                            paste("-num_threads", n_core),
                            paste("-out", blastout_fn))

        if(!is.null(word_size)){
            blast_args <- paste(blast_args,
                                paste("-word_size", word_size))
        }

        system2(command = file.path(blast_path, "blastn"),
                args = blast_args,
                stdout = FALSE)
    }
    return(blastout_fn)
}

#' @keywords internal
.parse_blastout <- function(blastout_fn, outfmt){
    if(outfmt == "short"){
        out_names <- c("qseqid", "sseqid", "qstart", "qend",
                       "sstart", "send", "length", "pident")

    } else {
        out_names <- c("qseqid", "sseqid", "qstart", "qend",
                       "qlen", "sstart", "send", "slen", "mismatch",
                       "gapopen", "sseq", "qseq")
    }
    out <- read.table(blastout_fn, sep = "\t")
    names(out) <- out_names
    out <- subset(out, subset = !is.na(qseqid) & !is.na(sseqid))
    return(out)
}
