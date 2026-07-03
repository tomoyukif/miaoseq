################################################################################
# Apptainer helpers for refless pipeline
################################################################################

#' Read FASTA or FASTQ sequences (format from file extension)
#' @keywords internal
.refless_read_seqs <- function(path){
    fmt <- if(grepl("\\.(fq|fastq)(\\.gz)?$", path, ignore.case = TRUE)){
        "fastq"
    } else {
        "fasta"
    }
    Biostrings::readDNAStringSet(path, format = fmt)
}

#' Resolve path to miaoseq-refless Apptainer image
#'
#' @param sif_path Optional explicit path to .sif file.
#' @return Character path to SIF, or NULL if not found.
#' @keywords internal
.refless_sif_path <- function(sif_path = NULL){
    if(!is.null(sif_path) && file.exists(sif_path)){
        return(normalizePath(sif_path, mustWork = TRUE))
    }
    candidates <- c(
        Sys.getenv("MIAOSEQ_REFLESS_SIF"),
        file.path("cursor_dev", "apptainer", "images", "miaoseq-refless.sif"),
        file.path("cursor_dev", "apptainer", "images", "miaoseq-refless.sandbox"),
        file.path(system.file(package = "miaoseq"), "..", "..", "cursor_dev",
                            "apptainer", "images", "miaoseq-refless.sif"),
        file.path(system.file(package = "miaoseq"), "..", "..", "cursor_dev",
                            "apptainer", "images", "miaoseq-refless.sandbox")
    )
    for(c in candidates){
        if(nzchar(c) && file.exists(c)){
            return(normalizePath(c, mustWork = TRUE))
        }
    }
    NULL
}

#' Path to bundled UMAP pre-cluster script
#' @keywords internal
.refless_umap_script <- function(){
    candidates <- c(
        Sys.getenv("MIAOSEQ_UMAP_SCRIPT"),
        file.path("cursor_dev", "apptainer", "scripts", "umap_kmer_cluster.py"),
        system.file("scripts", "umap_kmer_cluster.py", package = "miaoseq")
    )
    for(c in candidates){
        if(nzchar(c) && file.exists(c)){
            return(normalizePath(c, mustWork = TRUE))
        }
    }
    stop("umap_kmer_cluster.py not found. Rebuild cursor_dev/apptainer or set MIAOSEQ_UMAP_SCRIPT.")
}

#' Known host locations for external tools
#' @keywords internal
.refless_tool_candidates <- function(tool){
    env_key <- paste0("MIAOSEQ_", toupper(tool), "_PATH")
    env_path <- Sys.getenv(env_key)
    defaults <- list(
        mmseqs = c(
            "/home/ftom/07_tools/MMseqs2/build/bin/mmseqs",
            "/opt/conda/bin/mmseqs"
        ),
        minimap2 = c(
            "/home/ftom/07_tools/minimap2-2.28_x64-linux/minimap2",
            "/opt/conda/bin/minimap2"
        ),
        vsearch = c("/opt/conda/bin/vsearch"),
        python3 = c("/opt/conda/bin/python3", "/usr/bin/python3")
    )
    c(if(nzchar(env_path)) env_path, defaults[[tool]] %||% character())
}

#' Resolve executable for an external tool
#' @keywords internal
.refless_tool_bin <- function(tool){
    host <- Sys.which(tool)
    if(nzchar(host)){
        return(host)
    }
    for(path in .refless_tool_candidates(tool)){
        if(nzchar(path) && file.exists(path) && file.access(path, 1) == 0){
            return(normalizePath(path, mustWork = TRUE))
        }
    }
    NULL
}

`%||%` <- function(x, y) if(length(x)) x else y

#' Detect apptainer or singularity runner
#' @keywords internal
.refless_container_runner <- function(){
    if(nzchar(Sys.which("apptainer"))){
        return("apptainer")
    }
    if(nzchar(Sys.which("singularity"))){
        return("singularity")
    }
    NULL
}

#' Run a command inside the refless Apptainer image
#'
#' @param args Character vector of command arguments after the image.
#' @param bind Character vector of host:container bind paths.
#' @param sif_path Optional SIF path.
#' @return Exit status from system2.
#' @keywords internal
.refless_apptainer_exec <- function(args,
                                    bind = character(),
                                    sif_path = NULL,
                                    pwd = NULL,
                                    stdout = "",
                                    stderr = ""){
    runner <- .refless_container_runner()
    if(is.null(runner)){
        stop("apptainer or singularity not found in PATH")
    }
    sif <- .refless_sif_path(sif_path)
    if(is.null(sif)){
        stop("miaoseq-refless.sif not found. Build with: bash cursor_dev/apptainer/build.sh")
    }
    bind_args <- vapply(bind, function(b) paste0("--bind=", b), character(1))
    pwd_args <- if(!is.null(pwd) && nzchar(pwd)) c("--pwd", pwd) else character()
    cmd_args <- c("exec", bind_args, pwd_args, sif, args)
    system2(runner, cmd_args, stdout = stdout, stderr = stderr)
}

#' Run shell command via Apptainer or host
#' @keywords internal
.refless_run_shell <- function(cmd, work_dir = getwd(), sif_path = NULL){
    work_dir <- normalizePath(work_dir, mustWork = TRUE)
    sif <- .refless_sif_path(sif_path)
    if(!is.null(sif)){
        runner <- .refless_container_runner()
        full <- paste(
            runner, "exec",
            paste0("--bind=", work_dir, ":", work_dir),
            "--pwd", work_dir,
            sif,
            "bash", "-lc", shQuote(cmd)
        )
        return(system(full))
    }

    old_path <- Sys.getenv("PATH")
    on.exit(Sys.setenv(PATH = old_path), add = TRUE)
    extra <- unique(na.omit(dirname(c(
        .refless_tool_bin("minimap2"),
        .refless_tool_bin("mmseqs"),
        .refless_tool_bin("vsearch")
    ))))
    if(length(extra)){
        Sys.setenv(PATH = paste(c(extra, strsplit(old_path, ":", fixed = TRUE)[[1]]), collapse = ":"))
    }
    system(cmd, intern = FALSE)
}

#' Run external tool via Apptainer or host PATH
#'
#' @param tool Tool name (e.g. mmseqs, minimap2).
#' @param args Argument string.
#' @param work_dir Working directory to bind.
#' @param sif_path Optional SIF path.
#' @param prefer_host If TRUE, use host binary when available.
#' @keywords internal
.refless_run_tool <- function(tool,
                              args,
                              work_dir = getwd(),
                              sif_path = NULL,
                              prefer_host = TRUE){
    work_dir <- normalizePath(work_dir, mustWork = TRUE)
    run_in_dir <- function(status_fun){
        owd <- getwd()
        on.exit(setwd(owd), add = TRUE)
        setwd(work_dir)
        status_fun()
    }
    host_bin <- .refless_tool_bin(tool)
    if(prefer_host && !is.null(host_bin)){
        return(run_in_dir(function() system2(host_bin, args, stdout = "", stderr = "")))
    }
    sif <- .refless_sif_path(sif_path)
    if(!is.null(sif)){
        return(.refless_apptainer_exec(
            args = c(tool, strsplit(args, " ", fixed = TRUE)[[1]]),
            bind = paste0(work_dir, ":", work_dir),
            pwd = work_dir,
            sif_path = sif
        ))
    }
    if(!is.null(host_bin)){
        return(run_in_dir(function() system2(host_bin, args, stdout = "", stderr = "")))
    }
    stop("Tool ", tool, " not found. Build Apptainer image or set MIAOSEQ_", toupper(tool), "_PATH.")
}

#' Run UMAP pre-cluster script on host Python or in container
#' @keywords internal
.refless_run_umap_script <- function(fasta_fn,
                                     out_tsv,
                                     min_samples = 50,
                                     sif_path = NULL){
    script <- .refless_umap_script()
    wd <- normalizePath(dirname(fasta_fn), mustWork = TRUE)
    fa <- basename(fasta_fn)
    out <- basename(out_tsv)

    py <- .refless_tool_bin("python3")
    if(!is.null(py)){
        args <- c(
            script,
            file.path(wd, fa),
            file.path(wd, out),
            "--min-samples", min_samples
        )
        status <- system2(py, args, stdout = "", stderr = "")
        if(is.null(status) || status == 0){
            return(status)
        }
    }

    sif <- .refless_sif_path(sif_path)
    if(!is.null(sif)){
        bind <- unique(c(
            paste0(wd, ":", wd),
            paste0(dirname(script), ":", dirname(script))
        ))
        return(.refless_apptainer_exec(
            args = c(
                "python3", script,
                file.path(wd, fa),
                file.path(wd, out),
                "--min-samples", as.character(min_samples)
            ),
            bind = bind,
            pwd = wd,
            sif_path = sif
        ))
    }

    stop("UMAP pre-cluster requires python3 with umap-learn/hdbscan, or miaoseq-refless.sif")
}

#' Read primer sequences from CSV (col1=id, col2=seq)
#' @keywords internal
.read_primers <- function(primer_list){
    primers <- read.csv(primer_list, header = FALSE, stringsAsFactors = FALSE)
    f <- primers$V2[grep("_F$", primers$V1)]
    r <- primers$V2[grep("_R$", primers$V1)]
    if(length(f) == 0 || length(r) == 0){
        stop("primer_list must contain IDs ending with _F and _R")
    }
    list(f = f[1], r = r[1])
}

#' Parse index list CSV (5 columns)
#' @keywords internal
.read_index_list <- function(index_list){
    df <- read.csv(index_list, header = FALSE, stringsAsFactors = FALSE)
    if(ncol(df) < 5){
        stop("index_list requires 5 columns: id, f_index_id, f_index, r_index_id, r_index")
    }
    names(df) <- c("sample_id", "f_index_id", "f_index", "r_index_id", "r_index")
    df
}

#' Reverse complement DNA string
#' @keywords internal
.revcomp <- function(x){
    as.character(Biostrings::reverse(Biostrings::complement(Biostrings::DNAString(x))))
}
