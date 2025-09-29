source("R/functions.R")

working_dir <- "/path/to/your/working/directory"
out_dir <- file.path(working_dir, "output_directory_name")
in_dir <- "/path/to/pod5/directory" # MinION outputs raw sequence data files in a pod5 directory
genome_fn <- '/reference/genome/sequence.fa' # Reference genome sequence in fasta format
pam_list <- "inst/extdata/agr8_pam_list.csv"
index_list <- "inst/extdata/index_list.csv"
primer_list <- "inst/extdata/amplicon_primers.csv"

dorado_path <- "/path/to/dorado" # Path to dorado executable
samtools_path <- "/path/to/samtools" # Path to samtools executable
blast_path <- "/path/to/blast/bin" # Path to a directory containing blast executables (makeblastdb and blastn)

n_core <- 30 # Number of CPU cores to use

# Prepare amplicon database
amplicon_fn <- prepAmpliconDB(blast_path = blast_path,
                              primer_list = primer_list,
                              genome_fn = genome_fn,
                              out_dir = out_dir,
                              n_core = n_core)
amplicon_fn <- file.path(out_dir, "ref/amplicon.fa")

# Call edits
editcall_out <- miaoEditcall(in_dir = in_dir,
                             out_dir = out_dir,
                             dorado_path = dorado_path,
                             samtools_path = samtools_path,
                             blast_path = blast_path,
                             primer_list = primer_list,
                             pam_list = pam_list,
                             index_list = index_list,
                             genome_fn = genome_fn,
                             amplicon_fn = amplicon_fn,
                             size_sel = c(300, 450), # Set the valid range of read length for edit-calling in bp. Adjust based on your amplicon size.
                             check_window = 10, # Set the window size (in bp) around the expected cut site to search for edits. Adjust based on your experimental design.
                             n_core = n_core,
                             resume = FALSE) # Set resume = TRUE to resume from previous run

# Generate statistical summary of demultiplexing, read-alignment, and edit-calling
evalMiao(out_dir = out_dir,
         output_reads = FALSE) # Set output_reads = TRUE to output read sequences of undemultiplexed (erroneous) reads for debugging.
