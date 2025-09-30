# miaoseq

**miaoseq** is an R package for analyzing CRISPR-Cas9 editing outcomes from Oxford Nanopore MinION sequencing data. The package provides a comprehensive pipeline for basecalling, demultiplexing, read alignment, and edit-calling analysis.

## Overview

The miaoseq pipeline processes raw MinION sequencing data (pod5 format) through the following steps:

1. **Basecalling**: Converts raw electrical signals to DNA sequences using Dorado
2. **Demultiplexing**: Assigns reads to samples based on index sequences using BLAST
3. **Read Alignment**: Aligns reads to amplicon sequences using BLAST
4. **Edit-calling**: Identifies and characterizes CRISPR-Cas9 editing outcomes
5. **Evaluation**: Generates comprehensive statistical summaries

## Prerequisites

### External Software Requirements

miaoseq requires three external tools that must be installed separately from R:

#### 1. Dorado (Oxford Nanopore Technologies)
Dorado is used for basecalling raw MinION data.  
For installation instructions, please visit:  
https://github.com/nanoporetech/dorado

#### 2. BLAST (NCBI)
BLAST is used for sequence alignment and demultiplexing.  
For installation instructions, please visit:  
https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download

#### 3. Samtools
Samtools is used for BAM file processing.  
For installation instructions, please visit:  
https://github.com/samtools/samtools

### R Package Dependencies

miaoseq requires several R packages that will be installed automatically:

- **Biostrings**: For DNA sequence manipulation
- **dplyr**: For data manipulation
- **GenomicRanges**: For genomic interval operations
- **IRanges**: For interval operations
- **BiocGenerics**: For generic functions
- **pwalign**: For pairwise sequence alignment


## Installation

### Install miaoseq

```r
# Install miaoseq from GitHub using devtools:
devtools::install_github("tomoyukif/miaoseq")
```

### Load the package

```r
library(miaoseq)
```

## Usage Guide

### Step 1: Set up your analysis

First, define the paths and parameters for your analysis:

```r
# Load the package
source("R/functions.R")

# Define working directories
working_dir <- "/path/to/your/working/directory"
out_dir <- file.path(working_dir, "output_directory_name")
in_dir <- "/path/to/pod5/directory"  # MinION outputs raw sequence data files in a pod5 directory

# Reference files
genome_fn <- '/reference/genome/sequence.fa'  # Reference genome sequence in FASTA format
pam_list <- "inst/extdata/agr8_pam_list.csv"  # PAM site information
index_list <- "inst/extdata/index_list.csv"   # Index sequences for demultiplexing
primer_list <- "inst/extdata/amplicon_primers.csv"  # Primer sequences

# External tool paths
dorado_path <- "/path/to/dorado"              # Path to dorado executable
samtools_path <- "/path/to/samtools"          # Path to samtools executable
blast_path <- "/path/to/blast/bin"            # Path to a directory containing blast executables (makeblastdb and blastn)

# Analysis parameters
n_core <- 30  # Number of CPU cores to use
```

### Step 2: Prepare amplicon database

The `prepAmpliconDB()` function extracts amplicon sequences from your reference genome based on primer sequences:

```r
# Prepare amplicon database
amplicon_fn <- prepAmpliconDB(blast_path = blast_path,
                              primer_list = primer_list,
                              genome_fn = genome_fn,
                              out_dir = out_dir,
                              n_core = n_core)

# The function returns the path to the generated amplicon FASTA file
amplicon_fn <- file.path(out_dir, "ref/amplicon.fa")
```

**What this step does:**
- Reads primer sequences from the CSV file
- Uses BLAST to locate primers in the reference genome
- Extracts amplicon sequences between primer pairs
- Creates an indexed amplicon database for read alignment

### Step 3: Run the main analysis pipeline

The `miaoEditcall()` function orchestrates the entire analysis pipeline:

```r
# Call edits using the main pipeline
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
                             size_sel = c(300, 450),    # Valid range of read length for edit-calling (bp)
                             check_window = 10,         # Window size around expected cut site (bp)
                             n_core = n_core,
                             resume = FALSE)             # Set to TRUE to resume from previous run
```

**Pipeline steps performed by `miaoEditcall()`:**

1. **Basecalling**: Converts pod5 files to FASTQ using Dorado
2. **Quality filtering**: Removes low-quality reads and applies size selection
3. **Demultiplexing**: Assigns reads to samples based on index sequences
4. **Read alignment**: Aligns reads to amplicon sequences
5. **Edit-calling**: Identifies CRISPR-Cas9 editing outcomes

### Step 4: Generate evaluation report

The `evalMiao()` function creates comprehensive statistical summaries:

```r
# Generate statistical summary
evalMiao(out_dir = out_dir,
         output_reads = FALSE)  # Set to TRUE to output sequences of undemultiplexed reads for debugging
```

**Output files created:**
- `miao_summary/read_stats.tsv`: Overall read statistics
- `miao_summary/indexed_reads_per_gene.tsv`: Demultiplexing statistics
- `miao_summary/aligned_reads_per_gene.tsv`: Alignment statistics
- `editcall/editcall_summary.csv`: Final edit-calling results

## Parameter Configuration

### Key Parameters Explained

- **`size_sel`**: `c(min_length, max_length)` - Range of read lengths to retain (in bp)
  - Adjust based on your amplicon size
  - Example: `c(300, 450)` for 300-450 bp amplicons

- **`check_window`**: Window size around expected cut site (in bp)
  - Defines the region where edits are searched
  - Example: `10` searches 10 bp upstream and downstream of cut site

- **`n_core`**: Number of CPU cores for parallel processing
  - Higher values speed up analysis but require more memory
  - Recommended: 20-30 cores for typical datasets

- **`resume`**: Whether to resume from previous run
  - Set to `TRUE` if analysis was interrupted
  - Skips completed steps and continues from where it left off

### Input File Formats

#### Primer List (`primer_list`)
CSV file with two columns:
- Column 1: Primer ID (must end with "_F" for forward, "_R" for reverse)
- Column 2: Primer sequence

#### Index List (`index_list`)
CSV file with five columns:
- Column 1: Index pair ID
- Column 2: Forward index ID
- Column 3: Forward index sequence
- Column 4: Reverse index ID
- Column 5: Reverse index sequence

#### PAM List (`pam_list`)
CSV file with three columns:
- Column 1: Target gene name
- Column 2: Chromosome number
- Column 3: PAM position

## Output Files

The analysis creates several output directories:

- **`basecall/`**: Basecalling results and quality-filtered reads
- **`demultiplex/`**: Demultiplexing results and sample assignments
- **`align/`**: Read alignment results
- **`editcall/`**: Edit-calling results and summaries
- **`ref/`**: Reference files and amplicon database
- **`miao_summary/`**: Statistical summaries and reports

## Troubleshooting

### Common Issues

1. **External tool not found**
   - Ensure all external tools (dorado, blast, samtools) are installed and in your PATH
   - Check that the paths specified in the script are correct

2. **Memory issues**
   - Reduce `n_core` parameter if running out of memory
   - Ensure sufficient disk space for intermediate files

3. **BLAST errors**
   - Verify that BLAST executables (`makeblastdb`, `blastn`) are accessible
   - Check that input sequences are in proper FASTA format

4. **Demultiplexing failures**
   - Verify index sequences in your index list file
   - Check that index sequences are present in your reads

### Getting Help

For issues or questions:
1. Check the error messages carefully
2. Verify all input files are properly formatted
3. Ensure external tools are correctly installed
4. Check that you have sufficient computational resources

## Example Workflow

Here's a complete example workflow:

```r
# Load package
library("miaoseq")

# Set up paths
working_dir <- "/home/user/analysis"
out_dir <- file.path(working_dir, "miaoseq_results")
# The 'in_dir' variable should point to the directory containing your MinION raw data in pod5 format.
# Typically, after a MinION run, pod5 files are located in a subdirectory named 'pod5'.
in_dir <- "/data/minion_run/pod5"
genome_fn <- "/reference/genome.fa"

# External tool paths
dorado_path <- "/usr/local/bin/dorado"
samtools_path <- "/usr/local/bin/samtools"
blast_path <- "/usr/local/bin"

# Analysis parameters
n_core <- 30
size_sel <- c(300, 500)
check_window <- 10

# Run analysis
amplicon_fn <- prepAmpliconDB(blast_path = blast_path,
                              primer_list = "primers.csv",
                              genome_fn = genome_fn,
                              out_dir = out_dir,
                              n_core = n_core)

editcall_out <- miaoEditcall(in_dir = in_dir,
                             out_dir = out_dir,
                             dorado_path = dorado_path,
                             samtools_path = samtools_path,
                             blast_path = blast_path,
                             primer_list = "primers.csv",
                             pam_list = "pam_sites.csv",
                             index_list = "indices.csv",
                             genome_fn = genome_fn,
                             amplicon_fn = amplicon_fn,
                             size_sel = size_sel,
                             check_window = check_window,
                             n_core = n_core,
                             resume = FALSE)

# Generate summary
evalMiao(out_dir = out_dir, output_reads = FALSE)
```

> **Typical running time and memory usage:**  
> On a standard workstation using 30 CPU cores, processing the output read data from a single Oxford Nanopore Flongle cell typically takes about **1 hour**. Actual running time may vary depending on hardware, data size, and parameter settings.  
> **Memory usage:** In a typical run, miaoseq required approximately **20–30 GB of RAM**. Please ensure your system has sufficient memory available for large datasets.
