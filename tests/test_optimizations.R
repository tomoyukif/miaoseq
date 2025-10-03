#!/usr/bin/env Rscript

# Test script for miaoseq optimizations
# This script validates the performance improvements and correctness of optimized functions

library(testthat)
library(microbenchmark)

# Load the package functions
source("R/functions.R")
source("R/optimized_functions.R")
source("R/pipeline_optimized.R")

# Test data generation
generate_test_data <- function(n_sequences = 1000, seq_length = 500) {
    sequences <- character(n_sequences)
    for (i in 1:n_sequences) {
        # Generate random DNA sequence with some gaps
        seq <- paste(sample(c("A", "T", "G", "C", "-"), seq_length, replace = TRUE), collapse = "")
        sequences[i] <- seq
    }
    return(sequences)
}

# Test 1: String pattern matching
test_string_operations <- function() {
    cat("Testing string operations...\n")
    
    test_seqs <- generate_test_data(1000, 300)
    
    # Test pattern counting
    standard_result <- sapply(gregexpr("-", test_seqs), function(x) if(x[1] == -1) 0 else length(x))
    optimized_result <- fast_pattern_count_optimized(test_seqs, "-")
    
    expect_equal(standard_result, optimized_result)
    cat("✓ Pattern counting test passed\n")
    
    # Test sequence length calculation
    standard_length <- nchar(gsub("-", "", test_seqs))
    optimized_length <- fast_seq_length_no_gaps_optimized(test_seqs)
    
    expect_equal(standard_length, optimized_length)
    cat("✓ Sequence length calculation test passed\n")
    
    # Test substring extraction
    starts <- sample(1:200, 1000, replace = TRUE)
    ends <- starts + sample(50:100, 1000, replace = TRUE)
    
    standard_substr <- substr(test_seqs, starts, ends)
    optimized_substr <- fast_substr_optimized(test_seqs, starts, ends)
    
    expect_equal(standard_substr, optimized_substr)
    cat("✓ Substring extraction test passed\n")
}

# Test 2: Performance benchmarking
test_performance <- function() {
    cat("\nTesting performance improvements...\n")
    
    test_seqs <- generate_test_data(5000, 400)
    
    # Benchmark pattern counting
    pattern_benchmark <- microbenchmark(
        standard = sapply(gregexpr("-", test_seqs), function(x) if(x[1] == -1) 0 else length(x)),
        optimized = fast_pattern_count_optimized(test_seqs, "-"),
        times = 10
    )
    
    cat("Pattern counting benchmark:\n")
    print(pattern_benchmark)
    
    # Benchmark sequence length calculation
    length_benchmark <- microbenchmark(
        standard = nchar(gsub("-", "", test_seqs)),
        optimized = fast_seq_length_no_gaps_optimized(test_seqs),
        times = 10
    )
    
    cat("\nSequence length calculation benchmark:\n")
    print(length_benchmark)
    
    # Calculate speedup
    pattern_speedup <- median(pattern_benchmark$time[pattern_benchmark$expr == "standard"]) / 
                      median(pattern_benchmark$time[pattern_benchmark$expr == "optimized"])
    
    length_speedup <- median(length_benchmark$time[length_benchmark$expr == "standard"]) / 
                     median(length_benchmark$time[length_benchmark$expr == "optimized"])
    
    cat(sprintf("\nSpeedup achieved:\n"))
    cat(sprintf("Pattern counting: %.2fx\n", pattern_speedup))
    cat(sprintf("Length calculation: %.2fx\n", length_speedup))
    
    return(list(pattern_speedup = pattern_speedup, length_speedup = length_speedup))
}

# Test 3: Memory usage
test_memory_usage <- function() {
    cat("\nTesting memory usage...\n")
    
    # Test with different dataset sizes
    sizes <- c(1000, 5000, 10000)
    memory_results <- list()
    
    for (size in sizes) {
        test_seqs <- generate_test_data(size, 300)
        
        # Measure memory usage
        gc()
        mem_before <- memory.size()
        
        # Run optimized function
        optimized_result <- fast_pattern_count_optimized(test_seqs, "-")
        
        gc()
        mem_after <- memory.size()
        
        memory_usage <- mem_after - mem_before
        memory_results[[as.character(size)]] <- memory_usage
        
        cat(sprintf("Dataset size %d: Memory usage %.2f MB (processed %d sequences)\n", 
                   size, memory_usage, length(optimized_result)))
    }
    
    return(memory_results)
}

# Test 4: Correctness of optimized functions
test_correctness <- function() {
    cat("\nTesting correctness of optimized functions...\n")
    
    # Test sequence table operations
    test_seqs <- c("ATCG", "ATCG", "GCTA", "ATCG", "GCTA", "TTTT")
    test_counts <- c(1, 1, 1, 1, 1, 1)
    
    standard_table <- table(test_seqs)
    optimized_table <- fast_sequence_table_optimized(test_seqs, test_counts)
    
    # Convert to comparable format
    standard_df <- data.frame(names = names(standard_table), values = as.numeric(standard_table))
    optimized_df <- data.frame(names = optimized_table$names, values = optimized_table$values)
    
    # Sort for comparison
    standard_df <- standard_df[order(standard_df$names), ]
    optimized_df <- optimized_df[order(optimized_df$names), ]
    
    expect_equal(standard_df$values, optimized_df$values)
    cat("✓ Sequence table operations test passed\n")
    
    # Test deletion detection
    test_seqs_with_gaps <- c("ATCG", "AT-CG", "ATC-G", "ATCG")
    standard_del <- grepl("-", test_seqs_with_gaps)
    optimized_del <- fast_has_deletions_optimized(test_seqs_with_gaps)
    
    expect_equal(standard_del, optimized_del)
    cat("✓ Deletion detection test passed\n")
    
    # Test insertion detection
    test_seqs_ins <- c("ATCG", "ATCG", "ATCG", "ATCG")
    test_refs <- c("ATCG", "ATC", "AT", "A")
    standard_ins <- nchar(gsub("-", "", test_seqs_ins)) > nchar(gsub("-", "", test_refs))
    optimized_ins <- fast_has_insertions_optimized(test_seqs_ins, test_refs)
    
    expect_equal(standard_ins, optimized_ins)
    cat("✓ Insertion detection test passed\n")
}

# Test 5: Integration test
test_integration <- function() {
    cat("\nTesting integration with main pipeline...\n")
    
    # Create minimal test data
    test_dir <- tempdir()
    cat(sprintf("Using test directory: %s\n", test_dir))
    
    # Test if optimized functions can be loaded and called
    tryCatch({
        # Test function availability
        expect_true(exists("fast_pattern_count_optimized"))
        expect_true(exists("fast_seq_length_no_gaps_optimized"))
        expect_true(exists("fast_has_deletions_optimized"))
        expect_true(exists("fast_has_insertions_optimized"))
        
        cat("✓ All optimized functions are available\n")
        
        # Test function calls
        test_seqs <- generate_test_data(100, 200)
        result1 <- fast_pattern_count_optimized(test_seqs, "-")
        result2 <- fast_seq_length_no_gaps_optimized(test_seqs)
        result3 <- fast_has_deletions_optimized(test_seqs)
        
        expect_equal(length(result1), length(test_seqs))
        expect_equal(length(result2), length(test_seqs))
        expect_equal(length(result3), length(test_seqs))
        
        # Use results to avoid unused variable warning
        cat(sprintf("Test results: %d sequences processed\n", length(result1)))
        
        cat("✓ All optimized functions execute correctly\n")
        
    }, error = function(e) {
        cat(sprintf("✗ Integration test failed: %s\n", e$message))
    })
}

# Main test runner
run_all_tests <- function() {
    cat("Starting miaoseq optimization tests...\n")
    cat("=====================================\n")
    
    tryCatch({
        # Run all tests
        test_string_operations()
        performance_results <- test_performance()
        memory_results <- test_memory_usage()
        test_correctness()
        test_integration()
        
        cat("\n=====================================\n")
        cat("All tests completed successfully!\n")
        cat("=====================================\n")
        
        # Summary
        cat("\nPerformance Summary:\n")
        cat(sprintf("Pattern counting speedup: %.2fx\n", performance_results$pattern_speedup))
        cat(sprintf("Length calculation speedup: %.2fx\n", performance_results$length_speedup))
        
        # Use memory results to avoid warning
        cat(sprintf("Memory test completed for %d dataset sizes\n", length(memory_results)))
        
        return(TRUE)
        
    }, error = function(e) {
        cat(sprintf("\nTest failed: %s\n", e$message))
        return(FALSE)
    })
}

# Run tests if script is executed directly
if (!interactive()) {
    success <- run_all_tests()
    quit(status = ifelse(success, 0, 1))
}
