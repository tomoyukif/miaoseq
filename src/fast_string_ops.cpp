#include <Rcpp.h>
#include <string>
#include <vector>
#include <algorithm>
#include <unordered_map>

using namespace Rcpp;

// Fast string pattern matching and extraction
// [[Rcpp::export]]
CharacterVector fast_gregexpr(const CharacterVector& text, const std::string& pattern) {
    int n = text.size();
    CharacterVector result(n);
    
    for (int i = 0; i < n; i++) {
        std::string str = as<std::string>(text[i]);
        std::vector<int> positions;
        
        size_t pos = 0;
        while ((pos = str.find(pattern, pos)) != std::string::npos) {
            positions.push_back(pos + 1); // R uses 1-based indexing
            pos += pattern.length();
        }
        
        if (positions.empty()) {
            result[i] = NA_STRING;
        } else {
            result[i] = positions[0]; // Return first match position
        }
    }
    
    return result;
}

// Fast substring extraction with multiple positions
// [[Rcpp::export]]
CharacterVector fast_substr_multiple(const CharacterVector& text, 
                                   const IntegerVector& starts, 
                                   const IntegerVector& ends) {
    int n = text.size();
    int m = starts.size();
    CharacterVector result(m);
    
    for (int i = 0; i < m; i++) {
        if (i < n) {
            std::string str = as<std::string>(text[i]);
            int start = starts[i] - 1; // Convert to 0-based
            int end = ends[i] - 1;     // Convert to 0-based
            
            if (start >= 0 && end < str.length() && start <= end) {
                result[i] = str.substr(start, end - start + 1);
            } else {
                result[i] = "";
            }
        } else {
            result[i] = "";
        }
    }
    
    return result;
}

// Fast counting of pattern occurrences
// [[Rcpp::export]]
IntegerVector fast_pattern_count(const CharacterVector& text, const std::string& pattern) {
    int n = text.size();
    IntegerVector result(n);
    
    for (int i = 0; i < n; i++) {
        std::string str = as<std::string>(text[i]);
        int count = 0;
        size_t pos = 0;
        
        while ((pos = str.find(pattern, pos)) != std::string::npos) {
            count++;
            pos += pattern.length();
        }
        
        result[i] = count;
    }
    
    return result;
}

// Fast string replacement
// [[Rcpp::export]]
CharacterVector fast_gsub(const CharacterVector& text, 
                         const std::string& pattern, 
                         const std::string& replacement) {
    int n = text.size();
    CharacterVector result(n);
    
    for (int i = 0; i < n; i++) {
        std::string str = as<std::string>(text[i]);
        size_t pos = 0;
        
        while ((pos = str.find(pattern, pos)) != std::string::npos) {
            str.replace(pos, pattern.length(), replacement);
            pos += replacement.length();
        }
        
        result[i] = str;
    }
    
    return result;
}

// Fast sequence length calculation (excluding gaps)
// [[Rcpp::export]]
IntegerVector fast_seq_length_no_gaps(const CharacterVector& sequences) {
    int n = sequences.size();
    IntegerVector result(n);
    
    for (int i = 0; i < n; i++) {
        std::string seq = as<std::string>(sequences[i]);
        int count = 0;
        for (char c : seq) {
            if (c != '-') count++;
        }
        result[i] = count;
    }
    
    return result;
}

// Fast sequence comparison for edit detection
// [[Rcpp::export]]
LogicalVector fast_has_deletions(const CharacterVector& sequences) {
    int n = sequences.size();
    LogicalVector result(n);
    
    for (int i = 0; i < n; i++) {
        std::string seq = as<std::string>(sequences[i]);
        result[i] = seq.find('-') != std::string::npos;
    }
    
    return result;
}

// Fast sequence comparison for insertions
// [[Rcpp::export]]
LogicalVector fast_has_insertions(const CharacterVector& sequences, 
                                const CharacterVector& references) {
    int n = sequences.size();
    LogicalVector result(n);
    
    for (int i = 0; i < n; i++) {
        if (i < references.size()) {
            std::string seq = as<std::string>(sequences[i]);
            std::string ref = as<std::string>(references[i]);
            
            // Count non-gap characters
            int seq_len = 0, ref_len = 0;
            for (char c : seq) if (c != '-') seq_len++;
            for (char c : ref) if (c != '-') ref_len++;
            
            result[i] = seq_len > ref_len;
        } else {
            result[i] = false;
        }
    }
    
    return result;
}

// Fast table operation for sequence counting
// [[Rcpp::export]]
List fast_sequence_table(const CharacterVector& sequences, 
                        const IntegerVector& counts) {
    std::unordered_map<std::string, int> table;
    int n = sequences.size();
    
    for (int i = 0; i < n; i++) {
        std::string seq = as<std::string>(sequences[i]);
        int count = counts[i];
        table[seq] += count;
    }
    
    CharacterVector names(table.size());
    IntegerVector values(table.size());
    int idx = 0;
    
    for (const auto& pair : table) {
        names[idx] = pair.first;
        values[idx] = pair.second;
        idx++;
    }
    
    return List::create(
        Named("names") = names,
        Named("values") = values
    );
}

// Fast vectorized operations for alignment processing
// [[Rcpp::export]]
DataFrame fast_alignment_stats(const CharacterVector& query_seqs,
                              const CharacterVector& subject_seqs,
                              const IntegerVector& starts,
                              const IntegerVector& ends) {
    int n = query_seqs.size();
    
    CharacterVector extracted_query(n);
    CharacterVector extracted_subject(n);
    IntegerVector lengths(n);
    LogicalVector has_gaps(n);
    
    for (int i = 0; i < n; i++) {
        std::string qseq = as<std::string>(query_seqs[i]);
        std::string sseq = as<std::string>(subject_seqs[i]);
        
        int start = starts[i] - 1; // Convert to 0-based
        int end = ends[i] - 1;     // Convert to 0-based
        
        if (start >= 0 && end < qseq.length() && start <= end) {
            extracted_query[i] = qseq.substr(start, end - start + 1);
            extracted_subject[i] = sseq.substr(start, end - start + 1);
            
            std::string extracted = qseq.substr(start, end - start + 1);
            lengths[i] = extracted.length();
            has_gaps[i] = extracted.find('-') != std::string::npos;
        } else {
            extracted_query[i] = "";
            extracted_subject[i] = "";
            lengths[i] = 0;
            has_gaps[i] = false;
        }
    }
    
    return DataFrame::create(
        Named("query") = extracted_query,
        Named("subject") = extracted_subject,
        Named("length") = lengths,
        Named("has_gaps") = has_gaps
    );
}
