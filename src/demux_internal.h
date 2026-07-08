#ifndef MIAOSEQ_DEMUX_INTERNAL_H
#define MIAOSEQ_DEMUX_INTERNAL_H

#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "edlib.h"

namespace miaoseq {

inline char complement_base(char b) {
    switch (b) {
        case 'A': case 'a': return 'T';
        case 'T': case 't': return 'A';
        case 'G': case 'g': return 'C';
        case 'C': case 'c': return 'G';
        case 'N': case 'n': return 'N';
        default: return 'N';
    }
}

inline std::string reverse_complement(const std::string& s) {
    std::string out;
    out.resize(s.size());
    for (size_t i = 0; i < s.size(); ++i) {
        out[s.size() - 1 - i] = complement_base(s[i]);
    }
    return out;
}

inline std::string to_upper_acgt(const std::string& s) {
    std::string out = s;
    for (char& c : out) {
        if (c >= 'a' && c <= 'z') c = static_cast<char>(c - 'a' + 'A');
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T') c = 'N';
    }
    return out;
}

struct AnchorHit {
    bool found = false;
    bool flipped = false;
    int edit = -1;
    int start = -1; // 0-based inclusive start in the searched sequence
    int end = -1;   // 0-based inclusive end
    std::string seq; // sequence that was searched (native or RC window)
};

struct BarcodeHit {
    bool found = false;
    int parent_id = -1;
    int edit = -1;
    std::string barcode;
};

struct MutantEntry {
    int parent_id = -1;
    int edit = 0;
};

struct MutantDict {
    std::vector<std::string> parents;
    std::vector<std::string> parent_names;
    std::unordered_map<std::string, MutantEntry> mutant_map;
    std::vector<std::pair<std::string, std::string>> conflicts;
};

std::vector<std::string> mutate_barcode_once(const std::string& seq);

MutantDict build_mutant_dict(const std::vector<std::string>& barcodes,
                             const std::vector<std::string>& names,
                             int max_edit);

AnchorHit find_suffix_anchor(const std::string& window,
                             const std::string& suffix,
                             int max_edit,
                             bool allow_revcomp,
                             long long* n_edlib);

// All HW locations (native + optional RC), sorted by edit then start.
std::vector<AnchorHit> find_suffix_anchors(const std::string& window,
                                           const std::string& suffix,
                                           int max_edit,
                                           bool allow_revcomp,
                                           long long* n_edlib);

std::string extract_barcode(const std::string& seq, int suffix_start0, int barcode_len,
                            int offset = 0);

BarcodeHit lookup_barcode(const std::string& barcode, const MutantDict& dict);

int optional_prefix_edit(const std::string& seq,
                         const std::string& prefix,
                         int barcode_start0,
                         int max_prefix_edit,
                         long long* n_edlib);

} // namespace miaoseq

#endif
