#include "demux_internal.h"

#include <algorithm>
#include <set>

namespace miaoseq {

std::vector<std::string> mutate_barcode_once(const std::string& seq) {
    static const char bases[4] = {'A', 'C', 'G', 'T'};
    std::unordered_set<std::string> out;
    const size_t n = seq.size();

    for (size_t i = 0; i < n; ++i) {
        for (char b : bases) {
            if (b == seq[i]) continue;
            std::string m = seq;
            m[i] = b;
            out.insert(m);
        }
    }
    if (n > 1) {
        for (size_t i = 0; i < n; ++i) {
            std::string m = seq;
            m.erase(i, 1);
            out.insert(m);
        }
    }
    for (size_t i = 0; i <= n; ++i) {
        for (char b : bases) {
            std::string m = seq;
            m.insert(m.begin() + static_cast<long>(i), b);
            out.insert(m);
        }
    }
    return std::vector<std::string>(out.begin(), out.end());
}

MutantDict build_mutant_dict(const std::vector<std::string>& barcodes,
                             const std::vector<std::string>& names,
                             int max_edit) {
    MutantDict dict;
    dict.parents = barcodes;
    dict.parent_names = names;
    if (dict.parent_names.size() != dict.parents.size()) {
        dict.parent_names.resize(dict.parents.size());
    }

    // mutant -> (set of parents, min edit seen)
    struct Own {
        std::set<int> parents;
        int min_edit = 99;
    };
    std::unordered_map<std::string, Own> owners;

    auto add_owner = [&](const std::string& seq, int id, int edit) {
        Own& o = owners[seq];
        o.parents.insert(id);
        if (edit < o.min_edit) o.min_edit = edit;
    };

    for (int i = 0; i < static_cast<int>(barcodes.size()); ++i) {
        const std::string& parent = barcodes[static_cast<size_t>(i)];
        add_owner(parent, i, 0);
        if (max_edit < 1) continue;
        std::vector<std::string> level1 = mutate_barcode_once(parent);
        for (const auto& m : level1) add_owner(m, i, 1);
        if (max_edit >= 2) {
            for (const auto& m1 : level1) {
                for (const auto& m2 : mutate_barcode_once(m1)) {
                    if (m2 != parent) add_owner(m2, i, 2);
                }
            }
        }
    }

    for (const auto& kv : owners) {
        if (kv.second.parents.size() == 1) {
            MutantEntry e;
            e.parent_id = *kv.second.parents.begin();
            e.edit = kv.second.min_edit;
            dict.mutant_map.emplace(kv.first, e);
        } else {
            std::string parents;
            bool first = true;
            for (int id : kv.second.parents) {
                if (!first) parents += ",";
                first = false;
                if (id >= 0 && id < static_cast<int>(dict.parent_names.size()) &&
                    !dict.parent_names[static_cast<size_t>(id)].empty()) {
                    parents += dict.parent_names[static_cast<size_t>(id)];
                } else {
                    parents += std::to_string(id);
                }
            }
            dict.conflicts.emplace_back(kv.first, parents);
        }
    }
    return dict;
}

std::vector<AnchorHit> find_suffix_anchors(const std::string& window,
                                           const std::string& suffix,
                                           int max_edit,
                                           bool allow_revcomp,
                                           long long* n_edlib) {
    std::vector<AnchorHit> hits;

    auto consider = [&](const std::string& seq, bool flipped) {
        if (n_edlib) (*n_edlib)++;
        EdlibAlignConfig cfg = edlibNewAlignConfig(
            max_edit, EDLIB_MODE_HW, EDLIB_TASK_LOC, nullptr, 0);
        EdlibAlignResult r = edlibAlign(
            suffix.c_str(), static_cast<int>(suffix.size()),
            seq.c_str(), static_cast<int>(seq.size()),
            cfg);
        if (r.editDistance >= 0 && r.editDistance <= max_edit &&
            r.numLocations > 0 && r.startLocations != nullptr &&
            r.endLocations != nullptr) {
            for (int i = 0; i < r.numLocations; ++i) {
                AnchorHit cand;
                cand.found = true;
                cand.flipped = flipped;
                cand.edit = r.editDistance;
                cand.start = r.startLocations[i];
                cand.end = r.endLocations[i];
                cand.seq = seq;
                hits.push_back(cand);
            }
        }
        edlibFreeAlignResult(r);
    };

    consider(window, false);
    if (allow_revcomp) {
        consider(reverse_complement(window), true);
    }
    std::sort(hits.begin(), hits.end(), [](const AnchorHit& a, const AnchorHit& b) {
        if (a.edit != b.edit) return a.edit < b.edit;
        return a.start < b.start;
    });
    return hits;
}

AnchorHit find_suffix_anchor(const std::string& window,
                             const std::string& suffix,
                             int max_edit,
                             bool allow_revcomp,
                             long long* n_edlib) {
    std::vector<AnchorHit> hits =
        find_suffix_anchors(window, suffix, max_edit, allow_revcomp, n_edlib);
    if (hits.empty()) return AnchorHit();
    return hits.front();
}

std::string extract_barcode(const std::string& seq, int suffix_start0, int barcode_len,
                            int offset) {
    // barcode is immediately upstream of suffix in the searched sequence.
    int to = suffix_start0 - 1 + offset;
    int from = to - barcode_len + 1;
    if (from < 0 || to < from || to >= static_cast<int>(seq.size())) {
        return std::string();
    }
    return seq.substr(static_cast<size_t>(from), static_cast<size_t>(barcode_len));
}

BarcodeHit lookup_barcode(const std::string& barcode, const MutantDict& dict) {
    BarcodeHit hit;
    auto it = dict.mutant_map.find(barcode);
    if (it == dict.mutant_map.end()) return hit;
    hit.found = true;
    hit.parent_id = it->second.parent_id;
    hit.edit = it->second.edit;
    hit.barcode = barcode;
    return hit;
}

} // namespace miaoseq
