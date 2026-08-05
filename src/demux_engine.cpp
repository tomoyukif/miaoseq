#include "demux_engine.h"

#include <algorithm>
#include <unordered_set>

namespace miaoseq {
namespace {

struct EndHit {
    bool ok = false;
    int parent_id = -1;
    int barcode_edit = -1;
    int anchor_edit = -1;
};

bool better_hit(const EndHit& a, const EndHit& b) {
    if (!a.ok) return false;
    if (!b.ok) return true;
    if (a.barcode_edit != b.barcode_edit) return a.barcode_edit < b.barcode_edit;
    return a.anchor_edit < b.anchor_edit;
}

EndHit match_end(const std::string& window,
                 const std::string& suffix,
                 int barcode_len,
                 const MutantDict& dict,
                 int max_anchor_edit,
                 bool allow_revcomp,
                 bool* saw_suffix,
                 long long* n_edlib) {
    EndHit best;
    std::vector<AnchorHit> anchors = find_suffix_anchors(
        window, suffix, max_anchor_edit, allow_revcomp, n_edlib);
    if (!anchors.empty() && saw_suffix) *saw_suffix = true;
    if (anchors.empty()) return best;

    static const int offsets[] = {0, -1, 1, -2, 2, -3, 3, -4, 4};
    for (const AnchorHit& anchor : anchors) {
        for (int off : offsets) {
            std::string bc = extract_barcode(anchor.seq, anchor.start, barcode_len, off);
            if (bc.empty()) continue;
            BarcodeHit hit = lookup_barcode(bc, dict);
            if (!hit.found) continue;

            EndHit cand;
            cand.ok = true;
            cand.parent_id = hit.parent_id;
            cand.barcode_edit = hit.edit;
            cand.anchor_edit = anchor.edit;
            if (better_hit(cand, best)) best = cand;
            if (best.ok && best.barcode_edit == 0 && best.anchor_edit == 0) return best;
        }
    }
    return best;
}

struct ReadEnds {
    EndHit f_front, f_rear, r_front, r_rear;
    bool any_suffix = false;
};

ReadEnds scan_ends(const std::string& seq_raw,
                   const std::string& f_suffix,
                   const std::string& r_suffix,
                   int f_bc_len,
                   int r_bc_len,
                   const MutantDict& f_dict,
                   const MutantDict& r_dict,
                   int end_window,
                   int max_anchor_edit,
                   bool allow_revcomp,
                   long long* n_edlib) {
    ReadEnds ends;
    std::string seq = to_upper_acgt(seq_raw);
    const int n = static_cast<int>(seq.size());
    if (n < 1) return ends;
    const int w = std::min(end_window, n);
    std::string front = seq.substr(0, static_cast<size_t>(w));
    std::string rear = seq.substr(static_cast<size_t>(n - w), static_cast<size_t>(w));

    ends.f_front = match_end(front, f_suffix, f_bc_len, f_dict,
                             max_anchor_edit, allow_revcomp,
                             &ends.any_suffix, n_edlib);
    ends.f_rear = match_end(rear, f_suffix, f_bc_len, f_dict,
                            max_anchor_edit, allow_revcomp,
                            &ends.any_suffix, n_edlib);
    ends.r_front = match_end(front, r_suffix, r_bc_len, r_dict,
                             max_anchor_edit, allow_revcomp,
                             &ends.any_suffix, n_edlib);
    ends.r_rear = match_end(rear, r_suffix, r_bc_len, r_dict,
                            max_anchor_edit, allow_revcomp,
                            &ends.any_suffix, n_edlib);
    return ends;
}

int orientation_score(const EndHit& f, const EndHit& r) {
    return f.barcode_edit * 1000 + r.barcode_edit * 1000 +
           f.anchor_edit + r.anchor_edit;
}

int single_end_score(const EndHit& h) {
    return h.barcode_edit * 1000 + h.anchor_edit;
}

struct SingleCand {
    EndHit hit;
    bool is_f = true;
    int row = -1;
    int score = 0;
};

void consider_single_hit(const EndHit& hit,
                         bool is_f,
                         const std::unordered_map<std::string, int>& name_to_row,
                         const MutantDict& dict,
                         std::vector<SingleCand>& cands) {
    if (!hit.ok) return;
    const std::string& name =
        dict.parent_names[static_cast<size_t>(hit.parent_id)];
    auto it = name_to_row.find(name);
    if (it == name_to_row.end()) return;
    SingleCand c;
    c.hit = hit;
    c.is_f = is_f;
    c.row = it->second;
    c.score = single_end_score(hit);
    cands.push_back(c);
}

} // namespace

DemuxEngine::DemuxEngine(const std::string& f_suffix,
                         const std::string& r_suffix,
                         const std::vector<std::string>& f_barcodes,
                         const std::vector<std::string>& f_barcode_names,
                         const std::vector<std::string>& r_barcodes,
                         const std::vector<std::string>& r_barcode_names,
                         const std::vector<std::string>& pair_f_names,
                         const std::vector<std::string>& pair_r_names,
                         const std::vector<std::string>& pair_ids,
                         const std::vector<std::string>& sample_ids,
                         int end_window,
                         int max_anchor_edit,
                         int max_barcode_edit,
                         bool allow_revcomp,
                         bool allow_single_end)
    : f_suffix_(f_suffix),
      r_suffix_(r_suffix),
      f_dict_(build_mutant_dict(f_barcodes, f_barcode_names, max_barcode_edit)),
      r_dict_(build_mutant_dict(r_barcodes, r_barcode_names, max_barcode_edit)),
      end_window_(end_window),
      max_anchor_edit_(max_anchor_edit),
      allow_revcomp_(allow_revcomp),
      allow_single_end_(allow_single_end),
      pair_ids_(pair_ids),
      sample_ids_(sample_ids) {
    f_bc_len_ = f_barcodes.empty() ? 0 : static_cast<int>(f_barcodes[0].size());
    r_bc_len_ = r_barcodes.empty() ? 0 : static_cast<int>(r_barcodes[0].size());
    for (size_t i = 0; i < pair_ids_.size(); ++i) {
        pair_map_[pair_f_names[i] + "\t" + pair_r_names[i]] = static_cast<int>(i);
        f_name_to_row_[pair_f_names[i]] = static_cast<int>(i);
        r_name_to_row_[pair_r_names[i]] = static_cast<int>(i);
    }
}

DemuxOut DemuxEngine::classify(const std::string& seq, long long* n_edlib) const {
    DemuxOut out;
    ReadEnds ends = scan_ends(
        seq, f_suffix_, r_suffix_, f_bc_len_, r_bc_len_, f_dict_, r_dict_,
        end_window_, max_anchor_edit_, allow_revcomp_, n_edlib);

    struct Cand {
        bool ok = false;
        EndHit f, r;
        int row = -1;
        int score = 0;
    };

    Cand best, second;
    bool had_end_pair = false;

    auto consider = [&](const EndHit& f, const EndHit& r) {
        if (f.ok && r.ok) had_end_pair = true;
        if (!(f.ok && r.ok)) return;
        const std::string& f_name =
            f_dict_.parent_names[static_cast<size_t>(f.parent_id)];
        const std::string& r_name =
            r_dict_.parent_names[static_cast<size_t>(r.parent_id)];
        auto pit = pair_map_.find(f_name + "\t" + r_name);
        if (pit == pair_map_.end()) return;
        Cand c;
        c.ok = true;
        c.f = f;
        c.r = r;
        c.row = pit->second;
        c.score = orientation_score(f, r);
        if (!best.ok || c.score < best.score) {
            second = best;
            best = c;
        } else if (!second.ok || c.score < second.score) {
            second = c;
        }
    };

    consider(ends.f_front, ends.r_rear);
    consider(ends.f_rear, ends.r_front);

    if (best.ok && second.ok && second.score == best.score &&
        second.row != best.row) {
        out.reason = "ambiguous_pair";
        return out;
    }

    if (!best.ok) {
        if (allow_single_end_) {
            std::vector<SingleCand> single_cands;
            consider_single_hit(ends.f_front, true, f_name_to_row_, f_dict_, single_cands);
            consider_single_hit(ends.f_rear, true, f_name_to_row_, f_dict_, single_cands);
            consider_single_hit(ends.r_front, false, r_name_to_row_, r_dict_, single_cands);
            consider_single_hit(ends.r_rear, false, r_name_to_row_, r_dict_, single_cands);

            std::unordered_set<int> unique_rows;
            for (const SingleCand& c : single_cands) unique_rows.insert(c.row);

            if (unique_rows.size() == 1) {
                const int row = *unique_rows.begin();
                SingleCand best_single;
                bool have_best = false;
                for (const SingleCand& c : single_cands) {
                    if (c.row != row) continue;
                    if (!have_best || c.score < best_single.score ||
                        (c.score == best_single.score && c.is_f && !best_single.is_f)) {
                        best_single = c;
                        have_best = true;
                    }
                }

                out.status = 1;
                out.reason.clear();
                out.index_pair_id = pair_ids_[static_cast<size_t>(row)];
                out.sample_id = sample_ids_[static_cast<size_t>(row)];
                if (best_single.is_f) {
                    out.f_index_id =
                        f_dict_.parent_names[static_cast<size_t>(best_single.hit.parent_id)];
                    out.barcode_edit_f = best_single.hit.barcode_edit;
                    out.assign_mode = "single_f";
                } else {
                    out.r_index_id =
                        r_dict_.parent_names[static_cast<size_t>(best_single.hit.parent_id)];
                    out.barcode_edit_r = best_single.hit.barcode_edit;
                    out.assign_mode = "single_r";
                }
                out.match_class = (best_single.hit.barcode_edit == 0)
                    ? "complete_match" : "fuzzy_match";
                return out;
            }
            if (unique_rows.size() > 1) {
                out.reason = "ambiguous_ends";
                return out;
            }
        }

        if (had_end_pair) {
            out.reason = "invalid_pair";
        } else if (!ends.any_suffix) {
            out.reason = "no_suffix";
        } else {
            out.reason = "barcode_fail";
        }
        return out;
    }

    const int row = best.row;
    out.status = 1;
    out.reason.clear();
    out.index_pair_id = pair_ids_[static_cast<size_t>(row)];
    out.sample_id = sample_ids_[static_cast<size_t>(row)];
    out.f_index_id = f_dict_.parent_names[static_cast<size_t>(best.f.parent_id)];
    out.r_index_id = r_dict_.parent_names[static_cast<size_t>(best.r.parent_id)];
    out.barcode_edit_f = best.f.barcode_edit;
    out.barcode_edit_r = best.r.barcode_edit;
    out.match_class = (best.f.barcode_edit == 0 && best.r.barcode_edit == 0)
        ? "complete_match" : "fuzzy_match";
    out.assign_mode = "dual_end";
    return out;
}

} // namespace miaoseq
