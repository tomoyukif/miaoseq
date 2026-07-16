// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(Rcpp)]]

#include <Rcpp.h>
#include "demux_internal.h"

#include <algorithm>
#include <unordered_set>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace miaoseq;

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

//' Build barcode mutant dictionary (C++)
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List build_barcode_dict_cpp(CharacterVector barcodes,
                            CharacterVector names,
                            int max_barcode_edit = 2) {
    std::vector<std::string> bcs = as<std::vector<std::string>>(barcodes);
    std::vector<std::string> nms = as<std::vector<std::string>>(names);
    MutantDict d = build_mutant_dict(bcs, nms, max_barcode_edit);

    CharacterVector mutants(d.mutant_map.size());
    IntegerVector parent_ids(d.mutant_map.size());
    IntegerVector edits(d.mutant_map.size());
    int i = 0;
    for (const auto& kv : d.mutant_map) {
        mutants[i] = kv.first;
        parent_ids[i] = kv.second.parent_id + 1;
        edits[i] = kv.second.edit;
        ++i;
    }

    CharacterVector conflict_mut(d.conflicts.size());
    CharacterVector conflict_par(d.conflicts.size());
    for (size_t j = 0; j < d.conflicts.size(); ++j) {
        conflict_mut[j] = d.conflicts[j].first;
        conflict_par[j] = d.conflicts[j].second;
    }

    return List::create(
        Named("parents") = barcodes,
        Named("parent_names") = names,
        Named("mutants") = mutants,
        Named("parent_ids") = parent_ids,
        Named("edits") = edits,
        Named("conflict_mutants") = conflict_mut,
        Named("conflict_parents") = conflict_par,
        Named("n_mutants") = static_cast<int>(d.mutant_map.size()),
        Named("n_conflicts") = static_cast<int>(d.conflicts.size())
    );
}

//' Demultiplex reads with C++ edlib core
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List demux_reads_cpp(CharacterVector seqs,
                     CharacterVector ids,
                     std::string f_suffix,
                     std::string r_suffix,
                     CharacterVector f_barcodes,
                     CharacterVector f_barcode_names,
                     CharacterVector r_barcodes,
                     CharacterVector r_barcode_names,
                     CharacterVector pair_f_names,
                     CharacterVector pair_r_names,
                     CharacterVector pair_ids,
                     CharacterVector sample_ids,
                     int end_window = 120,
                     int max_anchor_edit = 10,
                     int max_barcode_edit = 2,
                     bool allow_revcomp = true,
                     bool allow_single_end = false,
                     int n_core = 1) {
    const int n = seqs.size();
    MutantDict f_dict = build_mutant_dict(
        as<std::vector<std::string>>(f_barcodes),
        as<std::vector<std::string>>(f_barcode_names),
        max_barcode_edit);
    MutantDict r_dict = build_mutant_dict(
        as<std::vector<std::string>>(r_barcodes),
        as<std::vector<std::string>>(r_barcode_names),
        max_barcode_edit);

    const int f_bc_len = f_barcodes.size() > 0 ?
        static_cast<int>(std::string(f_barcodes[0]).size()) : 0;
    const int r_bc_len = r_barcodes.size() > 0 ?
        static_cast<int>(std::string(r_barcodes[0]).size()) : 0;

    std::unordered_map<std::string, int> pair_map;
    std::unordered_map<std::string, int> f_name_to_row;
    std::unordered_map<std::string, int> r_name_to_row;
    for (int i = 0; i < pair_ids.size(); ++i) {
        pair_map[std::string(pair_f_names[i]) + "\t" +
                 std::string(pair_r_names[i])] = i;
        f_name_to_row[std::string(pair_f_names[i])] = i;
        r_name_to_row[std::string(pair_r_names[i])] = i;
    }

    std::vector<std::string> v_seqs = as<std::vector<std::string>>(seqs);

    std::vector<int> status(n, 0);
    std::vector<std::string> reason(n, "no_suffix");
    std::vector<std::string> out_pair(n), out_sample(n), out_f(n), out_r(n);
    std::vector<int> be_f(n, -1), be_r(n, -1);
    std::vector<std::string> match_class(n), assign_mode(n);
    long long n_edlib = 0;

    int cores = std::max(1, n_core);
#ifdef _OPENMP
    omp_set_num_threads(cores);
#endif

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 64) reduction(+:n_edlib)
#endif
    for (int i = 0; i < n; ++i) {
        long long local_edlib = 0;
        ReadEnds ends = scan_ends(
            v_seqs[static_cast<size_t>(i)],
            f_suffix, r_suffix,
            f_bc_len, r_bc_len, f_dict, r_dict,
            end_window, max_anchor_edit,
            allow_revcomp, &local_edlib);

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
                f_dict.parent_names[static_cast<size_t>(f.parent_id)];
            const std::string& r_name =
                r_dict.parent_names[static_cast<size_t>(r.parent_id)];
            auto pit = pair_map.find(f_name + "\t" + r_name);
            if (pit == pair_map.end()) return;
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
            reason[static_cast<size_t>(i)] = "ambiguous_pair";
            n_edlib += local_edlib;
            continue;
        }

        if (!best.ok) {
            bool rescued = false;
            if (allow_single_end) {
                std::vector<SingleCand> single_cands;
                consider_single_hit(ends.f_front, true, f_name_to_row, f_dict, single_cands);
                consider_single_hit(ends.f_rear, true, f_name_to_row, f_dict, single_cands);
                consider_single_hit(ends.r_front, false, r_name_to_row, r_dict, single_cands);
                consider_single_hit(ends.r_rear, false, r_name_to_row, r_dict, single_cands);

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

                    status[static_cast<size_t>(i)] = 1;
                    reason[static_cast<size_t>(i)] = "";
                    out_pair[static_cast<size_t>(i)] = std::string(pair_ids[row]);
                    out_sample[static_cast<size_t>(i)] = std::string(sample_ids[row]);
                    if (best_single.is_f) {
                        out_f[static_cast<size_t>(i)] =
                            f_dict.parent_names[static_cast<size_t>(best_single.hit.parent_id)];
                        be_f[static_cast<size_t>(i)] = best_single.hit.barcode_edit;
                        assign_mode[static_cast<size_t>(i)] = "single_f";
                    } else {
                        out_r[static_cast<size_t>(i)] =
                            r_dict.parent_names[static_cast<size_t>(best_single.hit.parent_id)];
                        be_r[static_cast<size_t>(i)] = best_single.hit.barcode_edit;
                        assign_mode[static_cast<size_t>(i)] = "single_r";
                    }
                    if (best_single.hit.barcode_edit == 0) {
                        match_class[static_cast<size_t>(i)] = "complete_match";
                    } else {
                        match_class[static_cast<size_t>(i)] = "fuzzy_match";
                    }
                    rescued = true;
                } else if (unique_rows.size() > 1) {
                    reason[static_cast<size_t>(i)] = "ambiguous_ends";
                    rescued = true;
                }
            }

            if (!rescued) {
                if (had_end_pair) {
                    reason[static_cast<size_t>(i)] = "invalid_pair";
                } else if (!ends.any_suffix) {
                    reason[static_cast<size_t>(i)] = "no_suffix";
                } else {
                    reason[static_cast<size_t>(i)] = "barcode_fail";
                }
            }
            n_edlib += local_edlib;
            continue;
        }

        int row = best.row;
        status[static_cast<size_t>(i)] = 1;
        reason[static_cast<size_t>(i)] = "";
        out_pair[static_cast<size_t>(i)] = std::string(pair_ids[row]);
        out_sample[static_cast<size_t>(i)] = std::string(sample_ids[row]);
        out_f[static_cast<size_t>(i)] =
            f_dict.parent_names[static_cast<size_t>(best.f.parent_id)];
        out_r[static_cast<size_t>(i)] =
            r_dict.parent_names[static_cast<size_t>(best.r.parent_id)];
        be_f[static_cast<size_t>(i)] = best.f.barcode_edit;
        be_r[static_cast<size_t>(i)] = best.r.barcode_edit;
        if (best.f.barcode_edit == 0 && best.r.barcode_edit == 0) {
            match_class[static_cast<size_t>(i)] = "complete_match";
        } else {
            match_class[static_cast<size_t>(i)] = "fuzzy_match";
        }
        assign_mode[static_cast<size_t>(i)] = "dual_end";
        n_edlib += local_edlib;
    }

    return List::create(
        Named("status") = status,
        Named("reason") = reason,
        Named("read_id") = ids,
        Named("index_pair_id") = out_pair,
        Named("sample_id") = out_sample,
        Named("f_index_id") = out_f,
        Named("r_index_id") = out_r,
        Named("barcode_edit_f") = be_f,
        Named("barcode_edit_r") = be_r,
        Named("match_class") = match_class,
        Named("assign_mode") = assign_mode,
        Named("n_edlib") = static_cast<double>(n_edlib)
    );
}
