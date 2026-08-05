// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(Rcpp)]]

#include <Rcpp.h>
#include "demux_internal.h"

#include <algorithm>
#include <limits>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace miaoseq;

namespace {

struct PrimerHit {
    bool found = false;
    int edit = -1;
    int start = -1;
    int end = -1;
};

PrimerHit find_primer_hw(const std::string& window,
                         const std::string& primer,
                         int max_edit,
                         long long* n_edlib) {
    PrimerHit best;
    if (primer.empty() || window.empty()) return best;

    if (n_edlib) (*n_edlib)++;
    EdlibAlignConfig cfg = edlibNewAlignConfig(
        max_edit, EDLIB_MODE_HW, EDLIB_TASK_LOC, nullptr, 0);
    EdlibAlignResult r = edlibAlign(
        primer.c_str(), static_cast<int>(primer.size()),
        window.c_str(), static_cast<int>(window.size()),
        cfg);
    if (r.editDistance >= 0 && r.editDistance <= max_edit &&
        r.numLocations > 0 && r.startLocations != nullptr &&
        r.endLocations != nullptr) {
        best.found = true;
        best.edit = r.editDistance;
        best.start = r.startLocations[0];
        best.end = r.endLocations[0];
    }
    edlibFreeAlignResult(r);
    return best;
}

struct GeneOrientationScore {
    bool ok = false;
    int gene_idx = -1;
    bool flipped = false;
    int total_edit = std::numeric_limits<int>::max();
    int f_edit = -1;
    int r_edit = -1;
};

GeneOrientationScore score_gene(const std::string& front,
                                const std::string& rear,
                                int gene_idx,
                                const std::string& f_primer,
                                const std::string& r_primer,
                                const std::string& f_rc,
                                const std::string& r_rc,
                                int max_primer_edit,
                                int max_useful_total,
                                long long* n_edlib) {
    GeneOrientationScore best;
    best.gene_idx = gene_idx;
    if (max_useful_total < 0) return best;

    // native: F at front, R_RC at rear
    {
        PrimerHit hf = find_primer_hw(front, f_primer, max_primer_edit, n_edlib);
        if (hf.found && hf.edit <= max_useful_total) {
            PrimerHit hr = find_primer_hw(rear, r_rc, max_primer_edit, n_edlib);
            if (hr.found) {
                const int total = hf.edit + hr.edit;
                if (total <= max_useful_total) {
                    best.ok = true;
                    best.flipped = false;
                    best.total_edit = total;
                    best.f_edit = hf.edit;
                    best.r_edit = hr.edit;
                    max_useful_total = total;
                    if (total == 0) return best;
                }
            }
        }
    }
    // flipped: R at front, F_RC at rear
    {
        PrimerHit hf = find_primer_hw(front, r_primer, max_primer_edit, n_edlib);
        if (!hf.found || hf.edit > max_useful_total) return best;
        PrimerHit hr = find_primer_hw(rear, f_rc, max_primer_edit, n_edlib);
        if (!hr.found) return best;
        const int total = hf.edit + hr.edit;
        if (!best.ok || total < best.total_edit) {
            if (total <= max_useful_total) {
                best.ok = true;
                best.flipped = true;
                best.total_edit = total;
                best.f_edit = hf.edit;
                best.r_edit = hr.edit;
            }
        }
    }
    return best;
}

int edlib_nw_distance(const std::string& a, const std::string& b, int max_edit) {
    if (a.empty() && b.empty()) return 0;
    if (a.empty() || b.empty()) return max_edit + 1;

    EdlibAlignConfig cfg = edlibNewAlignConfig(
        max_edit, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, nullptr, 0);
    EdlibAlignResult r = edlibAlign(
        a.c_str(), static_cast<int>(a.size()),
        b.c_str(), static_cast<int>(b.size()),
        cfg);
    const int d = r.editDistance;
    edlibFreeAlignResult(r);
    if (d < 0) return max_edit + 1;
    return d;
}

struct TrimSpan {
    std::string seq;
    int start = -1; // 0-based inclusive insert start
    int end = -1;   // 0-based exclusive insert end
};

TrimSpan trim_between_primers(const std::string& seq,
                              const std::string& f_primer,
                              const std::string& r_primer,
                              int max_edit,
                              bool include_primers = false) {
    TrimSpan out;
    if (seq.empty()) return out;
    const std::string r_rc = reverse_complement(r_primer);
    PrimerHit fh = find_primer_hw(seq, f_primer, max_edit, nullptr);
    PrimerHit rh = find_primer_hw(seq, r_rc, max_edit, nullptr);
    if (!fh.found || !rh.found) return out;

    int from = -1;
    int end_excl = -1;
    if (include_primers) {
        // Keep F..R primer span (edlib endLocations are inclusive).
        if (fh.start > rh.end) return out;
        from = fh.start;
        end_excl = rh.end + 1;
    } else {
        // Insert only: after F primer through before R primer.
        if (fh.end >= rh.start) return out;
        from = fh.end + 1;
        end_excl = rh.start;
    }
    const int len = end_excl - from;
    if (len < 1) return out;
    out.start = from;
    out.end = end_excl;
    out.seq = seq.substr(static_cast<size_t>(from), static_cast<size_t>(len));
    return out;
}

} // namespace

//' Trim reads to the amplicon span defined by primer pairs (C++)
//'
//' Returns a list with trimmed `seq`, and 0-based half-open
//' coordinates `start` / `end` (NA when trim fails).
//'
//' When `include_primers` is `FALSE` (default), the span is the insert
//' between primers. When `TRUE`, the span includes both primer matches
//' (F start through R end), so only outer adapters are removed.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List trim_amplicon_insert_cpp(CharacterVector seqs,
                              std::string f_primer,
                              std::string r_primer,
                              int max_edit = 5,
                              int n_core = 1,
                              bool include_primers = false) {
    f_primer = to_upper_acgt(f_primer);
    r_primer = to_upper_acgt(r_primer);
    const int n = seqs.size();
    std::vector<std::string> v_seqs = as<std::vector<std::string>>(seqs);
    std::vector<std::string> trimmed_v(static_cast<size_t>(n));
    std::vector<int> start_v(static_cast<size_t>(n), -1);
    std::vector<int> end_v(static_cast<size_t>(n), -1);

    int cores = std::max(1, n_core);
#ifdef _OPENMP
    omp_set_num_threads(cores);
#endif

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 64)
#endif
    for (int i = 0; i < n; ++i) {
        const std::string seq = to_upper_acgt(v_seqs[static_cast<size_t>(i)]);
        TrimSpan span = trim_between_primers(
            seq, f_primer, r_primer, max_edit, include_primers);
        trimmed_v[static_cast<size_t>(i)] = span.seq;
        start_v[static_cast<size_t>(i)] = span.start;
        end_v[static_cast<size_t>(i)] = span.end;
    }

    CharacterVector out_seq(n);
    IntegerVector out_start(n, NA_INTEGER);
    IntegerVector out_end(n, NA_INTEGER);
    for (int i = 0; i < n; ++i) {
        out_seq[i] = trimmed_v[static_cast<size_t>(i)];
        if (start_v[static_cast<size_t>(i)] >= 0) {
            out_start[i] = start_v[static_cast<size_t>(i)];
            out_end[i] = end_v[static_cast<size_t>(i)];
        }
    }
    return List::create(
        Named("seq") = out_seq,
        Named("start") = out_start,
        Named("end") = out_end
    );
}

//' Find best primer hit in a sequence window (C++)
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List find_primer_hit_cpp(std::string window,
                         std::string primer,
                         int max_edit = 3) {
    window = to_upper_acgt(window);
    primer = to_upper_acgt(primer);
    PrimerHit hit = find_primer_hw(window, primer, max_edit, nullptr);
    return List::create(
        Named("found") = hit.found,
        Named("edit") = hit.edit,
        Named("start") = hit.start,
        Named("end") = hit.end
    );
}

//' Edit distance with edlib NW (C++)
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
int edlib_edit_distance_cpp(std::string query,
                            std::string target,
                            int max_edit = 8) {
    query = to_upper_acgt(query);
    target = to_upper_acgt(target);
    return edlib_nw_distance(query, target, max_edit);
}

//' Assign genes by primer pairs across reads (C++)
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List assign_genes_primers_cpp(CharacterVector seqs,
                             CharacterVector gene_ids,
                             CharacterVector f_primers,
                             CharacterVector r_primers,
                             int end_window = 80,
                             int max_primer_edit = 3,
                             int n_core = 1) {
    const int n = seqs.size();
    const int g = gene_ids.size();
    std::vector<std::string> genes = as<std::vector<std::string>>(gene_ids);
    std::vector<std::string> f_pr = as<std::vector<std::string>>(f_primers);
    std::vector<std::string> r_pr = as<std::vector<std::string>>(r_primers);
    std::vector<std::string> f_rc(static_cast<size_t>(g));
    std::vector<std::string> r_rc(static_cast<size_t>(g));

    for (int i = 0; i < g; ++i) {
        f_pr[static_cast<size_t>(i)] = to_upper_acgt(f_pr[static_cast<size_t>(i)]);
        r_pr[static_cast<size_t>(i)] = to_upper_acgt(r_pr[static_cast<size_t>(i)]);
        f_rc[static_cast<size_t>(i)] = reverse_complement(f_pr[static_cast<size_t>(i)]);
        r_rc[static_cast<size_t>(i)] = reverse_complement(r_pr[static_cast<size_t>(i)]);
    }

    std::vector<std::string> v_seqs = as<std::vector<std::string>>(seqs);

    std::vector<int> gene_idx_v(static_cast<size_t>(n), -1);
    std::vector<int> flipped_v(static_cast<size_t>(n), 0);
    std::vector<int> total_edit_v(static_cast<size_t>(n), -1);
    std::vector<int> f_edit_v(static_cast<size_t>(n), -1);
    std::vector<int> r_edit_v(static_cast<size_t>(n), -1);
    std::vector<std::string> status_v(static_cast<size_t>(n), "no_primer_hit");

    long long n_edlib = 0;
    int cores = std::max(1, n_core);
#ifdef _OPENMP
    omp_set_num_threads(cores);
#endif

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 32) reduction(+:n_edlib)
#endif
    for (int i = 0; i < n; ++i) {
        long long local_edlib = 0;
        std::string seq = to_upper_acgt(v_seqs[static_cast<size_t>(i)]);
        const int len = static_cast<int>(seq.size());
        if (len < 1) {
            status_v[static_cast<size_t>(i)] = "empty_read";
            n_edlib += local_edlib;
            continue;
        }
        const int w = std::min(end_window, len);
        const std::string front = seq.substr(0, static_cast<size_t>(w));
        const std::string rear = seq.substr(static_cast<size_t>(len - w),
                                            static_cast<size_t>(w));

        GeneOrientationScore best;
        int n_best = 0;
        int useful_cap = std::numeric_limits<int>::max();
        for (int gi = 0; gi < g; ++gi) {
            GeneOrientationScore cand = score_gene(
                front, rear, gi,
                f_pr[static_cast<size_t>(gi)],
                r_pr[static_cast<size_t>(gi)],
                f_rc[static_cast<size_t>(gi)],
                r_rc[static_cast<size_t>(gi)],
                max_primer_edit, useful_cap, &local_edlib);
            if (!cand.ok) continue;
            if (!best.ok || cand.total_edit < best.total_edit) {
                best = cand;
                n_best = 1;
                useful_cap = best.total_edit;
            } else if (best.ok && cand.total_edit == best.total_edit) {
                ++n_best;
            }
        }

        if (!best.ok) {
            status_v[static_cast<size_t>(i)] = "no_primer_hit";
            n_edlib += local_edlib;
            continue;
        }
        if (n_best > 1) {
            status_v[static_cast<size_t>(i)] = "ambiguous_gene";
            n_edlib += local_edlib;
            continue;
        }

        gene_idx_v[static_cast<size_t>(i)] = best.gene_idx;
        flipped_v[static_cast<size_t>(i)] = best.flipped ? 1 : 0;
        total_edit_v[static_cast<size_t>(i)] = best.total_edit;
        f_edit_v[static_cast<size_t>(i)] = best.f_edit;
        r_edit_v[static_cast<size_t>(i)] = best.r_edit;
        status_v[static_cast<size_t>(i)] = "assigned";
        n_edlib += local_edlib;
    }

    IntegerVector gene_idx(n, NA_INTEGER);
    LogicalVector flipped(n, false);
    IntegerVector total_edit(n, NA_INTEGER);
    IntegerVector f_edit(n, NA_INTEGER);
    IntegerVector r_edit(n, NA_INTEGER);
    CharacterVector status(n);
    CharacterVector assigned_gene(n);

    for (int i = 0; i < n; ++i) {
        status[i] = status_v[static_cast<size_t>(i)];
        if (gene_idx_v[static_cast<size_t>(i)] < 0) {
            assigned_gene[i] = NA_STRING;
            continue;
        }
        gene_idx[i] = gene_idx_v[static_cast<size_t>(i)] + 1;
        flipped[i] = flipped_v[static_cast<size_t>(i)] != 0;
        total_edit[i] = total_edit_v[static_cast<size_t>(i)];
        f_edit[i] = f_edit_v[static_cast<size_t>(i)];
        r_edit[i] = r_edit_v[static_cast<size_t>(i)];
        assigned_gene[i] = genes[static_cast<size_t>(gene_idx_v[static_cast<size_t>(i)])];
    }

    return List::create(
        Named("gene_id") = assigned_gene,
        Named("gene_idx") = gene_idx,
        Named("flipped") = flipped,
        Named("total_edit") = total_edit,
        Named("f_edit") = f_edit,
        Named("r_edit") = r_edit,
        Named("status") = status,
        Named("n_edlib") = static_cast<double>(n_edlib)
    );
}


