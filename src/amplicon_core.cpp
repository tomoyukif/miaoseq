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

std::string apply_gaps(const std::string& seq,
                       const unsigned char* alignment,
                       int alignment_len,
                       bool is_query) {
    std::string out;
    out.reserve(static_cast<size_t>(alignment_len));
    int pos = 0;
    for (int i = 0; i < alignment_len; ++i) {
        const unsigned char code = alignment[i];
        if (is_query) {
            if (code == 0 || code == 2) {
                if (pos < static_cast<int>(seq.size())) {
                    out.push_back(seq[static_cast<size_t>(pos)]);
                    ++pos;
                }
            } else if (code == 1) {
                out.push_back('-');
            }
        } else {
            if (code == 0 || code == 1) {
                if (pos < static_cast<int>(seq.size())) {
                    out.push_back(seq[static_cast<size_t>(pos)]);
                    ++pos;
                }
            } else if (code == 2) {
                out.push_back('-');
            }
        }
    }
    return out;
}

struct TrimSpan {
    std::string seq;
    int start = -1; // 0-based inclusive insert start
    int end = -1;   // 0-based exclusive insert end
};

std::string apply_gaps_qual(const std::string& qual,
                            const unsigned char* alignment,
                            int alignment_len) {
    std::string out;
    out.reserve(static_cast<size_t>(alignment_len));
    int pos = 0;
    for (int i = 0; i < alignment_len; ++i) {
        const unsigned char code = alignment[i];
        if (code == 0 || code == 2) {
            if (pos < static_cast<int>(qual.size())) {
                out.push_back(qual[static_cast<size_t>(pos)]);
                ++pos;
            } else {
                out.push_back('!');
            }
        } else if (code == 1) {
            out.push_back(' '); // gap in query: no quality
        }
    }
    return out;
}

std::string majority_consensus(const std::vector<std::string>& seqs, int max_edit) {
    if (seqs.empty()) return std::string();
    if (seqs.size() == 1) return seqs[0];

    size_t ref_idx = 0;
    size_t ref_len = seqs[0].size();
    for (size_t i = 1; i < seqs.size(); ++i) {
        if (seqs[i].size() > ref_len) {
            ref_len = seqs[i].size();
            ref_idx = i;
        }
    }
    const std::string& ref = seqs[ref_idx];

    std::vector<std::string> gapped;
    gapped.reserve(seqs.size());
    size_t max_width = 0;

    for (const std::string& seq : seqs) {
        EdlibAlignConfig cfg = edlibNewAlignConfig(
            max_edit, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0);
        EdlibAlignResult r = edlibAlign(
            ref.c_str(), static_cast<int>(ref.size()),
            seq.c_str(), static_cast<int>(seq.size()),
            cfg);
        std::string g = ref;
        if (r.editDistance >= 0 && r.alignment != nullptr && r.alignmentLength > 0) {
            g = apply_gaps(seq, r.alignment, r.alignmentLength, true);
        }
        edlibFreeAlignResult(r);
        gapped.push_back(g);
        if (g.size() > max_width) max_width = g.size();
    }

    std::string out;
    out.reserve(max_width);
    for (size_t col = 0; col < max_width; ++col) {
        int counts[5] = {0, 0, 0, 0, 0}; // A,C,G,T,-
        int total = 0;
        auto bump = [&](char c) {
            switch (c) {
                case 'A': case 'a': counts[0]++; total++; break;
                case 'C': case 'c': counts[1]++; total++; break;
                case 'G': case 'g': counts[2]++; total++; break;
                case 'T': case 't': counts[3]++; total++; break;
                case '-': counts[4]++; total++; break;
                default: break;
            }
        };
        for (const std::string& g : gapped) {
            if (col < g.size()) bump(g[col]);
        }
        if (total < 1) continue;
        int best = 0;
        for (int i = 1; i < 5; ++i) {
            if (counts[i] > counts[best]) best = i;
        }
        if (counts[4] == counts[best] && counts[4] > 0) continue;
        static const char bases[] = {'A', 'C', 'G', 'T', '-'};
        if (bases[best] != '-') out.push_back(bases[best]);
    }
    return out;
}

std::string quality_consensus(const std::vector<std::string>& seqs,
                              const std::vector<std::string>& quals,
                              int max_edit) {
    if (seqs.empty()) return std::string();
    if (seqs.size() == 1) return seqs[0];
    if (quals.size() != seqs.size()) {
        return majority_consensus(seqs, max_edit);
    }

    size_t ref_idx = 0;
    size_t ref_len = seqs[0].size();
    for (size_t i = 1; i < seqs.size(); ++i) {
        if (seqs[i].size() > ref_len) {
            ref_len = seqs[i].size();
            ref_idx = i;
        }
    }
    const std::string& ref = seqs[ref_idx];

    std::vector<std::string> gapped_seq;
    std::vector<std::string> gapped_qual;
    gapped_seq.reserve(seqs.size());
    gapped_qual.reserve(seqs.size());
    size_t max_width = 0;

    for (size_t i = 0; i < seqs.size(); ++i) {
        const std::string& seq = seqs[i];
        const std::string& qual = quals[i];
        EdlibAlignConfig cfg = edlibNewAlignConfig(
            max_edit, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0);
        EdlibAlignResult r = edlibAlign(
            ref.c_str(), static_cast<int>(ref.size()),
            seq.c_str(), static_cast<int>(seq.size()),
            cfg);
        std::string gs = ref;
        std::string gq(ref.size(), '!');
        if (r.editDistance >= 0 && r.alignment != nullptr && r.alignmentLength > 0) {
            gs = apply_gaps(seq, r.alignment, r.alignmentLength, true);
            gq = apply_gaps_qual(qual, r.alignment, r.alignmentLength);
        } else if (qual.size() == seq.size()) {
            gq = qual;
        }
        edlibFreeAlignResult(r);
        gapped_seq.push_back(gs);
        gapped_qual.push_back(gq);
        if (gs.size() > max_width) max_width = gs.size();
    }

    std::string out;
    out.reserve(max_width);
    for (size_t col = 0; col < max_width; ++col) {
        double scores[5] = {0, 0, 0, 0, 0}; // A,C,G,T,-
        auto bump = [&](char c, double w) {
            switch (c) {
                case 'A': case 'a': scores[0] += w; break;
                case 'C': case 'c': scores[1] += w; break;
                case 'G': case 'g': scores[2] += w; break;
                case 'T': case 't': scores[3] += w; break;
                case '-': scores[4] += w; break;
                default: break;
            }
        };
        for (size_t i = 0; i < gapped_seq.size(); ++i) {
            if (col >= gapped_seq[i].size()) continue;
            const char c = gapped_seq[i][col];
            double w = 1.0;
            if (col < gapped_qual[i].size()) {
                const char qc = gapped_qual[i][col];
                if (qc != ' ') {
                    w = std::max(1.0, static_cast<double>(static_cast<unsigned char>(qc) - 33));
                } else {
                    w = 0.0; // gap placeholder
                }
            }
            bump(c, w > 0 ? w : 1.0);
        }
        int best = 0;
        for (int i = 1; i < 5; ++i) {
            if (scores[i] > scores[best]) best = i;
        }
        if (scores[best] <= 0) continue;
        if (best == 4) continue;
        static const char bases[] = {'A', 'C', 'G', 'T', '-'};
        out.push_back(bases[best]);
    }
    return out;
}

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


//' Majority consensus for a set of sequences (C++)
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
std::string consensus_sequences_cpp(CharacterVector seqs,
                                    int max_edit = 8) {
    std::vector<std::string> v;
    v.reserve(static_cast<size_t>(seqs.size()));
    for (int i = 0; i < seqs.size(); ++i) {
        v.push_back(to_upper_acgt(as<std::string>(seqs[i])));
    }
    return majority_consensus(v, max_edit);
}

//' Quality-weighted consensus (Phred) for a set of sequences (C++)
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
std::string consensus_sequences_quality_cpp(CharacterVector seqs,
                                            CharacterVector quals,
                                            int max_edit = 8) {
    if (quals.size() != seqs.size()) {
        stop("seqs and quals must have equal length.");
    }
    std::vector<std::string> v;
    std::vector<std::string> q;
    v.reserve(static_cast<size_t>(seqs.size()));
    q.reserve(static_cast<size_t>(quals.size()));
    for (int i = 0; i < seqs.size(); ++i) {
        v.push_back(to_upper_acgt(as<std::string>(seqs[i])));
        q.push_back(as<std::string>(quals[i]));
    }
    return quality_consensus(v, q, max_edit);
}

//' Assign genes by amplicon reference edit distance (C++)
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List assign_genes_amplicon_ref_cpp(CharacterVector seqs,
                                   CharacterVector ref_ids,
                                   CharacterVector ref_seqs,
                                   int max_edit = 12,
                                   int n_core = 1) {
    const int n = seqs.size();
    const int r = ref_ids.size();
    if (ref_seqs.size() != r) {
        stop("ref_ids and ref_seqs must have equal length.");
    }

    std::vector<std::string> v_seqs = as<std::vector<std::string>>(seqs);
    std::vector<std::string> v_refs(static_cast<size_t>(r));
    std::vector<std::string> v_ids = as<std::vector<std::string>>(ref_ids);
    for (int i = 0; i < r; ++i) {
        v_refs[static_cast<size_t>(i)] = to_upper_acgt(as<std::string>(ref_seqs[i]));
    }

    std::vector<std::string> gene_v(static_cast<size_t>(n));
    std::vector<std::string> status_v(static_cast<size_t>(n), "no_primer_hit");
    std::vector<int> edit_v(static_cast<size_t>(n), -1);

    int cores = std::max(1, n_core);
#ifdef _OPENMP
    omp_set_num_threads(cores);
#endif

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 16)
#endif
    for (int i = 0; i < n; ++i) {
        const std::string seq = to_upper_acgt(v_seqs[static_cast<size_t>(i)]);
        int best = max_edit + 1;
        int n_best = 0;
        int best_idx = -1;
        for (int j = 0; j < r; ++j) {
            const int d = edlib_nw_distance(
                seq, v_refs[static_cast<size_t>(j)], max_edit);
            if (d > max_edit) continue;
            if (d < best) {
                best = d;
                best_idx = j;
                n_best = 1;
            } else if (d == best) {
                ++n_best;
            }
        }
        if (best_idx < 0) {
            gene_v[static_cast<size_t>(i)] = "";
            status_v[static_cast<size_t>(i)] = "no_primer_hit";
            edit_v[static_cast<size_t>(i)] = -1;
        } else if (n_best > 1) {
            gene_v[static_cast<size_t>(i)] = "";
            status_v[static_cast<size_t>(i)] = "ambiguous_amplicon";
            edit_v[static_cast<size_t>(i)] = best;
        } else {
            gene_v[static_cast<size_t>(i)] = v_ids[static_cast<size_t>(best_idx)];
            status_v[static_cast<size_t>(i)] = "amplicon_fallback";
            edit_v[static_cast<size_t>(i)] = best;
        }
    }

    CharacterVector gene_id(n);
    CharacterVector status(n);
    IntegerVector total_edit(n, NA_INTEGER);
    for (int i = 0; i < n; ++i) {
        if (gene_v[static_cast<size_t>(i)].empty()) {
            gene_id[i] = NA_STRING;
        } else {
            gene_id[i] = gene_v[static_cast<size_t>(i)];
        }
        status[i] = status_v[static_cast<size_t>(i)];
        if (edit_v[static_cast<size_t>(i)] >= 0) {
            total_edit[i] = edit_v[static_cast<size_t>(i)];
        }
    }
    return List::create(
        Named("gene_id") = gene_id,
        Named("status") = status,
        Named("total_edit") = total_edit
    );
}

//' Greedy cluster + consensus for one gene bucket (C++)
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
DataFrame cluster_consensus_gene_cpp(CharacterVector seqs,
                                     std::string sample_id,
                                     std::string gene_id,
                                     std::string method = "both",
                                     int min_cluster_reads = 5,
                                     int max_clusters = 20,
                                     int max_edit = 12,
                                     Nullable<CharacterVector> quals = R_NilValue,
                                     std::string consensus_backend = "majority") {
    const int n = seqs.size();
    if (n < 1) {
        return DataFrame::create();
    }

    std::vector<std::string> v_seqs(static_cast<size_t>(n));
    std::vector<std::string> v_quals;
    const bool use_quality = (consensus_backend == "quality");
    if (use_quality) {
        if (quals.isNull()) {
            stop("consensus_backend='quality' requires quals.");
        }
        CharacterVector qv = quals.get();
        if (qv.size() != n) {
            stop("quals must have the same length as seqs.");
        }
        v_quals.resize(static_cast<size_t>(n));
        for (int i = 0; i < n; ++i) {
            v_quals[static_cast<size_t>(i)] = as<std::string>(qv[i]);
        }
    }

    for (int i = 0; i < n; ++i) {
        v_seqs[static_cast<size_t>(i)] = to_upper_acgt(as<std::string>(seqs[i]));
    }

    std::vector<int> cluster_of(static_cast<size_t>(n), -1);
    std::vector<std::string> reps;
    reps.reserve(static_cast<size_t>(max_clusters));

    for (int i = 0; i < n; ++i) {
        const std::string& seq = v_seqs[static_cast<size_t>(i)];
        int best_cluster = -1;
        int best_dist = max_edit + 1;
        const int n_reps = static_cast<int>(reps.size());
        for (int c = 0; c < n_reps; ++c) {
            const int d = edlib_nw_distance(seq, reps[static_cast<size_t>(c)], max_edit);
            if (d <= max_edit && d < best_dist) {
                best_dist = d;
                best_cluster = c;
            }
        }
        if (best_cluster >= 0) {
            cluster_of[static_cast<size_t>(i)] = best_cluster;
        } else if (n_reps < max_clusters) {
            reps.push_back(seq);
            cluster_of[static_cast<size_t>(i)] = n_reps;
        }
    }

    const int n_rep = static_cast<int>(reps.size());
    std::vector<int> counts(static_cast<size_t>(n_rep), 0);
    for (int i = 0; i < n; ++i) {
        if (cluster_of[static_cast<size_t>(i)] >= 0) {
            counts[static_cast<size_t>(cluster_of[static_cast<size_t>(i)])]++;
        }
    }

    std::vector<int> order;
    order.reserve(static_cast<size_t>(n_rep));
    for (int c = 0; c < n_rep; ++c) {
        if (counts[static_cast<size_t>(c)] >= min_cluster_reads) {
            order.push_back(c);
        }
    }
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        return counts[static_cast<size_t>(a)] > counts[static_cast<size_t>(b)];
    });
    if (order.empty()) {
        return DataFrame::create();
    }

    int n_gene_reads = 0;
    for (int c : order) n_gene_reads += counts[static_cast<size_t>(c)];

    const bool do_consensus = (method == "both" || method == "consensus");
    const int k = static_cast<int>(order.size());

    CharacterVector out_sample(k);
    CharacterVector out_gene(k);
    IntegerVector out_cid(k);
    CharacterVector out_seq(k);
    IntegerVector out_n(k);
    IntegerVector out_n_gene(k);
    NumericVector out_frac(k);
    NumericVector out_frac_sample(k);
    CharacterVector out_method(k);

    for (int j = 0; j < k; ++j) {
        const int c = order[static_cast<size_t>(j)];
        std::vector<std::string> members;
        std::vector<std::string> member_quals;
        members.reserve(static_cast<size_t>(counts[static_cast<size_t>(c)]));
        if (use_quality) {
            member_quals.reserve(static_cast<size_t>(counts[static_cast<size_t>(c)]));
        }
        for (int i = 0; i < n; ++i) {
            if (cluster_of[static_cast<size_t>(i)] == c) {
                members.push_back(v_seqs[static_cast<size_t>(i)]);
                if (use_quality) {
                    member_quals.push_back(v_quals[static_cast<size_t>(i)]);
                }
            }
        }
        std::string cons;
        std::string method_label;
        if (!do_consensus) {
            cons = members.front();
            method_label = "greedy";
        } else if (use_quality) {
            cons = quality_consensus(members, member_quals, max_edit);
            method_label = "quality";
        } else {
            cons = majority_consensus(members, max_edit);
            method_label = "majority";
        }

        out_sample[j] = sample_id;
        out_gene[j] = gene_id;
        out_cid[j] = j + 1;
        out_seq[j] = cons;
        out_n[j] = counts[static_cast<size_t>(c)];
        out_n_gene[j] = n_gene_reads;
        out_frac[j] = static_cast<double>(counts[static_cast<size_t>(c)]) /
                      static_cast<double>(n_gene_reads);
        out_frac_sample[j] = static_cast<double>(counts[static_cast<size_t>(c)]) /
                             static_cast<double>(n);
        out_method[j] = method_label;
    }

    return DataFrame::create(
        Named("sample_id") = out_sample,
        Named("gene_id") = out_gene,
        Named("cluster_id") = out_cid,
        Named("seq") = out_seq,
        Named("n_reads") = out_n,
        Named("n_reads_gene") = out_n_gene,
        Named("fraction") = out_frac,
        Named("fraction_sample") = out_frac_sample,
        Named("method") = out_method
    );
}
