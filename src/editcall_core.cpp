// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(Rcpp)]]

#include <Rcpp.h>
#include "demux_internal.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

namespace {

struct AlnPack {
    std::string ref_aln;
    std::string query_aln;
    bool ok = false;
};

struct WinSeq {
    std::string read_seq;
    std::string ref_seq;
    bool ok = false;
};

struct PamRow {
    std::string target_gene;
    std::string guide_id;
    int cut_insert = 0;
    int win_start = 0;
    int win_end = 0;
};

struct JointEv {
    std::string guide_i;
    std::string guide_j;
    std::string target_gene_i;
    std::string target_gene_j;
    std::string event_class;
    int del_span = 0;
    int expected_span = 0;
    WinSeq junction;
};

struct EditRow {
    std::string sample_id;
    std::string index_pair_id;
    std::string target_gene;
    std::string read_seq;
    std::string ref_seq;
    bool intact = false;
    std::string allele_source;
};

struct JointRow {
    std::string sample_id;
    std::string index_pair_id;
    std::string read_id;
    std::string gene_id;
    std::string guide_i;
    std::string guide_j;
    std::string target_gene_i;
    std::string target_gene_j;
    std::string event_class;
    int del_span = 0;
    int expected_span = 0;
};

std::string ungapped(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
        if (c != '-') out.push_back(c);
    }
    return out;
}

AlnPack nw_gapped_alignment(const std::string& query, const std::string& target) {
    AlnPack out;
    if (query.empty() || target.empty()) return out;
    EdlibAlignConfig cfg = edlibNewAlignConfig(
        -1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0);
    EdlibAlignResult r = edlibAlign(
        query.c_str(), static_cast<int>(query.size()),
        target.c_str(), static_cast<int>(target.size()),
        cfg);
    if (r.status != EDLIB_STATUS_OK || r.alignment == nullptr ||
        r.alignmentLength < 1) {
        edlibFreeAlignResult(r);
        return out;
    }
    out.ref_aln.reserve(static_cast<size_t>(r.alignmentLength));
    out.query_aln.reserve(static_cast<size_t>(r.alignmentLength));
    int qi = 0;
    int ti = 0;
    for (int i = 0; i < r.alignmentLength; ++i) {
        const unsigned char code = r.alignment[i];
        // edlib NW path (empirical): INSERT advances query / gaps target;
        // DELETE advances target / gaps query.
        if (code == EDLIB_EDOP_MATCH || code == EDLIB_EDOP_MISMATCH) {
            if (qi >= static_cast<int>(query.size()) ||
                ti >= static_cast<int>(target.size())) break;
            out.query_aln.push_back(query[static_cast<size_t>(qi++)]);
            out.ref_aln.push_back(target[static_cast<size_t>(ti++)]);
        } else if (code == EDLIB_EDOP_INSERT) {
            if (qi >= static_cast<int>(query.size())) break;
            out.query_aln.push_back(query[static_cast<size_t>(qi++)]);
            out.ref_aln.push_back('-');
        } else if (code == EDLIB_EDOP_DELETE) {
            if (ti >= static_cast<int>(target.size())) break;
            out.query_aln.push_back('-');
            out.ref_aln.push_back(target[static_cast<size_t>(ti++)]);
        }
    }
    edlibFreeAlignResult(r);
    out.ok = (out.ref_aln.size() == out.query_aln.size() &&
              !out.ref_aln.empty());
    return out;
}

struct Bounds {
    int start = 0; // 1-based inclusive
    int end = 0;
};

Bounds gapped_window_bounds(const std::string& ref_aln, int start_u, int end_u) {
    const int target_len = end_u - start_u;
    int n_ins = 0;
    const int prefix_len = start_u - 1;
    for (int i = 0; i < prefix_len && i < static_cast<int>(ref_aln.size()); ++i) {
        if (ref_aln[static_cast<size_t>(i)] == '-') ++n_ins;
    }
    const int target_start = start_u + n_ins;
    std::string target_rest;
    if (target_start >= 1 && target_start <= static_cast<int>(ref_aln.size())) {
        target_rest = ref_aln.substr(static_cast<size_t>(target_start - 1));
    }
    int target_end = target_len;
    int detected_ins = 0;
    while (true) {
        int n = 0;
        const int lim = std::min(target_end, static_cast<int>(target_rest.size()));
        for (int i = 0; i < lim; ++i) {
            if (target_rest[static_cast<size_t>(i)] == '-') ++n;
        }
        if (n > detected_ins) {
            target_end += (n - detected_ins);
            detected_ins = n;
        } else {
            break;
        }
    }
    target_end = target_start + target_end;
    Bounds b;
    b.start = target_start;
    b.end = target_end;
    return b;
}

bool alignment_end_has_anchor(const std::string& ref_aln,
                              const std::string& query_aln,
                              int gstart,
                              int gend,
                              int anchor_bp,
                              bool left) {
    if (gend < gstart || gstart < 1 ||
        gend > static_cast<int>(ref_aln.size()) ||
        gend > static_cast<int>(query_aln.size())) {
        return false;
    }
    const int n = gend - gstart + 1;
    if (n < anchor_bp) return false;
    const int from = left ? 0 : (n - anchor_bp);
    for (int k = 0; k < anchor_bp; ++k) {
        const size_t idx = static_cast<size_t>(gstart - 1 + from + k);
        const char rc = ref_aln[idx];
        const char qc = query_aln[idx];
        if (rc == '-' || qc == '-' || rc != qc) return false;
    }
    return true;
}

WinSeq extract_adaptive_edit_window(const std::string& ref_aln,
                                    const std::string& query_aln,
                                    int left_u,
                                    int right_u,
                                    int anchor_bp,
                                    int max_expand) {
    WinSeq out;
    if (right_u < left_u) return out;
    const int ref_ungapped_len = static_cast<int>(ungapped(ref_aln).size());
    const int min_left = std::max(1, left_u - max_expand);
    const int max_right = std::min(ref_ungapped_len, right_u + max_expand);
    while (true) {
        Bounds b = gapped_window_bounds(ref_aln, left_u, right_u);
        if (b.end < b.start || b.start < 1 ||
            b.end > static_cast<int>(query_aln.size())) {
            return out;
        }
        const bool left_ok = alignment_end_has_anchor(
            ref_aln, query_aln, b.start, b.end, anchor_bp, true);
        const bool right_ok = alignment_end_has_anchor(
            ref_aln, query_aln, b.start, b.end, anchor_bp, false);
        if (left_ok && right_ok) {
            out.read_seq = query_aln.substr(
                static_cast<size_t>(b.start - 1),
                static_cast<size_t>(b.end - b.start + 1));
            out.ref_seq = ref_aln.substr(
                static_cast<size_t>(b.start - 1),
                static_cast<size_t>(b.end - b.start + 1));
            if (ungapped(out.read_seq).empty()) return WinSeq();
            out.ok = true;
            return out;
        }
        bool expanded = false;
        if (!left_ok) {
            if (left_u <= min_left) return WinSeq();
            --left_u;
            expanded = true;
        }
        if (!right_ok) {
            if (right_u >= max_right) return WinSeq();
            ++right_u;
            expanded = true;
        }
        if (!expanded) return WinSeq();
    }
}

int count_query_deletions(const std::string& ref_aln,
                          const std::string& query_aln,
                          int start_u,
                          int end_u) {
    if (end_u < start_u) return 0;
    Bounds b = gapped_window_bounds(ref_aln, start_u, end_u);
    if (b.end < b.start || b.start < 1 ||
        b.end > static_cast<int>(ref_aln.size()) ||
        b.end > static_cast<int>(query_aln.size())) {
        return 0;
    }
    int n = 0;
    for (int i = b.start; i <= b.end; ++i) {
        const char r = ref_aln[static_cast<size_t>(i - 1)];
        const char q = query_aln[static_cast<size_t>(i - 1)];
        if (r != '-' && q == '-') ++n;
    }
    return n;
}

WinSeq local_window_for_pam(const AlnPack& aln, const PamRow& pam,
                            int anchor_bp, int max_expand) {
    const int ref_len = static_cast<int>(ungapped(aln.ref_aln).size());
    const int start = std::max(1, pam.win_start);
    const int end = std::min(ref_len, pam.win_end);
    if (end < start) return WinSeq();
    return extract_adaptive_edit_window(
        aln.ref_aln, aln.query_aln, start, end, anchor_bp, max_expand);
}

std::vector<JointEv> classify_joint_events(
    std::vector<PamRow> pam_rows,
    const std::unordered_map<std::string, WinSeq>& local_wins,
    const std::string& ref_aln,
    const std::string& query_aln,
    int check_window,
    int anchor_bp,
    int max_expand,
    int min_span_bp,
    int excision_tol_bp) {
    std::vector<JointEv> events;
    if (pam_rows.size() < 2) return events;
    std::sort(pam_rows.begin(), pam_rows.end(),
              [](const PamRow& a, const PamRow& b) {
                  if (a.cut_insert != b.cut_insert) return a.cut_insert < b.cut_insert;
                  if (a.guide_id != b.guide_id) return a.guide_id < b.guide_id;
                  return a.target_gene < b.target_gene;
              });
    const int ref_len = static_cast<int>(ungapped(ref_aln).size());
    for (size_t k = 0; k + 1 < pam_rows.size(); ++k) {
        const PamRow& row_i = pam_rows[k];
        const PamRow& row_j = pam_rows[k + 1];
        const std::string& tg_i = row_i.target_gene;
        const std::string& tg_j = row_j.target_gene;
        const std::string guide_i = !row_i.guide_id.empty() ? row_i.guide_id : tg_i;
        const std::string guide_j = !row_j.guide_id.empty() ? row_j.guide_id : tg_j;
        const int c_i = row_i.cut_insert;
        const int c_j = row_j.cut_insert;
        const int expected = c_j - c_i;
        const int del_span = count_query_deletions(ref_aln, query_aln, c_i, c_j);
        const int left_u = std::max(1, c_i - check_window);
        int right_u = c_j + check_window;
        right_u = std::min(ref_len, right_u);
        WinSeq junction;
        if (right_u >= left_u) {
            junction = extract_adaptive_edit_window(
                ref_aln, query_aln, left_u, right_u, anchor_bp, max_expand);
        }
        const bool outer_ok = junction.ok;
        WinSeq win_i;
        WinSeq win_j;
        auto it_i = local_wins.find(tg_i);
        auto it_j = local_wins.find(tg_j);
        if (it_i != local_wins.end()) win_i = it_i->second;
        if (it_j != local_wins.end()) win_j = it_j->second;
        const bool i_ok = win_i.ok && !win_i.read_seq.empty();
        const bool j_ok = win_j.ok && !win_j.read_seq.empty();
        const bool i_wt = i_ok && win_i.read_seq == win_i.ref_seq;
        const bool j_wt = j_ok && win_j.read_seq == win_j.ref_seq;

        std::string event_class;
        if (outer_ok && expected >= min_span_bp &&
            std::abs(del_span - expected) <= excision_tol_bp) {
            event_class = "both_cut_excision";
        } else if (!i_ok && !j_ok && !outer_ok) {
            continue;
        } else if (i_wt && j_wt) {
            event_class = "wt";
        } else if (i_ok && !i_wt && (!j_ok || j_wt)) {
            event_class = "g_i_only";
        } else if (j_ok && !j_wt && (!i_ok || i_wt)) {
            event_class = "g_j_only";
        } else if (i_ok && j_ok && !i_wt && !j_wt) {
            event_class = "both_local";
        } else if (outer_ok && (i_ok || j_ok)) {
            event_class = "both_local";
        } else {
            continue;
        }

        JointEv ev;
        ev.guide_i = guide_i;
        ev.guide_j = guide_j;
        ev.target_gene_i = tg_i;
        ev.target_gene_j = tg_j;
        ev.event_class = event_class;
        ev.del_span = del_span;
        ev.expected_span = expected;
        if (event_class == "both_cut_excision") ev.junction = junction;
        events.push_back(std::move(ev));
    }
    return events;
}

std::unordered_map<std::string, WinSeq> resolve_excision_plan_a_alleles(
    const std::vector<JointEv>& events) {
    std::unordered_map<std::string, WinSeq> allele;
    std::unordered_map<std::string, int> n_hit;
    for (const JointEv& ev : events) {
        if (ev.event_class != "both_cut_excision" || !ev.junction.ok) continue;
        const std::string tgs[2] = {ev.target_gene_i, ev.target_gene_j};
        for (const std::string& tg : tgs) {
            const int n = n_hit.count(tg) ? n_hit[tg] : 0;
            n_hit[tg] = n + 1;
            if (n >= 1) {
                WinSeq tok;
                tok.read_seq = "---";
                tok.ref_seq = "---";
                tok.ok = true;
                allele[tg] = tok;
            } else {
                allele[tg] = ev.junction;
            }
        }
    }
    return allele;
}

std::string as_upper_acgt_keep_other(const std::string& s) {
    std::string out = s;
    for (char& c : out) {
        if (c >= 'a' && c <= 'z') c = static_cast<char>(c - 'a' + 'A');
    }
    return out;
}

} // namespace

//' Global NW gapped alignment via edlib (query vs target/ref).
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List edlib_nw_align_cpp(std::string query, std::string target) {
    AlnPack aln = nw_gapped_alignment(
        as_upper_acgt_keep_other(query),
        as_upper_acgt_keep_other(target));
    return List::create(
        Named("ok") = aln.ok,
        Named("ref_aln") = aln.ref_aln,
        Named("query_aln") = aln.query_aln);
}

//' Map ungapped ref coordinates to gapped alignment bounds (1-based inclusive).
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
IntegerVector gapped_window_bounds_cpp(std::string ref_aln,
                                       int start_ungapped,
                                       int end_ungapped) {
    Bounds b = gapped_window_bounds(ref_aln, start_ungapped, end_ungapped);
    IntegerVector out = IntegerVector::create(
        Named("start") = b.start, Named("end") = b.end);
    return out;
}

//' Adaptive edit-window extraction (C++ port of the R helper).
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List extract_adaptive_edit_window_cpp(std::string ref_aln,
                                      std::string query_aln,
                                      int start,
                                      int end,
                                      int anchor_bp = 5,
                                      int max_expand = 50) {
    WinSeq win = extract_adaptive_edit_window(
        ref_aln, query_aln, start, end, anchor_bp, max_expand);
    if (!win.ok) return R_NilValue;
    return List::create(
        Named("read_seq") = win.read_seq,
        Named("ref_seq") = win.ref_seq);
}

//' Per-read editcall extraction (trim assumed done): align + windows + joint.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List editcall_process_reads_cpp(CharacterVector sample_id,
                                CharacterVector index_pair_id,
                                CharacterVector read_id,
                                CharacterVector gene_id,
                                CharacterVector seqs,
                                CharacterVector pam_gene,
                                CharacterVector pam_target_gene,
                                CharacterVector pam_guide_id,
                                IntegerVector pam_cut_insert,
                                IntegerVector pam_win_start,
                                IntegerVector pam_win_end,
                                CharacterVector ref_gene_id,
                                CharacterVector ref_seq,
                                int check_window = 10,
                                int anchor_bp = 5,
                                int max_expand = 50,
                                int min_span_bp = 30,
                                int excision_tol_bp = 20,
                                int n_core = 1) {
    const int n = seqs.size();
    std::unordered_map<std::string, std::string> refs;
    refs.reserve(static_cast<size_t>(ref_gene_id.size()));
    for (int i = 0; i < ref_gene_id.size(); ++i) {
        if (ref_gene_id[i] == NA_STRING || ref_seq[i] == NA_STRING) continue;
        refs[as<std::string>(ref_gene_id[i])] =
            as_upper_acgt_keep_other(as<std::string>(ref_seq[i]));
    }

    std::unordered_map<std::string, std::vector<PamRow>> pams;
    for (int i = 0; i < pam_gene.size(); ++i) {
        if (pam_gene[i] == NA_STRING || pam_target_gene[i] == NA_STRING) continue;
        PamRow row;
        row.target_gene = as<std::string>(pam_target_gene[i]);
        if (pam_guide_id[i] != NA_STRING) {
            row.guide_id = as<std::string>(pam_guide_id[i]);
            if (row.guide_id == "NA") row.guide_id.clear();
        }
        row.cut_insert = pam_cut_insert[i];
        row.win_start = pam_win_start[i];
        row.win_end = pam_win_end[i];
        pams[as<std::string>(pam_gene[i])].push_back(std::move(row));
    }

    std::vector<std::string> v_sample(static_cast<size_t>(n));
    std::vector<std::string> v_pair(static_cast<size_t>(n));
    std::vector<std::string> v_read(static_cast<size_t>(n));
    std::vector<std::string> v_gene(static_cast<size_t>(n));
    std::vector<std::string> v_seq(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        v_sample[static_cast<size_t>(i)] =
            sample_id[i] == NA_STRING ? "" : as<std::string>(sample_id[i]);
        v_pair[static_cast<size_t>(i)] =
            index_pair_id[i] == NA_STRING ? "" : as<std::string>(index_pair_id[i]);
        v_read[static_cast<size_t>(i)] =
            read_id[i] == NA_STRING ? "" : as<std::string>(read_id[i]);
        v_gene[static_cast<size_t>(i)] =
            gene_id[i] == NA_STRING ? "" : as<std::string>(gene_id[i]);
        v_seq[static_cast<size_t>(i)] =
            seqs[i] == NA_STRING ? "" :
            as_upper_acgt_keep_other(as<std::string>(seqs[i]));
    }

    int cores = std::max(1, n_core);
#ifdef _OPENMP
    omp_set_num_threads(cores);
#endif

    std::vector<EditRow> all_edits;
    std::vector<JointRow> all_joints;
    std::vector<std::string> plan_a_keys;
    int n_discard_expand = 0;
    all_edits.reserve(static_cast<size_t>(n));
    all_joints.reserve(static_cast<size_t>(n));
    plan_a_keys.reserve(static_cast<size_t>(n));

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        std::vector<EditRow> local_edits;
        std::vector<JointRow> local_joints;
        std::vector<std::string> local_keys;
        int local_discard = 0;
        std::unordered_map<std::string, AlnPack> cache;
        local_edits.reserve(64);
        local_joints.reserve(64);
        local_keys.reserve(64);

#ifdef _OPENMP
#pragma omp for schedule(dynamic, 32)
#endif
        for (int i = 0; i < n; ++i) {
            const std::string& gid = v_gene[static_cast<size_t>(i)];
            const std::string& seq = v_seq[static_cast<size_t>(i)];
            if (gid.empty() || seq.empty()) continue;
            auto pref = refs.find(gid);
            auto ppam = pams.find(gid);
            if (pref == refs.end() || ppam == pams.end() || ppam->second.empty()) {
                continue;
            }
            const std::string cache_key = gid + "\t" + seq;
            AlnPack aln;
            auto cit = cache.find(cache_key);
            if (cit != cache.end()) {
                aln = cit->second;
            } else {
                aln = nw_gapped_alignment(seq, pref->second);
                cache.emplace(cache_key, aln);
            }
            if (!aln.ok) continue;

            std::unordered_map<std::string, WinSeq> local_wins;
            local_wins.reserve(ppam->second.size());
            for (const PamRow& prow : ppam->second) {
                local_wins[prow.target_gene] = local_window_for_pam(
                    aln, prow, anchor_bp, max_expand);
            }
            std::vector<JointEv> events = classify_joint_events(
                ppam->second, local_wins, aln.ref_aln, aln.query_aln,
                check_window, anchor_bp, max_expand, min_span_bp,
                excision_tol_bp);
            for (const JointEv& ev : events) {
                JointRow jr;
                jr.sample_id = v_sample[static_cast<size_t>(i)];
                jr.index_pair_id = v_pair[static_cast<size_t>(i)];
                jr.read_id = v_read[static_cast<size_t>(i)];
                jr.gene_id = gid;
                jr.guide_i = ev.guide_i;
                jr.guide_j = ev.guide_j;
                jr.target_gene_i = ev.target_gene_i;
                jr.target_gene_j = ev.target_gene_j;
                jr.event_class = ev.event_class;
                jr.del_span = ev.del_span;
                jr.expected_span = ev.expected_span;
                local_joints.push_back(std::move(jr));
            }
            std::unordered_map<std::string, WinSeq> excision =
                resolve_excision_plan_a_alleles(events);

            bool emitted_any = false;
            for (const PamRow& prow : ppam->second) {
                WinSeq win;
                std::string allele_source = "local";
                bool intact = false;
                auto eit = excision.find(prow.target_gene);
                if (eit != excision.end()) {
                    win = eit->second;
                    allele_source = "excision";
                    intact = false;
                } else {
                    auto wit = local_wins.find(prow.target_gene);
                    if (wit != local_wins.end()) win = wit->second;
                    if (win.ok && !win.read_seq.empty()) {
                        intact = (win.read_seq == win.ref_seq);
                        allele_source = "local";
                    } else {
                        win = WinSeq();
                    }
                }
                if (!win.ok) {
                    if (eit == excision.end()) ++local_discard;
                    continue;
                }
                EditRow er;
                er.sample_id = v_sample[static_cast<size_t>(i)];
                er.index_pair_id = v_pair[static_cast<size_t>(i)];
                er.target_gene = prow.target_gene;
                er.read_seq = win.read_seq;
                er.ref_seq = win.ref_seq;
                er.intact = intact;
                er.allele_source = allele_source;
                local_edits.push_back(std::move(er));
                emitted_any = true;
            }
            if (emitted_any) {
                local_keys.push_back(
                    v_sample[static_cast<size_t>(i)] + "\t" +
                    v_pair[static_cast<size_t>(i)] + "\t" + gid);
            }
        }

#ifdef _OPENMP
#pragma omp critical(miaoseq_editcall_merge)
#endif
        {
            all_edits.insert(all_edits.end(),
                             local_edits.begin(), local_edits.end());
            all_joints.insert(all_joints.end(),
                              local_joints.begin(), local_joints.end());
            plan_a_keys.insert(plan_a_keys.end(),
                               local_keys.begin(), local_keys.end());
            n_discard_expand += local_discard;
        }
    }

    const int ne = static_cast<int>(all_edits.size());
    CharacterVector e_sample(ne), e_pair(ne), e_tg(ne), e_read(ne), e_ref(ne),
        e_src(ne);
    LogicalVector e_intact(ne);
    IntegerVector e_count(ne, 1);
    for (int i = 0; i < ne; ++i) {
        e_sample[i] = all_edits[static_cast<size_t>(i)].sample_id;
        e_pair[i] = all_edits[static_cast<size_t>(i)].index_pair_id;
        e_tg[i] = all_edits[static_cast<size_t>(i)].target_gene;
        e_read[i] = all_edits[static_cast<size_t>(i)].read_seq;
        e_ref[i] = all_edits[static_cast<size_t>(i)].ref_seq;
        e_intact[i] = all_edits[static_cast<size_t>(i)].intact;
        e_src[i] = all_edits[static_cast<size_t>(i)].allele_source;
    }
    DataFrame edit_df = DataFrame::create(
        Named("sample_id") = e_sample,
        Named("index_pair_id") = e_pair,
        Named("target_gene") = e_tg,
        Named("read_seq") = e_read,
        Named("ref_seq") = e_ref,
        Named("count") = e_count,
        Named("intact") = e_intact,
        Named("allele_source") = e_src);

    const int nj = static_cast<int>(all_joints.size());
    CharacterVector j_sample(nj), j_pair(nj), j_read(nj), j_gene(nj),
        j_gi(nj), j_gj(nj), j_tgi(nj), j_tgj(nj), j_cls(nj);
    IntegerVector j_del(nj), j_exp(nj);
    for (int i = 0; i < nj; ++i) {
        j_sample[i] = all_joints[static_cast<size_t>(i)].sample_id;
        j_pair[i] = all_joints[static_cast<size_t>(i)].index_pair_id;
        j_read[i] = all_joints[static_cast<size_t>(i)].read_id;
        j_gene[i] = all_joints[static_cast<size_t>(i)].gene_id;
        j_gi[i] = all_joints[static_cast<size_t>(i)].guide_i;
        j_gj[i] = all_joints[static_cast<size_t>(i)].guide_j;
        j_tgi[i] = all_joints[static_cast<size_t>(i)].target_gene_i;
        j_tgj[i] = all_joints[static_cast<size_t>(i)].target_gene_j;
        j_cls[i] = all_joints[static_cast<size_t>(i)].event_class;
        j_del[i] = all_joints[static_cast<size_t>(i)].del_span;
        j_exp[i] = all_joints[static_cast<size_t>(i)].expected_span;
    }
    DataFrame joint_df = DataFrame::create(
        Named("sample_id") = j_sample,
        Named("index_pair_id") = j_pair,
        Named("read_id") = j_read,
        Named("gene_id") = j_gene,
        Named("guide_i") = j_gi,
        Named("guide_j") = j_gj,
        Named("target_gene_i") = j_tgi,
        Named("target_gene_j") = j_tgj,
        Named("event_class") = j_cls,
        Named("del_span") = j_del,
        Named("expected_span") = j_exp);

    return List::create(
        Named("edit_records") = edit_df,
        Named("joint_records") = joint_df,
        Named("plan_a_keys") = wrap(plan_a_keys),
        Named("n_discard_expand") = n_discard_expand);
}
