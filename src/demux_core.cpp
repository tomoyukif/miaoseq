// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(Rcpp)]]

#include <Rcpp.h>
#include "demux_engine.h"
#include "demux_internal.h"

#include <algorithm>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace miaoseq;

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
    DemuxEngine engine(
        f_suffix, r_suffix,
        as<std::vector<std::string>>(f_barcodes),
        as<std::vector<std::string>>(f_barcode_names),
        as<std::vector<std::string>>(r_barcodes),
        as<std::vector<std::string>>(r_barcode_names),
        as<std::vector<std::string>>(pair_f_names),
        as<std::vector<std::string>>(pair_r_names),
        as<std::vector<std::string>>(pair_ids),
        as<std::vector<std::string>>(sample_ids),
        end_window, max_anchor_edit, max_barcode_edit,
        allow_revcomp, allow_single_end);

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
        DemuxOut hit = engine.classify(v_seqs[static_cast<size_t>(i)], &local_edlib);
        status[static_cast<size_t>(i)] = hit.status;
        reason[static_cast<size_t>(i)] = hit.reason;
        out_pair[static_cast<size_t>(i)] = hit.index_pair_id;
        out_sample[static_cast<size_t>(i)] = hit.sample_id;
        out_f[static_cast<size_t>(i)] = hit.f_index_id;
        out_r[static_cast<size_t>(i)] = hit.r_index_id;
        be_f[static_cast<size_t>(i)] = hit.barcode_edit_f;
        be_r[static_cast<size_t>(i)] = hit.barcode_edit_r;
        match_class[static_cast<size_t>(i)] = hit.match_class;
        assign_mode[static_cast<size_t>(i)] = hit.assign_mode;
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
