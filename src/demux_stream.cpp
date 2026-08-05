// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(Rcpp)]]

#include <Rcpp.h>
#include "demux_engine.h"
#include "fastq_io.h"

#include <algorithm>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace miaoseq;

namespace {

struct SampleStats {
    long long n_reads = 0;
    long long n_complete = 0;
    long long n_fuzzy = 0;
    long long n_single_end = 0;
};

struct FastqRec {
    std::string header;
    std::string seq;
    std::string plus;
    std::string qual;
    std::string id;
};

std::string escape_tsv(const std::string& s) {
    // Fields are barcode IDs / reasons; no tabs expected. Keep simple.
    return s;
}

void write_assignment_line(TsvWriter& w,
                           const std::string& rid,
                           const DemuxOut& hit,
                           const std::string& source_file) {
    std::ostringstream oss;
    oss << escape_tsv(rid) << '\t'
        << escape_tsv(hit.index_pair_id) << '\t'
        << escape_tsv(hit.f_index_id) << '\t'
        << escape_tsv(hit.r_index_id) << '\t'
        << escape_tsv(hit.sample_id) << '\t'
        << escape_tsv(source_file) << '\t'
        << hit.barcode_edit_f << '\t'
        << hit.barcode_edit_r << '\t'
        << escape_tsv(hit.match_class) << '\t'
        << escape_tsv(hit.assign_mode);
    w.write_line(oss.str());
}

void write_unassigned_line(TsvWriter& w,
                           const std::string& rid,
                           const std::string& reason,
                           const std::string& source_file) {
    std::ostringstream oss;
    oss << escape_tsv(rid) << '\t'
        << escape_tsv(reason) << '\t'
        << escape_tsv(source_file);
    w.write_line(oss.str());
}

} // namespace

//' Stream-demultiplex FASTQ files; write TSV (+ optional per-index-pair FASTQ)
//'
//' Reads FASTQ in chunks, demultiplexes with OpenMP, and streams
//' `assignments.tsv` / `unassigned.tsv`. When `split_reads` is TRUE, also
//' writes one FASTQ per `index_pair_id` in the same pass.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List demux_fastq_stream_cpp(CharacterVector fastq_files,
                            std::string demult_dir,
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
                            int n_core = 1,
                            int chunk_size = 20000,
                            bool split_reads = false,
                            bool compress = true,
                            bool include_unassigned_fastq = false) {
    if (chunk_size < 1) chunk_size = 1;
    const int n_pairs = pair_ids.size();
    if (n_pairs < 1) stop("No index pairs provided.");
    if (sample_ids.size() != n_pairs ||
        pair_f_names.size() != n_pairs ||
        pair_r_names.size() != n_pairs) {
        stop("pair/sample vectors must have equal length.");
    }

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

    std::vector<std::string> v_sample_ids = as<std::vector<std::string>>(sample_ids);
    std::vector<std::string> v_pair_ids = as<std::vector<std::string>>(pair_ids);
    std::unordered_map<std::string, int> pair_to_idx;
    pair_to_idx.reserve(static_cast<size_t>(n_pairs) * 2);
    std::vector<SampleStats> stats(static_cast<size_t>(n_pairs));
    for (int i = 0; i < n_pairs; ++i) {
        pair_to_idx.emplace(v_pair_ids[static_cast<size_t>(i)], i);
    }

    const std::string assign_path = demult_dir + "/assignments.tsv";
    const std::string unassign_path = demult_dir + "/unassigned.tsv";
    TsvWriter assign_w;
    TsvWriter unassign_w;
    assign_w.open(assign_path, false);
    unassign_w.open(unassign_path, false);
    assign_w.write_line(
        "read_id\tindex_pair_id\tf_index_id\tr_index_id\tsample_id\t"
        "source_file\tbarcode_edit_f\tbarcode_edit_r\tmatch_class\tassign_mode");
    unassign_w.write_line("read_id\treason\tsource_file");

    std::string by_sample_dir = demult_dir + "/by_sample";
    const std::string fq_ext = compress ? ".fq.gz" : ".fq";
    std::vector<FastqWriter> sample_writers;
    std::vector<char> sample_opened;
    FastqWriter unassigned_writer;
    bool unassigned_opened = false;
    std::unordered_map<std::string, int> reason_counts;

    if (split_reads) {
        sample_writers.resize(static_cast<size_t>(n_pairs));
        sample_opened.assign(static_cast<size_t>(n_pairs), 0);
    }

    int cores = std::max(1, n_core);
#ifdef _OPENMP
    omp_set_num_threads(cores);
#endif

    long long n_edlib = 0;
    long long n_records = 0;
    long long n_assigned = 0;
    long long n_unassigned = 0;
    long long n_incomplete = 0;

    for (int fi = 0; fi < fastq_files.size(); ++fi) {
        const std::string fq = as<std::string>(fastq_files[fi]);
        LineReader reader(fq);
        std::vector<FastqRec> chunk;
        chunk.reserve(static_cast<size_t>(chunk_size));

        auto flush_chunk = [&]() {
            const int n = static_cast<int>(chunk.size());
            if (n < 1) return;
            std::vector<DemuxOut> hits(static_cast<size_t>(n));
            long long local_sum = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 64) reduction(+:local_sum)
#endif
            for (int i = 0; i < n; ++i) {
                long long local = 0;
                hits[static_cast<size_t>(i)] =
                    engine.classify(chunk[static_cast<size_t>(i)].seq, &local);
                local_sum += local;
            }
            n_edlib += local_sum;

            for (int i = 0; i < n; ++i) {
                const FastqRec& rec = chunk[static_cast<size_t>(i)];
                const DemuxOut& hit = hits[static_cast<size_t>(i)];
                if (hit.status == 1) {
                    ++n_assigned;
                    write_assignment_line(assign_w, rec.id, hit, fq);
                    auto sit = pair_to_idx.find(hit.index_pair_id);
                    if (sit != pair_to_idx.end()) {
                        SampleStats& st = stats[static_cast<size_t>(sit->second)];
                        ++st.n_reads;
                        if (hit.match_class == "complete_match") ++st.n_complete;
                        if (hit.match_class == "fuzzy_match") ++st.n_fuzzy;
                        if (hit.assign_mode == "single_f" ||
                            hit.assign_mode == "single_r") {
                            ++st.n_single_end;
                        }
                        if (split_reads) {
                            const size_t si = static_cast<size_t>(sit->second);
                            if (!sample_opened[si]) {
                                const std::string path =
                                    by_sample_dir + "/" +
                                    safe_filename(v_pair_ids[si]) + fq_ext;
                                sample_writers[si].open(path, compress);
                                sample_opened[si] = 1;
                            }
                            sample_writers[si].write_record(
                                rec.header, rec.seq, rec.plus, rec.qual);
                        }
                    }
                } else {
                    ++n_unassigned;
                    write_unassigned_line(unassign_w, rec.id, hit.reason, fq);
                    reason_counts[hit.reason]++;
                    if (split_reads && include_unassigned_fastq) {
                        if (!unassigned_opened) {
                            unassigned_writer.open(
                                by_sample_dir + "/unassigned" + fq_ext, compress);
                            unassigned_opened = true;
                        }
                        unassigned_writer.write_record(
                            rec.header, rec.seq, rec.plus, rec.qual);
                    }
                }
            }
            chunk.clear();
        };

        std::string h, s, p, q;
        while (reader.getline(h)) {
            if (!reader.getline(s) || !reader.getline(p) || !reader.getline(q)) {
                ++n_incomplete;
                break;
            }
            FastqRec rec;
            rec.header = h;
            rec.seq = s;
            rec.plus = p;
            rec.qual = q;
            rec.id = parse_fastq_id(h);
            chunk.push_back(std::move(rec));
            ++n_records;
            if (static_cast<int>(chunk.size()) >= chunk_size) flush_chunk();
        }
        flush_chunk();
    }

    assign_w.close();
    unassign_w.close();
    if (split_reads) {
        for (size_t i = 0; i < sample_writers.size(); ++i) {
            if (sample_opened[i]) sample_writers[i].close();
        }
        if (unassigned_opened) unassigned_writer.close();
    }

    CharacterVector summary_sample(n_pairs);
    CharacterVector summary_pair(n_pairs);
    IntegerVector summary_n(n_pairs);
    IntegerVector summary_complete(n_pairs);
    IntegerVector summary_fuzzy(n_pairs);
    IntegerVector summary_single(n_pairs);
    for (int i = 0; i < n_pairs; ++i) {
        summary_sample[i] = v_sample_ids[static_cast<size_t>(i)];
        summary_pair[i] = v_pair_ids[static_cast<size_t>(i)];
        summary_n[i] = static_cast<int>(stats[static_cast<size_t>(i)].n_reads);
        summary_complete[i] =
            static_cast<int>(stats[static_cast<size_t>(i)].n_complete);
        summary_fuzzy[i] =
            static_cast<int>(stats[static_cast<size_t>(i)].n_fuzzy);
        summary_single[i] =
            static_cast<int>(stats[static_cast<size_t>(i)].n_single_end);
    }

    const int n_reasons = static_cast<int>(reason_counts.size());
    CharacterVector reason_names(n_reasons);
    IntegerVector reason_n(n_reasons);
    int ri = 0;
    for (const auto& kv : reason_counts) {
        reason_names[ri] = kv.first;
        reason_n[ri] = static_cast<int>(kv.second);
        ++ri;
    }

    return List::create(
        Named("n_records") = n_records,
        Named("n_assigned") = n_assigned,
        Named("n_unassigned") = n_unassigned,
        Named("n_incomplete") = n_incomplete,
        Named("n_edlib") = static_cast<double>(n_edlib),
        Named("assignments_tsv") = assign_path,
        Named("unassigned_tsv") = unassign_path,
        Named("summary_sample_id") = summary_sample,
        Named("summary_index_pair_id") = summary_pair,
        Named("summary_n_reads") = summary_n,
        Named("summary_n_complete") = summary_complete,
        Named("summary_n_fuzzy") = summary_fuzzy,
        Named("summary_n_single_end") = summary_single,
        Named("reason_names") = reason_names,
        Named("reason_n") = reason_n
    );
}

//' Split FASTQ by assignments TSV (C++ streaming)
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List split_fastq_by_assignment_cpp(CharacterVector fastq_files,
                                   std::string assignments_tsv,
                                   std::string out_dir,
                                   bool compress = true,
                                   bool include_unassigned = false,
                                   std::string unassigned_tsv = "") {
    std::unordered_map<std::string, std::string> id2sample;
    bool keyed_with_source = false;
    {
        LineReader ar(assignments_tsv);
        std::string line;
        bool first = true;
        int source_col = -1;
        while (ar.getline(line)) {
            if (first) {
                first = false;
                if (line.find("read_id") != std::string::npos) {
                    std::vector<std::string> cols;
                    size_t a = 0;
                    while (true) {
                        size_t b = line.find('\t', a);
                        cols.push_back(b == std::string::npos ? line.substr(a)
                                                             : line.substr(a, b - a));
                        if (b == std::string::npos) break;
                        a = b + 1;
                    }
                    for (int ci = 0; ci < static_cast<int>(cols.size()); ++ci) {
                        if (cols[static_cast<size_t>(ci)] == "source_file") {
                            source_col = ci;
                            keyed_with_source = true;
                            break;
                        }
                    }
                    continue;
                }
            }
            if (line.empty()) continue;
            std::vector<std::string> cols;
            size_t a = 0;
            while (true) {
                size_t b = line.find('\t', a);
                cols.push_back(b == std::string::npos ? line.substr(a)
                                                     : line.substr(a, b - a));
                if (b == std::string::npos) break;
                a = b + 1;
            }
            if (cols.size() < 2 || cols[0].empty() || cols[1].empty()) continue;
            if (keyed_with_source && source_col >= 0 &&
                source_col < static_cast<int>(cols.size()) &&
                !cols[static_cast<size_t>(source_col)].empty()) {
                id2sample[cols[0] + "\t" + cols[static_cast<size_t>(source_col)]] =
                    cols[1];
            } else {
                id2sample[cols[0]] = cols[1];
            }
        }
    }

    std::unordered_map<std::string, bool> un_ids;
    if (include_unassigned && !unassigned_tsv.empty()) {
        LineReader ur(unassigned_tsv);
        std::string line;
        bool first = true;
        while (ur.getline(line)) {
            if (first) {
                first = false;
                if (line.find("read_id") != std::string::npos) continue;
            }
            if (line.empty()) continue;
            size_t p = line.find('\t');
            const std::string rid =
                (p == std::string::npos) ? line : line.substr(0, p);
            if (!rid.empty()) un_ids[rid] = true;
        }
    }

    // Collect unique samples and open writers
    std::vector<std::string> samples;
    std::unordered_map<std::string, int> sample_idx;
    for (const auto& kv : id2sample) {
        if (sample_idx.find(kv.second) == sample_idx.end()) {
            sample_idx[kv.second] = static_cast<int>(samples.size());
            samples.push_back(kv.second);
        }
    }
    std::sort(samples.begin(), samples.end());
    sample_idx.clear();
    for (size_t i = 0; i < samples.size(); ++i) {
        sample_idx[samples[i]] = static_cast<int>(i);
    }

    const std::string ext = compress ? ".fq.gz" : ".fq";
    std::vector<FastqWriter> writers(samples.size());
    std::vector<long long> counts(samples.size(), 0);
    for (size_t i = 0; i < samples.size(); ++i) {
        writers[i].open(out_dir + "/" + safe_filename(samples[i]) + ext, compress);
    }
    FastqWriter un_writer;
    long long n_un = 0;
    if (include_unassigned) {
        un_writer.open(out_dir + "/unassigned" + ext, compress);
    }

    long long n_records = 0;
    long long n_written = 0;
    long long n_incomplete = 0;
    for (int fi = 0; fi < fastq_files.size(); ++fi) {
        const std::string fq = as<std::string>(fastq_files[fi]);
        LineReader reader(fq);
        std::string h, s, p, q;
        while (reader.getline(h)) {
            if (!reader.getline(s) || !reader.getline(p) || !reader.getline(q)) {
                ++n_incomplete;
                break;
            }
            ++n_records;
            const std::string rid = parse_fastq_id(h);
            auto it = keyed_with_source ? id2sample.find(rid + "\t" + fq)
                                        : id2sample.find(rid);
            if (it == id2sample.end() && keyed_with_source) {
                it = id2sample.find(rid);
            }
            if (it != id2sample.end()) {
                auto sit = sample_idx.find(it->second);
                if (sit != sample_idx.end()) {
                    writers[static_cast<size_t>(sit->second)]
                        .write_record(h, s, p, q);
                    counts[static_cast<size_t>(sit->second)]++;
                    ++n_written;
                }
            } else if (include_unassigned && un_ids.count(rid)) {
                un_writer.write_record(h, s, p, q);
                ++n_un;
            }
        }
    }

    for (auto& w : writers) w.close();
    if (include_unassigned) un_writer.close();

    CharacterVector out_samples(samples.size());
    IntegerVector out_counts(samples.size());
    for (size_t i = 0; i < samples.size(); ++i) {
        out_samples[static_cast<int>(i)] = samples[i];
        out_counts[static_cast<int>(i)] = static_cast<int>(counts[i]);
    }

    return List::create(
        Named("index_pair_id") = out_samples,
        Named("n_reads") = out_counts,
        Named("n_records") = n_records,
        Named("n_written") = n_written,
        Named("n_unassigned_written") = n_un,
        Named("n_incomplete") = n_incomplete
    );
}
