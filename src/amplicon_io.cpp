// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(Rcpp)]]

#include <Rcpp.h>
#include "demux_internal.h"

#include <zlib.h>

#include <algorithm>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace miaoseq;

namespace {

class LineReader {
public:
    explicit LineReader(const std::string& path) : path_(path) {
        if (path.size() >= 3 &&
            (path.compare(path.size() - 3, 3, ".gz") == 0 ||
             path.compare(path.size() - 3, 3, ".GZ") == 0)) {
            gz_ = gzopen(path.c_str(), "rb");
            if (!gz_) {
                stop("Failed to open gzipped FASTQ: " + path);
            }
            gzipped_ = true;
        } else {
            plain_.open(path.c_str(), std::ios::binary);
            if (!plain_) {
                stop("Failed to open FASTQ: " + path);
            }
        }
        buf_.resize(1 << 20); // 1 MiB
    }

    ~LineReader() {
        if (gz_) gzclose(gz_);
        if (plain_.is_open()) plain_.close();
    }

    bool getline(std::string& out) {
        out.clear();
        while (true) {
            while (pos_ < size_) {
                char c = buf_[static_cast<size_t>(pos_++)];
                if (c == '\n') {
                    if (!out.empty() && out.back() == '\r') out.pop_back();
                    return true;
                }
                out.push_back(c);
            }
            if (!refill()) {
                if (!out.empty()) {
                    if (out.back() == '\r') out.pop_back();
                    return true;
                }
                return false;
            }
        }
    }

private:
    bool refill() {
        pos_ = 0;
        size_ = 0;
        if (gzipped_) {
            int n = gzread(gz_, buf_.data(), static_cast<unsigned>(buf_.size()));
            if (n <= 0) return false;
            size_ = n;
            return true;
        }
        plain_.read(buf_.data(), static_cast<std::streamsize>(buf_.size()));
        std::streamsize n = plain_.gcount();
        if (n <= 0) return false;
        size_ = static_cast<int>(n);
        return true;
    }

    std::string path_;
    bool gzipped_ = false;
    gzFile gz_ = nullptr;
    std::ifstream plain_;
    std::vector<char> buf_;
    int pos_ = 0;
    int size_ = 0;
};

inline std::string parse_fastq_id(const std::string& header) {
    size_t start = 0;
    if (!header.empty() && header[0] == '@') start = 1;
    size_t end = start;
    while (end < header.size() && header[end] != ' ' && header[end] != '\t') {
        ++end;
    }
    return header.substr(start, end - start);
}

inline std::string reverse_string(const std::string& s) {
    std::string out(s.rbegin(), s.rend());
    return out;
}

} // namespace

//' Bucket FASTQ reads by demultiplex assignments (C++)
//'
//' Streams one or more FASTQ files with a buffered reader and returns
//' per-sample `read_id` / `seq` / `qual` character vectors for target samples.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List bucket_fastq_assignments_cpp(CharacterVector fastq_files,
                                  CharacterVector read_ids,
                                  CharacterVector sample_ids,
                                  CharacterVector target_samples) {
    if (read_ids.size() != sample_ids.size()) {
        stop("read_ids and sample_ids must have equal length.");
    }

    std::unordered_map<std::string, int> sample_index;
    sample_index.reserve(static_cast<size_t>(target_samples.size()) * 2);
    for (int i = 0; i < target_samples.size(); ++i) {
        sample_index[as<std::string>(target_samples[i])] = i;
    }

    std::unordered_map<std::string, int> id2sample;
    id2sample.reserve(static_cast<size_t>(read_ids.size()) * 2);
    for (int i = 0; i < read_ids.size(); ++i) {
        const std::string sid = as<std::string>(sample_ids[i]);
        auto sit = sample_index.find(sid);
        if (sit == sample_index.end()) continue;
        id2sample[as<std::string>(read_ids[i])] = sit->second;
    }

    const int n_samples = target_samples.size();
    std::vector<std::vector<std::string>> out_ids(static_cast<size_t>(n_samples));
    std::vector<std::vector<std::string>> out_seqs(static_cast<size_t>(n_samples));
    std::vector<std::vector<std::string>> out_quals(static_cast<size_t>(n_samples));

    long long n_records = 0;
    long long n_kept = 0;

    for (int fi = 0; fi < fastq_files.size(); ++fi) {
        const std::string path = as<std::string>(fastq_files[fi]);
        LineReader reader(path);
        std::string h, s, p, q;
        while (reader.getline(h)) {
            if (!reader.getline(s) || !reader.getline(p) || !reader.getline(q)) break;
            ++n_records;
            const std::string rid = parse_fastq_id(h);
            auto it = id2sample.find(rid);
            if (it == id2sample.end()) continue;
            const int si = it->second;
            out_ids[static_cast<size_t>(si)].push_back(rid);
            out_seqs[static_cast<size_t>(si)].push_back(s);
            out_quals[static_cast<size_t>(si)].push_back(q);
            ++n_kept;
        }
    }

    List out(n_samples);
    CharacterVector names(n_samples);
    for (int i = 0; i < n_samples; ++i) {
        names[i] = target_samples[i];
        out[i] = List::create(
            Named("read_id") = wrap(out_ids[static_cast<size_t>(i)]),
            Named("seq") = wrap(out_seqs[static_cast<size_t>(i)]),
            Named("qual") = wrap(out_quals[static_cast<size_t>(i)])
        );
    }
    out.attr("names") = names;
    out.attr("n_records") = n_records;
    out.attr("n_kept") = n_kept;
    return out;
}

//' Read all FASTQ records into id/seq/qual vectors (C++)
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List read_fastq_seqs_cpp(std::string fastq) {
    LineReader reader(fastq);
    std::vector<std::string> ids;
    std::vector<std::string> seqs;
    std::vector<std::string> quals;
    ids.reserve(1024);
    seqs.reserve(1024);
    quals.reserve(1024);

    std::string h, s, p, q;
    while (reader.getline(h)) {
        if (!reader.getline(s) || !reader.getline(p) || !reader.getline(q)) break;
        ids.push_back(parse_fastq_id(h));
        seqs.push_back(s);
        quals.push_back(q);
    }
    return List::create(
        Named("read_id") = wrap(ids),
        Named("seq") = wrap(seqs),
        Named("qual") = wrap(quals)
    );
}

//' Reverse-complement sequences optionally per-read (C++)
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
CharacterVector reverse_complement_seqs_cpp(CharacterVector seqs,
                                            LogicalVector flip) {
    const int n = seqs.size();
    CharacterVector out(n);
    if (flip.size() != 1 && flip.size() != n) {
        stop("flip must be length 1 or length(seqs).");
    }
    for (int i = 0; i < n; ++i) {
        const bool do_flip = flip.size() == 1 ? flip[0] : flip[i];
        std::string s = to_upper_acgt(as<std::string>(seqs[i]));
        out[i] = do_flip ? reverse_complement(s) : s;
    }
    return out;
}

//' Reverse quality strings optionally per-read (C++)
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
CharacterVector reverse_quals_cpp(CharacterVector quals, LogicalVector flip) {
    const int n = quals.size();
    CharacterVector out(n);
    if (flip.size() != 1 && flip.size() != n) {
        stop("flip must be length 1 or length(quals).");
    }
    for (int i = 0; i < n; ++i) {
        const bool do_flip = flip.size() == 1 ? flip[0] : flip[i];
        const std::string q = as<std::string>(quals[i]);
        out[i] = do_flip ? reverse_string(q) : q;
    }
    return out;
}
