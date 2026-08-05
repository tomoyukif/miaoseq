#ifndef MIAOSEQ_FASTQ_IO_H
#define MIAOSEQ_FASTQ_IO_H

#include <zlib.h>

#include <fstream>
#include <string>
#include <utility>
#include <vector>

#include <Rcpp.h>

namespace miaoseq {

inline std::string parse_fastq_id(const std::string& header) {
    size_t start = 0;
    if (!header.empty() && header[0] == '@') start = 1;
    size_t end = start;
    while (end < header.size() && header[end] != ' ' && header[end] != '\t') {
        ++end;
    }
    return header.substr(start, end - start);
}

inline std::string safe_filename(const std::string& x) {
    std::string out = x;
    for (char& c : out) {
        const bool ok =
            (c >= 'A' && c <= 'Z') ||
            (c >= 'a' && c <= 'z') ||
            (c >= '0' && c <= '9') ||
            c == '.' || c == '_' || c == '-';
        if (!ok) c = '_';
    }
    return out;
}

class LineReader {
public:
    explicit LineReader(const std::string& path) : path_(path) {
        if (path.size() >= 3 &&
            (path.compare(path.size() - 3, 3, ".gz") == 0 ||
             path.compare(path.size() - 3, 3, ".GZ") == 0)) {
            gz_ = gzopen(path.c_str(), "rb");
            if (!gz_) {
                Rcpp::stop("Failed to open gzipped FASTQ: " + path);
            }
            gzipped_ = true;
        } else {
            plain_.open(path.c_str(), std::ios::binary);
            if (!plain_) {
                Rcpp::stop("Failed to open FASTQ: " + path);
            }
        }
        buf_.resize(1 << 20);
    }

    ~LineReader() {
        if (gz_) gzclose(gz_);
        if (plain_.is_open()) plain_.close();
    }

    LineReader(const LineReader&) = delete;
    LineReader& operator=(const LineReader&) = delete;

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

class FastqWriter {
public:
    FastqWriter() = default;

    void open(const std::string& path, bool compress) {
        close();
        path_ = path;
        compress_ = compress;
        buf_.clear();
        buf_.reserve(1 << 20);
        if (compress_) {
            gz_ = gzopen(path.c_str(), "wb6");
            if (!gz_) Rcpp::stop("Failed to open for writing: " + path);
        } else {
            plain_.open(path.c_str(), std::ios::binary);
            if (!plain_) Rcpp::stop("Failed to open for writing: " + path);
        }
        open_ = true;
    }

    ~FastqWriter() { close(); }

    FastqWriter(const FastqWriter&) = delete;
    FastqWriter& operator=(const FastqWriter&) = delete;

    FastqWriter(FastqWriter&& other) noexcept {
        *this = std::move(other);
    }

    FastqWriter& operator=(FastqWriter&& other) noexcept {
        if (this == &other) return *this;
        close();
        path_ = std::move(other.path_);
        compress_ = other.compress_;
        open_ = other.open_;
        gz_ = other.gz_;
        plain_ = std::move(other.plain_);
        buf_ = std::move(other.buf_);
        other.gz_ = nullptr;
        other.open_ = false;
        return *this;
    }

    void write_record(const std::string& header,
                      const std::string& seq,
                      const std::string& plus,
                      const std::string& qual) {
        if (!open_) return;
        if (!header.empty() && header[0] == '@') {
            buf_.append(header);
        } else {
            buf_.push_back('@');
            buf_.append(header);
        }
        buf_.push_back('\n');
        buf_.append(seq);
        buf_.push_back('\n');
        if (plus.empty()) {
            buf_.push_back('+');
        } else {
            buf_.append(plus);
        }
        buf_.push_back('\n');
        buf_.append(qual);
        buf_.push_back('\n');
        if (buf_.size() >= (1u << 20)) flush();
    }

    void flush() {
        if (!open_ || buf_.empty()) return;
        if (compress_) {
            if (gzwrite(gz_, buf_.data(), static_cast<unsigned>(buf_.size())) == 0) {
                Rcpp::stop("Failed writing gzipped FASTQ: " + path_);
            }
        } else {
            plain_.write(buf_.data(), static_cast<std::streamsize>(buf_.size()));
            if (!plain_) Rcpp::stop("Failed writing FASTQ: " + path_);
        }
        buf_.clear();
    }

    void close() {
        if (!open_) return;
        flush();
        if (gz_) {
            gzclose(gz_);
            gz_ = nullptr;
        }
        if (plain_.is_open()) plain_.close();
        open_ = false;
    }

private:
    std::string path_;
    bool compress_ = false;
    bool open_ = false;
    gzFile gz_ = nullptr;
    std::ofstream plain_;
    std::string buf_;
};

class TsvWriter {
public:
    TsvWriter() = default;

    void open(const std::string& path, bool append = false) {
        close();
        path_ = path;
        plain_.open(
            path.c_str(),
            append ? (std::ios::binary | std::ios::app)
                   : (std::ios::binary | std::ios::trunc));
        if (!plain_) Rcpp::stop("Failed to open TSV for writing: " + path);
        buf_.clear();
        buf_.reserve(1 << 20);
        open_ = true;
    }

    ~TsvWriter() { close(); }

    TsvWriter(const TsvWriter&) = delete;
    TsvWriter& operator=(const TsvWriter&) = delete;

    void write_line(const std::string& line) {
        if (!open_) return;
        buf_.append(line);
        if (line.empty() || line.back() != '\n') buf_.push_back('\n');
        if (buf_.size() >= (1u << 20)) flush();
    }

    void flush() {
        if (!open_ || buf_.empty()) return;
        plain_.write(buf_.data(), static_cast<std::streamsize>(buf_.size()));
        if (!plain_) Rcpp::stop("Failed writing TSV: " + path_);
        buf_.clear();
    }

    void close() {
        if (!open_) return;
        flush();
        if (plain_.is_open()) plain_.close();
        open_ = false;
    }

private:
    std::string path_;
    bool open_ = false;
    std::ofstream plain_;
    std::string buf_;
};

} // namespace miaoseq

#endif
