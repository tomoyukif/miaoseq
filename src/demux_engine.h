#ifndef MIAOSEQ_DEMUX_ENGINE_H
#define MIAOSEQ_DEMUX_ENGINE_H

#include "demux_internal.h"

#include <string>
#include <unordered_map>
#include <vector>

namespace miaoseq {

struct DemuxOut {
    int status = 0; // 1 = assigned
    std::string reason = "no_suffix";
    std::string index_pair_id;
    std::string sample_id;
    std::string f_index_id;
    std::string r_index_id;
    int barcode_edit_f = -1;
    int barcode_edit_r = -1;
    std::string match_class;
    std::string assign_mode;
};

class DemuxEngine {
public:
    DemuxEngine(const std::string& f_suffix,
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
                bool allow_single_end);

    DemuxOut classify(const std::string& seq, long long* n_edlib) const;

private:
    std::string f_suffix_;
    std::string r_suffix_;
    MutantDict f_dict_;
    MutantDict r_dict_;
    int f_bc_len_ = 0;
    int r_bc_len_ = 0;
    int end_window_ = 120;
    int max_anchor_edit_ = 10;
    bool allow_revcomp_ = true;
    bool allow_single_end_ = false;
    std::vector<std::string> pair_ids_;
    std::vector<std::string> sample_ids_;
    std::unordered_map<std::string, int> pair_map_;
    std::unordered_map<std::string, int> f_name_to_row_;
    std::unordered_map<std::string, int> r_name_to_row_;
};

} // namespace miaoseq

#endif
