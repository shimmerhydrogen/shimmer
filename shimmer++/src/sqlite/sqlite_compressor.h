#pragma once

#include <iostream>
#include <vector>
#include "sqlite.hpp"


namespace shimmer {
struct compressor_profile_sample {
    compressor_mode     mode;
    double              time;
    double              power;
    double              outpress;
    double              inpress;
    double              ratio;
    double              massflow;

    bool operator<(const compressor_profile_sample& other) const {
        return time < other.time;
    }
};

struct setting_compr_stat {
    int         i_sfrom;
    int         i_sto;

    double      max_power;
    double      max_outpress;
    double      min_inpress;
    double      max_ratio;
    double      min_ratio;
    double      max_massflow;

    std::vector<compressor_profile_sample> profile;

    bool operator<(const setting_compr_stat& other) const {
        return std::pair{i_sfrom, i_sto} < std::pair{other.i_sfrom, other.i_sto};
    }
};

} // namespace shimmer