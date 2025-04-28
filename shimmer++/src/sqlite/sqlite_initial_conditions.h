#pragma once

#include <vector>
#include <algorithm>
#include <concepts>

namespace shimmer {

struct station_initial_condition {
    int     i_snum;     // Internal number of the station
    double  init_P;
    double  init_L;

    bool operator<(const station_initial_condition& other) const {
        return i_snum < other.i_snum;
    }
};

struct pipe_initial_condition {
    int     i_sfrom;
    int     i_sto;
    double  init_G;

    bool operator<(const pipe_initial_condition& other) const {
        return std::pair{i_sfrom, i_sto} < std::pair{other.i_sfrom, other.i_sto};
    }
};

}
