#pragma once

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
        return (i_sfrom < other.i_sfrom) or
            (i_sfrom == other.i_sfrom and i_sto < other.i_sto);
    }
};

}

template<typename T>
    requires std::is_enum_v<T>
constexpr auto operator+(T e) {
    return std::underlying_type_t<T>(e);
}