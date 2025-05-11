#pragma once

#include <vector>
#include <algorithm>

namespace shimmer {

struct station_initial_condition {
    int     i_snum;
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

namespace database {

int load(sqlite3 *db, const optvector<int>& s_u2i,
    std::vector<station_initial_condition>& sics);
int store(sqlite3 *db, const std::vector<int>& s_i2u,
    const std::vector<station_initial_condition>& sics);
int load(sqlite3 *db, const optvector<int>& s_u2i,
    std::vector<pipe_initial_condition>& pics);
int store(sqlite3 *db, const std::vector<int>& s_i2u,
    const std::vector<pipe_initial_condition>& pics);

} //namespace database

}
