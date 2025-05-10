#pragma once

#include <iostream>
#include <vector>

namespace shimmer {


struct setting_pipe {
    int         i_sfrom;
    int         i_sto;

    double      length;
    double      diameter;
    double      roughness;

    bool operator<(const setting_pipe& other) const {
        return std::pair{i_sfrom, i_sto} < std::pair{other.i_sfrom, other.i_sto};
    }
};

namespace database {

int load(sqlite3 *db, const optvector<int>& s_u2i,
    std::vector<setting_pipe>& settings);
int store(sqlite3 *db, const std::vector<int>& s_i2u,
    const std::vector<setting_pipe>& settings);

} //namespace database


} // namespace shimmer