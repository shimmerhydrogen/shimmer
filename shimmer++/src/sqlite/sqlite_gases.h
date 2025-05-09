#pragma once

#include <iostream>
#include <vector>
#include "sqlite.hpp"

#define NUM_GASES 21
namespace shimmer {

enum class gas : int {
    CH4 = 0, N2, CO2, C2H6, C3H8, iC4H10, nC4H10, iC5H12, nC5H12,
    C6H14, C7H16, C8H18, C9H20, C10H22, H2, O2, CO, H2O, H2S, He, Ar
};
struct gas_mass_fractions {
    int                             i_snum;
    std::array<double, NUM_GASES>   fractions;

    bool operator<(const gas_mass_fractions& other) const {
        return i_snum < other.i_snum;
    }
};

namespace database {

int load(sqlite3 *db, const optvector<int>& s_u2i,
    std::vector<gas_mass_fractions>& fracs);
int store(sqlite3 *db, const std::vector<int>& s_i2u,
    const std::vector<gas_mass_fractions>& fracs);

} //namespace database

} // namespace shimmer