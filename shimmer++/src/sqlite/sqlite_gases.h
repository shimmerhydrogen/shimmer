/*
 * This is the SHIMMER gas network simulator.
 * Copyright (C) 2023-2024-2025 Politecnico di Torino
 * 
 * Dipartimento di Matematica "G. L. Lagrange" - DISMA
 * Dipartimento di Energia "G. Ferraris" - DENERG
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 * 
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <iostream>
#include <vector>

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