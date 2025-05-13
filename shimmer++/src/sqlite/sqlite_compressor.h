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

    std::pair<compressor_mode, double> value_bymode() const {
        switch (mode) {
            case compressor_mode::ON_POWER:
                return {mode, power};
            case compressor_mode::ON_OPRESS:
                return {mode, outpress};
            case compressor_mode::ON_IPRESS:
                return {mode, inpress};
            case compressor_mode::ON_RATIO:
                return {mode, ratio};
            case compressor_mode::ON_MASSFLOW:
                return {mode, massflow};
            case compressor_mode::OFF_BYPASS:
            case compressor_mode::OFF_CLOSED:
                return {mode, 0.0};
        }
        throw std::invalid_argument("bad compressor mode");
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
    double      ramp_coeff;
    double      efficiency;

    std::vector<compressor_profile_sample> profile;

    bool operator<(const setting_compr_stat& other) const {
        return std::pair{i_sfrom, i_sto} < std::pair{other.i_sfrom, other.i_sto};
    }
};

namespace database {
int
load(sqlite3 *db, const optvector<int>& s_u2i,
    std::vector<setting_compr_stat>& settings);
}

} // namespace shimmer