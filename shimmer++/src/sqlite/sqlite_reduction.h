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
struct reduction_profile_sample {
    reduction_mode      mode;
    double              time;
    double              power;
    double              outpress;
    double              inpress;
    double              ratio;
    double              massflow;

    bool operator<(const reduction_profile_sample& other) const {
        return time < other.time;
    }

    std::pair<reduction_mode, double> value_bymode() const {
        switch (mode) {
            case reduction_mode::ON_OPRESS:
                return {mode, outpress};
            case reduction_mode::ON_IPRESS:
                return {mode, inpress};
            case reduction_mode::ON_RATIO:
                return {mode, ratio};
            case reduction_mode::ON_MASSFLOW:
                return {mode, massflow};
            case reduction_mode::OFF_BYPASS:
            case reduction_mode::OFF_CLOSED:
                return {mode, 0.0};
        }
        throw std::invalid_argument("bad compressor mode");
    }
};

struct setting_red_stat {
    int         i_sfrom;
    int         i_sto;

    double      max_outpress;
    double      min_inpress;
    double      max_ratio;
    double      min_ratio;;
    double      max_massflow;

    std::vector<reduction_profile_sample> profile;

    bool operator<(const setting_red_stat& other) const {
        return std::pair{i_sfrom, i_sto} < std::pair{other.i_sfrom, other.i_sto};
    }
};

} // namespace shimmer