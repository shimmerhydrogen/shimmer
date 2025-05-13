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
#include <concepts>

namespace shimmer {

enum class station_type : int {
    ENTRY_P_REG = 1,        /* ReMi w/o backflow */
    ENTRY_L_REG = 2,        /* Injection w/ pressure control */
    EXIT_L_REG  = 3,        /* Consumption point w/o pressure control */
    JUNCTION    = 4,        /* Junction */
    PRIVATE_INLET  = 10,    /* Inlet, internal use only */
    PRIVATE_OUTLET = 11     /* Outlet, internal use only */
};

enum class pipe_type : int {
    PIPE        = 0,    /* Plain pipe */
    COMPR_STAT  = 1,    /* Compressor station */
    RED_STAT    = 2,    /* Reduction station */
    VALVE       = 3     /* Valve */
};

enum class compressor_mode : int {
    ON_POWER    = 0,    /* Compressor on, Control mode power driver */
    ON_OPRESS   = 1,    /* Compressor on, Control mode outlet pressure */
    ON_IPRESS   = 2,    /* Compressor on, Control mode inlet pressure */
    ON_RATIO    = 3,    /* Compressor on, Control mode compression ratio */
    ON_MASSFLOW = 4,    /* Compressor on, Control mode mass flow */
    OFF_BYPASS  = 10,   /* Compressor off, bypass */
    OFF_CLOSED  = 11    /* Compressor off, closed */
};

enum class reduction_mode : int {
    ON_OPRESS   = 0,    /* Reduction station on, Control mode outlet pressure */
    ON_IPRESS   = 1,    /* Reduction station on, Control mode inlet pressure */
    ON_RATIO    = 2,    /* Reduction station on, Control mode compression ratio */
    ON_MASSFLOW = 3,    /* Reduction station on, Control mode mass flow */
    OFF_BYPASS  = 10,   /* Reduction station off, bypass */
    OFF_CLOSED  = 11    /* Reduction station off, closed */
};

enum class gas : int {
    CH4 = 0, N2, CO2, C2H6, C3H8, iC4H10, nC4H10, iC5H12, nC5H12,
    C6H14, C7H16, C8H18, C9H20, C10H22, H2, O2, CO, H2O, H2S, He, Ar
};

} // namespace shimmer

template<typename T>
    requires std::is_enum_v<T>
constexpr auto operator+(T e) {
    return std::underlying_type_t<T>(e);
}
