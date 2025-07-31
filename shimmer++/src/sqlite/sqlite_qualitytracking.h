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


struct setting_pipe_QT {
    int         i_sfrom;
    int         i_sto;

    int         nsegs;

    bool operator<(const setting_pipe_QT& other) const {
        return std::pair{i_sfrom, i_sto} < std::pair{other.i_sfrom, other.i_sto};
    }
};

namespace database {

int load(sqlite3 *db, const optvector<int>& s_u2i,
    std::vector<setting_pipe_QT>& settings);
int store(sqlite3 *db, const std::vector<int>& s_i2u,
    const std::vector<setting_pipe_QT>& settings);

} //namespace database


} // namespace shimmer