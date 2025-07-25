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

#include <optional>
#include <cassert>
#include <sqlite3.h>
#include "errors.h"
#include "sqlite.hpp"

namespace shimmer {

namespace database {

int load(sqlite3 *db, const optvector<int>& s_u2i,
    std::vector<gas_molar_fractions>& fracs)
{
    sqlite3_stmt *stmt = nullptr;

    std::string q = "SELECT * FROM gas_molar_fractions";
    int rc = sqlite3_prepare_v2(db, q.c_str(), q.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << q << "': ";
        std::cerr << sqlite3_errmsg(db) << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    }

    if ( sqlite3_column_count(stmt) != 1+NUM_GASES ) {
        std::cerr << "ERROR: invalid table column number while importing ";
        std::cerr << "gas molar fractions" << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    }

    /* Import gas molar fractions for all the stations */
    while (sqlite3_step(stmt) == SQLITE_ROW) {
        gas_molar_fractions frac;
        int u_snum = sqlite3_column_int(stmt, /*column*/ 0);
        auto i_snum_opt = convert_u2i(s_u2i, u_snum);
        if (not i_snum_opt) {
            std::cerr << "ERROR: while importing gas molar fractions for station ";
            std::cerr << u_snum << ": s_u2i got an invalid station number.";
            std::cerr << std::endl;
            return SHIMMER_DATABASE_PROBLEM;
        }

        frac.i_snum = i_snum_opt.value();
        assert(frac.fractions.size() == NUM_GASES);
        for (int i = 0; i < NUM_GASES; i++) {
            frac.fractions[i] = sqlite3_column_double(stmt, 1+i);
        }
        
        fracs.push_back( std::move(frac) );
    }
    rc = sqlite3_finalize(stmt);

    std::sort(fracs.begin(), fracs.end());
    return SHIMMER_SUCCESS;
}

int store(sqlite3 *db, const std::vector<int>& s_i2u,
    const std::vector<gas_molar_fractions>& fracs)
{
    char *errmsg;
    sqlite3_stmt *stmt = nullptr;

    std::string q = "INSERT INTO gas_molar_fractions VALUES "
        "(?, ?, ?, ?, ?, ?, ?, ?, "
        "    ?, ?, ?, ?, ?, ?, ?, "
        "    ?, ?, ?, ?, ?, ?, ?)";
    
    int rc = sqlite3_prepare_v2(db, q.c_str(), q.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << q << "': ";
        std::cerr << sqlite3_errmsg(db) << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    }

    rc = sqlite3_exec(db, "BEGIN TRANSACTION", nullptr, nullptr, &errmsg);

    for (int i = 0; i < fracs.size(); i++) {
        int snum = convert_i2u(s_i2u, fracs[i].i_snum);
        rc = sqlite3_bind_int(stmt, 1, snum);
        for (size_t ig = 0; ig < NUM_GASES; ig++) {
            rc = sqlite3_bind_double(stmt, ig+2, fracs[i].fractions[ig]);
        }
        rc = sqlite3_step(stmt);
        rc = sqlite3_clear_bindings(stmt);
        rc = sqlite3_reset(stmt);
    }

    rc = sqlite3_exec(db, "COMMIT", nullptr, nullptr, &errmsg);

    return SHIMMER_SUCCESS;
}

} //namespace database



} //namespace shimmer