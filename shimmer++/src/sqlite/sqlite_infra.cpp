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

 /* In the database the numbering of the stations (user numbering) can be
 * non-contiguous, so we want to make it contiguous (internal numbering).
 * The mapping is stored in the two vectors `s_u2i` (user-to-internal) and
 * `s_i2u` (internal-to-user). */

#include "sqlite.hpp"
#include "errors.h"

namespace shimmer::database {

int
renumber_stations(sqlite3 *db, optvector<int>& s_u2i, std::vector<int>& s_i2u)
{
    int i_num = 0;

    /* Callback for aggregate functions COUNT() and MAX() */
    auto af_callback = [](void *res, int, char **data, char **) {
        int *ires = (int *)res;
        *ires = (data[0] != nullptr) ? atoi(data[0]) : 0;
        return 0;
    };

    int max_station_number; /* Get the maximum assigned station number */
    std::string q = "SELECT MAX(s_number) FROM stations";
    int rc = sqlite3_exec(db, q.c_str(), af_callback, &max_station_number, nullptr);
    if (rc) { goto db_problem; }

    int station_count; /* Get the number of stations in the database */
    q = "SELECT COUNT(s_number) FROM stations";
    rc = sqlite3_exec(db, q.c_str(), af_callback, &station_count, nullptr);
    if (rc) { goto db_problem; }

    if (station_count == 0) {
        return SHIMMER_SUCCESS;
    }

    /* Resize the mapping arrays */
    s_u2i.resize(max_station_number+1);
    s_i2u.resize(station_count);

    /* Compute the actual mapping */
    q = "SELECT * FROM stations ORDER BY s_number";
    sqlite3_stmt *stmt;
    rc = sqlite3_prepare_v2(db, q.c_str(), q.length(), &stmt, nullptr);
    if (rc) { goto db_problem; }

    while (sqlite3_step(stmt) != SQLITE_DONE) {
        int u_num = sqlite3_column_int(stmt, 0);
        s_u2i.at(u_num) = i_num;
        s_i2u.at(i_num) = u_num;
        i_num++;
    }
    
    sqlite3_finalize(stmt);
    return SHIMMER_SUCCESS;

db_problem:
    std::cerr << "SQL error on query '" << q << "': ";
    std::cerr << sqlite3_errmsg(db) << std::endl;
    return SHIMMER_DATABASE_PROBLEM;
}

} //namespace shimmer::database