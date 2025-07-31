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
#include <sqlite3.h>
#include "errors.h"
#include "sqlite.hpp"

/* Database I/O functions for pipe type 'Plain pipe' */

namespace shimmer {

namespace pipe {

/* Define indices of columns in parameters table */
enum class params_col : int {
    p_name      = 0,
    s_from      = 1,
    s_to        = 2,
    diameter    = 3,
    length      = 4,
    roughness   = 5,
    ref_nsegs   = 6,
};

} // namespace pipe

namespace database {

int load(sqlite3 *db, const optvector<int>& s_u2i, 
    std::vector<setting_pipe>& settings) {
    using namespace pipe;

    sqlite3_stmt *stmt = nullptr;

    std::string qlim = "SELECT pipe_parameters.* "
	    "FROM pipe_parameters INNER JOIN pipelines "
	    "ON pipe_parameters.p_name = pipelines.p_name "
		"   AND pipe_parameters.s_from = pipelines.s_from "
		"   AND pipe_parameters.s_to = pipelines.s_to "
	    "WHERE pipelines.p_type = ?";

    int rc = sqlite3_prepare_v2(db, qlim.c_str(), qlim.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << qlim << "': ";
        std::cerr << sqlite3_errmsg(db) << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    }

    rc = sqlite3_bind_int(stmt, 1, +pipe_type::PIPE);

    int table_cols = sqlite3_column_count(stmt);

    /* Import limits for all the stations */
    while (sqlite3_step(stmt) == SQLITE_ROW) {
        setting_pipe setting;
        int u_from = sqlite3_column_int(stmt, +params_col::s_from);
        int u_to = sqlite3_column_int(stmt, +params_col::s_to);
        
        auto i_sfrom_opt = convert_u2i(s_u2i, u_from);
        auto i_sto_opt = convert_u2i(s_u2i, u_to);
        if (not i_sfrom_opt or not i_sto_opt) {
            std::cerr << "s_u2i: invalid station numbers. Inconsistent data in DB?" << std::endl;
            return SHIMMER_DATABASE_PROBLEM;
        }

        setting.i_sfrom = i_sfrom_opt.value();
        setting.i_sto = i_sto_opt.value();
        setting.diameter = sqlite3_column_double(stmt, +params_col::diameter);
        setting.length = sqlite3_column_double(stmt, +params_col::length);
        setting.roughness = sqlite3_column_double(stmt, +params_col::roughness);
        setting.ref_nsegs = 0;
        if (table_cols > 6) {
            setting.ref_nsegs = sqlite3_column_int(stmt, +params_col::ref_nsegs);
        }
        settings.push_back( std::move(setting) );
    }
    rc = sqlite3_finalize(stmt);

    std::sort(settings.begin(), settings.end());
    return SHIMMER_SUCCESS;
}

int store(sqlite3 *db, const std::vector<int>& s_i2u,
    const std::vector<setting_pipe>& settings)
{
    sqlite3_stmt *stmt = nullptr;
    using namespace pipe;

    std::string q =
        "INSERT INTO pipe_parameters "
        "VALUES (?, ?, ?, ?, ?, ?, ?)";
        
    int rc = sqlite3_prepare_v2(db, q.c_str(), q.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << q << "': ";
        std::cerr << sqlite3_errmsg(db) << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    }

    rc = sqlite3_exec(db, "BEGIN TRANSACTION", nullptr, nullptr, nullptr);

    for (int i = 0; i < settings.size(); i++) {
        auto from = convert_i2u(s_i2u, settings[i].i_sfrom);
        auto to = convert_i2u(s_i2u, settings[i].i_sto);
        std::string pipe_name = "pipe-" + std::to_string(from) + "-" + std::to_string(to);

        rc = sqlite3_bind_text(stmt, 1, pipe_name.c_str(), pipe_name.length(), nullptr);
        rc = sqlite3_bind_int(stmt, 2, from);
        rc = sqlite3_bind_int(stmt, 3, to);
        rc = sqlite3_bind_double(stmt, 4, settings[i].diameter);
        rc = sqlite3_bind_double(stmt, 5, settings[i].length);
        rc = sqlite3_bind_double(stmt, 6, settings[i].roughness);
        rc = sqlite3_bind_int(stmt, 7, settings[i].ref_nsegs);
        rc = sqlite3_step(stmt);
        rc = sqlite3_clear_bindings(stmt);
        rc = sqlite3_reset(stmt);
    }
    
    rc = sqlite3_exec(db, "COMMIT", nullptr, nullptr, nullptr);

    rc = sqlite3_finalize(stmt);

    return 0; 
}

} // namespace database

} // namespace shimmer
