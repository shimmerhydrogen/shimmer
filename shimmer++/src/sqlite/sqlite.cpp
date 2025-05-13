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

#include <iostream>
#include <cstdio>
#include "errors.h"
#include "sqlite.hpp"

namespace shimmer {

int check_sqlite(sqlite3 *db, int rc)
{
#ifdef DEBUG_SQL
    if (rc) {
        std::cerr << "SQL problem: " << sqlite3_errmsg(db) << std::endl;
        throw std::invalid_argument("");
    }
#endif
    return rc;
}

std::optional<table_name_pair_t>
limits_and_profile_table_names(sqlite3 *db, station_type stat_type)
{
    sqlite3_stmt *stmt = nullptr;

    enum class col : int {
        t_limits_table = 0,
        t_profile_table = 1
    };

    std::string q =
        "SELECT t_limits_table, t_profile_table "
        "FROM station_types "
        "WHERE t_type = ?";

    int rc = sqlite3_prepare_v2(db, q.c_str(), q.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << q << "': " << sqlite3_errmsg(db) << std::endl;
        return {};
    }

    rc = sqlite3_bind_int(stmt, 1, +stat_type);

    if ( sqlite3_step(stmt) != SQLITE_ROW ) {
        std::cerr << "Shimmer DB: Invalid station type " << +stat_type << std::endl;
        return {};
    }
    
    std::string limits_table_name;
    if ( sqlite3_column_type(stmt, +col::t_limits_table) == SQLITE3_TEXT ) {
        limits_table_name = (const char *) sqlite3_column_text(stmt, +col::t_limits_table);
    }

    std::string profile_table_name;
    if ( sqlite3_column_type(stmt, +col::t_profile_table) == SQLITE3_TEXT ) {
        profile_table_name = (const char *) sqlite3_column_text(stmt, +col::t_profile_table);
    }

    if ( (limits_table_name == "") or (profile_table_name == "") ) {
        return {};
    }

    rc = sqlite3_clear_bindings(stmt);
    rc = sqlite3_reset(stmt);
    rc = sqlite3_finalize(stmt);

    return std::pair(limits_table_name, profile_table_name);
}

std::optional<std::string>
table_name(sqlite3 *db, setting_table st, station_type stat_type)
{
    sqlite3_stmt *stmt = nullptr;

    std::string q;
    if ( st == setting_table::limits ) {
        q = "SELECT t_limits_table "
        "FROM station_types "
        "WHERE t_type = ?";
    } else if ( st == setting_table::profiles ) {
        q = "SELECT t_profile_table "
        "FROM station_types "
        "WHERE t_type = ?";
    } else {
        std::cerr << "Invalid setting table specified" << std::endl;
        return {};
    }

    int rc = sqlite3_prepare_v2(db, q.c_str(), q.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << q << "': ";
        std::cerr << sqlite3_errmsg(db) << std::endl;
        return {};
    }

    rc = sqlite3_bind_int(stmt, 1, +stat_type);

    if ( sqlite3_step(stmt) != SQLITE_ROW ) {
        std::cerr << "Shimmer DB: Invalid station type " << +stat_type << std::endl;
        return {};
    }
    
    std::string tname;
    if ( sqlite3_column_type(stmt, /*column*/ 0) == SQLITE3_TEXT ) {
        tname = (const char *) sqlite3_column_text(stmt, /*column*/ 0);
    }

    if (tname == "") {
        return {};
    }

    rc = sqlite3_clear_bindings(stmt);
    rc = sqlite3_reset(stmt);
    rc = sqlite3_finalize(stmt);

    return tname;
}



} // namespace shimmer


