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

/* Database I/O functions for pipe type 'Compressor' */

namespace shimmer {

namespace compr_stat {

/* Define indices of columns in profile table */
enum class limits_col : int {
    s_from          = 1,
    s_to            = 2,
    max_power       = 3,
    max_outpress    = 4,
    min_inpress     = 5,
    max_ratio       = 6,
    min_ratio       = 7,
    max_massflow    = 8,
};

/* Define indices of columns in profile table */
enum class profile_col : int {
    s_from      = 1,
    s_to        = 2,
    prf_time    = 3,
    controlmode = 4,
    power       = 5,
    outpress    = 6,
    inpress     = 7,
    ratio       = 8,
    massflow    = 9,
};

} // namespace compr_stat

namespace database {

int
load(sqlite3 *db, const optvector<int>& s_u2i,
    std::vector<setting_compr_stat>& settings)
{
    using namespace compr_stat;

    sqlite3_stmt *stmt = nullptr;

    std::string qlim = "SELECT compressor_limits.* "
	    "FROM compressor_limits INNER JOIN pipelines "
	    "ON compressor_limits.p_name = pipelines.p_name "
		"   AND compressor_limits.s_from = pipelines.s_from "
		"   AND compressor_limits.s_to = pipelines.s_to "
	    "WHERE pipelines.p_type = ?";
    int rc = sqlite3_prepare_v2(db, qlim.c_str(), qlim.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << qlim << "': ";
        std::cerr << sqlite3_errmsg(db) << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    }

    rc = sqlite3_bind_int(stmt, 1, +pipe_type::COMPR_STAT);

    /* Import limits for all the stations */
    while (sqlite3_step(stmt) == SQLITE_ROW) {
        setting_compr_stat setting;
        int u_from = sqlite3_column_int(stmt, +limits_col::s_from);
        int u_to = sqlite3_column_int(stmt, +limits_col::s_to);
        
        auto i_sfrom_opt = s_u2i.at(u_from);
        auto i_sto_opt = s_u2i.at(u_to);
        if (not i_sfrom_opt or not i_sto_opt) {
            std::cerr << "s_u2i: invalid station numbers. Inconsistent data in DB?" << std::endl;
            return SHIMMER_DATABASE_PROBLEM;
        }
        setting.u_sfrom = u_from;
        setting.u_sto = u_to;
        setting.i_sfrom = i_sfrom_opt.value();
        setting.i_sto = i_sto_opt.value();
        setting.max_power = sqlite3_column_double(stmt, +limits_col::max_power);
        setting.max_outpress = sqlite3_column_double(stmt, +limits_col::max_outpress);
        setting.min_inpress = sqlite3_column_double(stmt, +limits_col::min_inpress);
        setting.max_ratio = sqlite3_column_double(stmt, +limits_col::max_ratio);
        setting.min_ratio = sqlite3_column_double(stmt, +limits_col::min_ratio);
        setting.max_massflow = sqlite3_column_double(stmt, +limits_col::max_massflow);
        setting.efficiency = 1.0;
        settings.push_back( std::move(setting) );
    }
    rc = sqlite3_finalize(stmt);

    std::string qprof = "SELECT * FROM compressor_profile WHERE s_from = ? AND s_to = ?";
    rc = sqlite3_prepare_v2(db, qprof.c_str(), qprof.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << qlim << "': ";
        std::cerr << sqlite3_errmsg(db) << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    }

    /* Import profiles for all the stations */
    int missing_settings = 0;
    for (auto& setting : settings) {
        rc = sqlite3_bind_int(stmt, 1, setting.u_sfrom);
        rc = sqlite3_bind_int(stmt, 2, setting.u_sto);
        std::vector<compressor_profile_sample> profile;
        while (sqlite3_step(stmt) == SQLITE_ROW) {
            compressor_profile_sample s;
            s.mode = static_cast<compressor_mode>(
                sqlite3_column_int(stmt, +profile_col::controlmode)
            );
            s.time  = sqlite3_column_double(stmt, +profile_col::prf_time);
            s.power =  sqlite3_column_double(stmt, +profile_col::power);
            s.outpress = sqlite3_column_double(stmt, +profile_col::outpress);
            s.inpress = sqlite3_column_double(stmt, +profile_col::inpress);
            s.ratio = sqlite3_column_double(stmt, +profile_col::ratio);
            s.massflow = sqlite3_column_double(stmt, +profile_col::massflow);
            profile.push_back(s);
        }

        if (profile.size() == 0) {
            std::cout << "Warning: compressor (" << setting.u_sfrom << ", ";
            std::cerr << setting.u_sto << ") has ";
            std::cout << "no pressure profile data defined." << std::endl;
            missing_settings++;
        }

        rc = sqlite3_clear_bindings(stmt);
        rc = sqlite3_reset(stmt);
        std::sort(profile.begin(), profile.end());
        setting.profile = std::move(profile);
    }
    
    sqlite3_finalize(stmt);
    std::sort(settings.begin(), settings.end());

    if (missing_settings > 0)
        return SHIMMER_MISSING_DATA;

    return SHIMMER_SUCCESS;
}

int
store(sqlite3 *db, const std::vector<int>& s_i2u,
    std::vector<setting_compr_stat>& settings)
{
    return SHIMMER_MISSING_DATA;
}

} //namespace database

} // namespace shimmer