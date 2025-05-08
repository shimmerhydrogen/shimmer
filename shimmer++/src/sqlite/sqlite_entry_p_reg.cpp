#include <optional>
#include <cassert>
#include <sqlite3.h>
#include "errors.h"
#include "sqlite.hpp"

/* Database I/O functions for station type 'ReMi station w/o backflow' */

namespace shimmer {
namespace entry_p_reg_priv {

/* Define indices of columns in limits table */
enum class limits_col : int {
    s_number = 0,
    lim_Lmin = 1,
    lim_Lmax = 2,
    lim_Pmin = 3,
    lim_Pmax = 4,
};

/* Define indices of columns in profile table */
enum class profile_col : int {
    s_number = 0,
    prf_time = 1,
    prf_Pset = 2
};

} // namespace entry_p_reg_priv

int
network_database::import_entry_p_reg(std::vector<setting_entry_p_reg>& settings)
{
    using namespace entry_p_reg_priv;

    auto tabnames_opt = limits_and_profile_table_names(station_type::ENTRY_P_REG);
    if ( not tabnames_opt ) {
        std::cerr << "Shimmer DB: cannot retrieve table names for 'entry_l_reg' station" << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    }

    auto [limits_tab, profile_tab] = tabnames_opt.value();

    char *zErrMsg = nullptr;
    sqlite3_stmt *stmt = nullptr;

    std::string qlim = "SELECT * FROM " + limits_tab;
    int rc = sqlite3_prepare_v2(db_, qlim.c_str(), qlim.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << qlim << "': " << zErrMsg << std::endl;
        sqlite3_free(zErrMsg);
        return SHIMMER_DATABASE_PROBLEM;
    }

    settings.clear();

    /* Import limits for all the stations */
    while (sqlite3_step(stmt) == SQLITE_ROW) {
        setting_entry_p_reg setting;
        setting.u_snum = sqlite3_column_int(stmt, +limits_col::s_number);
        
        auto i_snum_opt = s_u2i.at(setting.u_snum);
        if (not i_snum_opt) {
            std::cerr << "s_u2i: invalid station number. Inconsistent data in DB?" << std::endl;
            return SHIMMER_DATABASE_PROBLEM;
        }

        setting.i_snum = i_snum_opt.value();
        setting.Lmin = sqlite3_column_double(stmt, +limits_col::lim_Lmin);
        setting.Lmax = sqlite3_column_double(stmt, +limits_col::lim_Lmax);
        setting.Pmin = sqlite3_column_double(stmt, +limits_col::lim_Pmin);
        setting.Pmax = sqlite3_column_double(stmt, +limits_col::lim_Pmax);
        settings.push_back( std::move(setting) );
    }
    rc = sqlite3_finalize(stmt);

    std::string qprof = "SELECT * FROM " + profile_tab + " WHERE s_number = ?";
    rc = sqlite3_prepare_v2(db_, qprof.c_str(), qprof.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << qprof << "': " << zErrMsg << std::endl;
        sqlite3_free(zErrMsg);
        return SHIMMER_DATABASE_PROBLEM;
    }

    /* Import profiles for all the stations */
    for (auto& setting : settings) {
        rc = sqlite3_bind_int(stmt, 1, setting.u_snum);
        std::vector<sample> profile;
        while (sqlite3_step(stmt) == SQLITE_ROW) {
            sample s;
            s.time  = sqlite3_column_double(stmt, +profile_col::prf_time);
            s.value = sqlite3_column_double(stmt, +profile_col::prf_Pset);
            profile.push_back(s);
        }

        if (profile.size() == 0) {
            std::cout << "Warning: node " << setting.u_snum << " has ";
            std::cout << "no pressure profile data defined." << std::endl;
        }

        rc = sqlite3_clear_bindings(stmt);
        rc = sqlite3_reset(stmt);
        std::sort(profile.begin(), profile.end());
        setting.Pprofile = std::move(profile);
    }
    
    sqlite3_finalize(stmt);
    std::sort(settings.begin(), settings.end());
    return SHIMMER_SUCCESS;
}

} // namespace shimmer