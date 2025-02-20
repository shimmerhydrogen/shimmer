#include <optional>
#include <cassert>
#include <sqlite3.h>

#include "sqlite.hpp"

/* Database I/O functions for station type 'consumption point w/o pressure control' */

namespace shimmer {

namespace exit_l_reg_priv {

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

} // namespace exit_l_reg_priv

int
network_database::import_exit_l_reg(std::vector<setting_exit_l_reg>& settings)
{
    using namespace exit_l_reg_priv;

    auto tabnames_opt = limits_and_profile_table_names(0);
    if ( not tabnames_opt ) {
        return -1;
    }

    auto [limits_tab, profile_tab] = tabnames_opt.value();

    char *zErrMsg = nullptr;
    sqlite3_stmt *stmt = nullptr;

    std::string qlim = "SELECT * FROM " + limits_tab;
    int rc = sqlite3_prepare_v2(db_, qlim.c_str(), qlim.length(), &stmt, nullptr);
    if (rc) {
        fprintf(stderr, "%s:%d SQL error: %s\n", __FILE__, __LINE__, zErrMsg);
        sqlite3_free(zErrMsg);
        return 1;
    }

    settings.clear();

    /* Import limits for all the stations */
    while (sqlite3_step(stmt) == SQLITE_ROW) {
        setting_exit_l_reg setting;
        setting.u_snum = sqlite3_column_int(stmt, +limits_col::s_number);
        
        auto i_snum_opt = s_u2i.at(setting.u_snum);
        assert(i_snum_opt && "s_u2i: invalid station number. Inconsistent data in DB?");
        setting.i_snum = i_snum_opt.value();

        setting.Lmin = sqlite3_column_double(stmt, +limits_col::lim_Lmin);
        setting.Lmax = sqlite3_column_double(stmt, +limits_col::lim_Lmax);
        setting.Pmin = sqlite3_column_double(stmt, +limits_col::lim_Pmin);
        setting.Pmax = sqlite3_column_double(stmt, +limits_col::lim_Pmax);
        settings.push_back( std::move(setting) );
    }
    rc = sqlite3_finalize(stmt);

    std::string qprof = "SELECT * FROM " + profile_tab + " WHERE s_number = ?";
    rc = sqlite3_prepare_v2(db_, qlim.c_str(), qlim.length(), &stmt, nullptr);
    if (rc) {
        fprintf(stderr, "%s:%d SQL error: %s\n", __FILE__, __LINE__, zErrMsg);
        sqlite3_free(zErrMsg);
        return 1;
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
            std::cout << "no mass flow profile data defined." << std::endl;
        }

        rc = sqlite3_clear_bindings(stmt);
        rc = sqlite3_reset(stmt);
        std::sort(profile.begin(), profile.end());
        setting.Lprofile = std::move(profile);
    }
    
    sqlite3_finalize(stmt);
    std::sort(settings.begin(), settings.end());
    return 0;
}

} // namespace shimmer