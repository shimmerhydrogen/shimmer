#include <optional>
#include <cassert>
#include <sqlite3.h>
#include "errors.h"
#include "sqlite.hpp"

/* Database I/O functions for station type 'consumption point w/o pressure control' */

namespace shimmer {

namespace outlet_priv {

/* Define indices of columns in profile table */
enum class profile_col : int {
    s_number = 0,
    prf_time = 1,
    prf_Lset = 2
};

} // namespace outlet_priv

int
network_database::import_outlet(std::vector<setting_outlet>& settings)
{
    using namespace outlet_priv;

    auto tabnames_opt = limits_and_profile_table_names(station_type::PRIVATE_OUTLET);
    if ( not tabnames_opt ) {
        std::cerr << "Shimmer DB: cannot retrieve table names for 'outlet' station" << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    }

    auto [limits_tab, profile_tab] = tabnames_opt.value();

    sqlite3_stmt *stmt = nullptr;

    std::string qprfs = "SELECT s_number FROM " + profile_tab + " GROUP BY s_number";
    int rc = sqlite3_prepare_v2(db_, qprfs.c_str(), qprfs.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << qprfs << "': ";
        std::cerr << sqlite3_errmsg(db_) << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    }

    while (sqlite3_step(stmt) == SQLITE_ROW) {
        setting_outlet setting;
        setting.u_snum = sqlite3_column_int(stmt, 0);
        auto i_snum_opt = s_u2i.at(setting.u_snum);
        if (not i_snum_opt) {
            std::cerr << "s_u2i: invalid station number. Inconsistent data in DB?" << std::endl;
            return SHIMMER_DATABASE_PROBLEM;
        }

        setting.i_snum = i_snum_opt.value();
        settings.push_back( std::move(setting) );
    }

    std::string qprof = "SELECT * FROM " + profile_tab + " WHERE s_number = ?";
    rc = sqlite3_prepare_v2(db_, qprof.c_str(), qprof.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << qprof << "': ";
        std::cerr << sqlite3_errmsg(db_) << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    }

    /* Import profiles for all the stations */
    for (auto& setting : settings) {
        rc = sqlite3_bind_int(stmt, 1, setting.u_snum);
        std::vector<sample> profile;
        while (sqlite3_step(stmt) == SQLITE_ROW) {
            sample s;
            s.time  = sqlite3_column_double(stmt, +profile_col::prf_time);
            s.value = sqlite3_column_double(stmt, +profile_col::prf_Lset);
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
    return SHIMMER_SUCCESS;
}

} // namespace shimmer