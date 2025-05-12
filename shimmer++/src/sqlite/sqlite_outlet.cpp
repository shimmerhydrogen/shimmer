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

namespace database {
int load(sqlite3 *db, const optvector<int>& s_u2i,
    std::vector<setting_outlet>& settings)
{
    using namespace outlet_priv;

    auto tname_opt = table_name(db, setting_table::profiles,
        station_type::PRIVATE_OUTLET);
    if ( not tname_opt ) {
        std::cerr << "Shimmer DB: cannot retrieve profile table name ";
        std::cerr << "for 'outlet' station" << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    }

    sqlite3_stmt *stmt = nullptr;

    std::string q = "SELECT s_number FROM " + tname_opt.value() +
        " GROUP BY s_number";
    int rc = sqlite3_prepare_v2(db, q.c_str(), q.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << q << "': ";
        std::cerr << sqlite3_errmsg(db) << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    }

    while (sqlite3_step(stmt) == SQLITE_ROW) {
        setting_outlet setting;
        setting.u_snum = sqlite3_column_int(stmt, /* column */ 0);
        auto i_snum_opt = convert_u2i(s_u2i, setting.u_snum);
        if (not i_snum_opt) {
            std::cerr << "s_u2i: invalid station number. Inconsistent data in DB?" << std::endl;
            return SHIMMER_DATABASE_PROBLEM;
        }

        setting.i_snum = i_snum_opt.value();
        settings.push_back( std::move(setting) );
    }

    std::string qprof = "SELECT * FROM " + tname_opt.value() +
        " WHERE s_number = ?";
    rc = sqlite3_prepare_v2(db, qprof.c_str(), qprof.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << qprof << "': ";
        std::cerr << sqlite3_errmsg(db) << std::endl;
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

int store(sqlite3 *db, const std::vector<int>& s_i2u,
    const std::vector<setting_outlet>& settings)
{
    using namespace outlet_priv;

    auto tname_opt = table_name(db, setting_table::profiles,
        station_type::PRIVATE_OUTLET);
    if ( not tname_opt ) {
        std::cerr << "Shimmer DB: cannot retrieve profile table name ";
        std::cerr << "for 'outlet' station" << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    }

    sqlite3_stmt *stmt = nullptr;

    std::string q =
        "INSERT INTO " + tname_opt.value() + " (s_number, prf_time, prf_Lset) "
        "VALUES (?, ?, ?)";
        
    int rc = sqlite3_prepare_v2(db, q.c_str(), q.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << q << "': ";
        std::cerr << sqlite3_errmsg(db) << std::endl;
        return -1;
    }

    rc = sqlite3_exec(db, "BEGIN TRANSACTION", nullptr, nullptr, nullptr);

    for (auto& setting : settings) {
        for (auto& sample : setting.Lprofile) {
            rc = sqlite3_bind_int(stmt, 1, convert_i2u(s_i2u, setting.i_snum));
            rc = sqlite3_bind_int(stmt, 2, sample.time);
            rc = sqlite3_bind_int(stmt, 3, sample.value);
            rc = sqlite3_step(stmt);
            rc = sqlite3_clear_bindings(stmt);
            rc = sqlite3_reset(stmt);
        }
    }
    
    rc = sqlite3_exec(db, "COMMIT", nullptr, nullptr, nullptr);

    rc = sqlite3_finalize(stmt);

    return SHIMMER_SUCCESS;
}

} // namespace database



} // namespace shimmer