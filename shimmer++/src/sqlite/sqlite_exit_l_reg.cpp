#include <optional>
#include <cassert>
#include <sqlite3.h>
#include "errors.h"
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
    prf_Lset = 2
};

} // namespace exit_l_reg_priv

namespace database {

static int
load_limits(sqlite3 *db, const optvector<int>& s_u2i,
    std::vector<setting_exit_l_reg>& settings)
{
    using namespace exit_l_reg_priv;
    auto tname_opt = table_name(db, setting_table::limits,
        station_type::EXIT_L_REG);
    if ( not tname_opt ) {
        std::cerr << "Shimmer DB: cannot retrieve limits table name ";
        std::cerr << "for 'exit_l_reg' station" << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    }

    sqlite3_stmt *stmt = nullptr;

    std::string qlim = "SELECT * FROM " + tname_opt.value();
    int rc = sqlite3_prepare_v2(db, qlim.c_str(), qlim.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << qlim << "': ";
        std::cerr << sqlite3_errmsg(db) << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    }

    /* Import limits for all the stations */
    while (sqlite3_step(stmt) == SQLITE_ROW) {
        setting_exit_l_reg setting;
        setting.u_snum = sqlite3_column_int(stmt, +limits_col::s_number);
        
        auto i_snum_opt = convert_u2i(s_u2i, setting.u_snum);
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

    return SHIMMER_SUCCESS;
}

static int
load_profiles(sqlite3 *db, const optvector<int>& s_u2i,
    std::vector<setting_exit_l_reg>& settings)
{
    using namespace exit_l_reg_priv;
    auto tname_opt = table_name(db, setting_table::profiles,
        station_type::EXIT_L_REG);
    if ( not tname_opt ) {
        std::cerr << "Shimmer DB: cannot retrieve profile table name ";
        std::cerr << "for 'exit_l_reg' station" << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    }

    sqlite3_stmt *stmt = nullptr;

    std::string qprof = "SELECT * FROM " + tname_opt.value() + " WHERE s_number = ?";
    int rc = sqlite3_prepare_v2(db, qprof.c_str(), qprof.length(), &stmt, nullptr);
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
            std::cout << "no pressure profile data defined." << std::endl;
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

int load(sqlite3 *db, const optvector<int>& s_u2i,
    std::vector<setting_exit_l_reg>& settings)
{
    load_limits(db, s_u2i, settings);
    load_profiles(db, s_u2i, settings);
    return SHIMMER_SUCCESS;
}

static int store_limits(sqlite3 *db, const std::vector<int>& s_i2u,
    const std::vector<setting_exit_l_reg>& settings)
{
    using namespace exit_l_reg_priv;
    auto tname_opt = table_name(db, setting_table::limits,
        station_type::EXIT_L_REG);
    if ( not tname_opt ) {
        std::cerr << "Shimmer DB: cannot retrieve limits table name ";
        std::cerr << "for 'exit_l_reg' station" << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    }

    sqlite3_stmt *stmt = nullptr;

    std::string q =
        "INSERT INTO " + tname_opt.value() + " VALUES (?, ?, ?, ?, ?)";
        
    int rc = sqlite3_prepare_v2(db, q.c_str(), q.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << q << "': ";
        std::cerr << sqlite3_errmsg(db) << std::endl;
        return SHIMMER_DATABASE_PROBLEM;
    }

    rc = sqlite3_exec(db, "BEGIN TRANSACTION", nullptr, nullptr, nullptr);

    for (auto& setting : settings) {
        rc = sqlite3_bind_int(stmt, 1, convert_i2u(s_i2u, setting.i_snum) );
        rc = sqlite3_bind_double(stmt, 2, setting.Lmin);
        rc = sqlite3_bind_double(stmt, 3, setting.Lmax);
        rc = sqlite3_bind_double(stmt, 4, setting.Pmin);
        rc = sqlite3_bind_double(stmt, 5, setting.Pmax);
        rc = sqlite3_step(stmt);
        rc = sqlite3_clear_bindings(stmt);
        rc = sqlite3_reset(stmt);
    }

    rc = sqlite3_exec(db, "COMMIT", nullptr, nullptr, nullptr);

    return SHIMMER_SUCCESS;
}

static int store_profiles(sqlite3 *db, const std::vector<int>& s_i2u,
    const std::vector<setting_exit_l_reg>& settings)
{
    using namespace exit_l_reg_priv;
    auto tname_opt = table_name(db, setting_table::profiles,
        station_type::EXIT_L_REG);
    if ( not tname_opt ) {
        std::cerr << "Shimmer DB: cannot retrieve profile table name ";
        std::cerr << "for 'exit_l_reg' station" << std::endl;
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
        return SHIMMER_DATABASE_PROBLEM;
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

int store(sqlite3 *db, const std::vector<int>& s_i2u,
    const std::vector<setting_exit_l_reg>& settings)
{
    store_limits(db, s_i2u, settings);
    store_profiles(db, s_i2u, settings);
    return SHIMMER_SUCCESS;
}

} // namespace database

} // namespace shimmer