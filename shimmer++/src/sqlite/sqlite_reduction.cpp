#include <optional>
#include <sqlite3.h>
#include "errors.h"
#include "sqlite.hpp"

/* Database I/O functions for pipe type 'Compressor' */

namespace shimmer {

namespace red_stat {

/* Define indices of columns in profile table */
enum class limits_col : int {
    s_from          = 1,
    s_to            = 2,
    max_outpress    = 3,
    min_inpress     = 4,
    max_ratio       = 5,
    min_ratio       = 6,
    max_massflow    = 7,
};

/* Define indices of columns in profile table */
enum class profile_col : int {
    s_from      = 1,
    s_to        = 2,
    prf_time    = 3,
    controlmode = 4,
    outpress    = 5,
    inpress     = 6,
    ratio       = 7,
    massflow    = 8,
};

} // namespace red_stat

int
network_database::import_red_stat(std::vector<setting_red_stat>& settings)
{
    using namespace red_stat;

    char *zErrMsg = nullptr;
    sqlite3_stmt *stmt = nullptr;

    std::string qlim = "SELECT * FROM reduction_limits";
    int rc = sqlite3_prepare_v2(db_, qlim.c_str(), qlim.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << qlim << "': " << zErrMsg << std::endl;
        sqlite3_free(zErrMsg);
        return SHIMMER_DATABASE_PROBLEM;
    }

    /* Import limits for all the stations */
    while (sqlite3_step(stmt) == SQLITE_ROW) {
        setting_red_stat setting;
        int u_from = sqlite3_column_int(stmt, +limits_col::s_from);
        int u_to = sqlite3_column_int(stmt, +limits_col::s_to);
        
        auto i_sfrom_opt = nd_.s_u2i.at(u_from);
        auto i_sto_opt = nd_.s_u2i.at(u_to);
        if (not i_sfrom_opt or not i_sto_opt) {
            std::cerr << "s_u2i: invalid station numbers. Inconsistent data in DB?" << std::endl;
            return SHIMMER_DATABASE_PROBLEM;
        }

        setting.i_sfrom = i_sfrom_opt.value();
        setting.i_sto = i_sto_opt.value();
        setting.max_outpress = sqlite3_column_double(stmt, +limits_col::max_outpress);
        setting.min_inpress = sqlite3_column_double(stmt, +limits_col::min_inpress);
        setting.max_ratio = sqlite3_column_double(stmt, +limits_col::max_ratio);
        setting.min_ratio = sqlite3_column_double(stmt, +limits_col::min_ratio);
        setting.max_massflow = sqlite3_column_double(stmt, +limits_col::max_massflow);
        settings.push_back( std::move(setting) );
    }
    rc = sqlite3_finalize(stmt);

    std::string qprof = "SELECT * FROM reduction_profile WHERE s_from = ? AND s_to = ?";
    rc = sqlite3_prepare_v2(db_, qlim.c_str(), qlim.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << qprof << "': " << zErrMsg << std::endl;
        sqlite3_free(zErrMsg);
        return SHIMMER_DATABASE_PROBLEM;
    }

    /* Import profiles for all the stations */
    for (auto& setting : settings) {
        rc = sqlite3_bind_int(stmt, 1, setting.i_sfrom);
        rc = sqlite3_bind_int(stmt, 2, setting.i_sto);
        std::vector<reduction_profile_sample> profile;
        while (sqlite3_step(stmt) == SQLITE_ROW) {
            reduction_profile_sample s;
            s.mode = static_cast<reduction_mode>(
                sqlite3_column_int(stmt, +profile_col::controlmode)
            );
            s.time  = sqlite3_column_double(stmt, +profile_col::prf_time);
            s.outpress = sqlite3_column_double(stmt, +profile_col::outpress);
            s.inpress = sqlite3_column_double(stmt, +profile_col::inpress);
            s.ratio = sqlite3_column_double(stmt, +profile_col::ratio);
            s.massflow = sqlite3_column_double(stmt, +profile_col::massflow);
            profile.push_back(s);
        }

        if (profile.size() == 0) {
            std::cout << "Warning: pipe from " << setting.i_sfrom << " to ";
            std::cerr << setting.i_sto << " has ";
            std::cout << "no pressure profile data defined." << std::endl;
        }

        rc = sqlite3_clear_bindings(stmt);
        rc = sqlite3_reset(stmt);
        std::sort(profile.begin(), profile.end());
        setting.profile = std::move(profile);
    }
    
    sqlite3_finalize(stmt);
    std::sort(settings.begin(), settings.end());
    return SHIMMER_SUCCESS;
}

} // namespace shimmer