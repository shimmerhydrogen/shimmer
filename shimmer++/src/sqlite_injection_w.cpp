#include <optional>
#include <cassert>
#include <sqlite3.h>

#include "sqlite.hpp"

/* Database I/O functions for station type 'ReMi station w/o backflow' */

namespace shimmer {

int
network_database::import_injection_w(std::vector<setting_injection_w>& settings)
{
    char *zErrMsg = nullptr;
    sqlite3_stmt *stmt = nullptr;

    std::string qlim = "SELECT * FROM limits_injection_w";
    int rc = sqlite3_prepare_v2(db_, qlim.c_str(), qlim.length(), &stmt, nullptr);
    if (rc) {
        fprintf(stderr, "%s:%d SQL error: %s\n", __FILE__, __LINE__, zErrMsg);
        sqlite3_free(zErrMsg);
        return 1;
    }

    /* Import limits for all the stations */
    while (sqlite3_step(stmt) == SQLITE_ROW) {
        setting_injection_w setting;
        setting.u_snum = sqlite3_column_int(stmt, 0);
        
        auto i_snum_opt = s_u2i.at(setting.u_snum);
        assert(i_snum_opt && "s_u2i: invalid station number. Inconsistent data in DB?");
        setting.i_snum = i_snum_opt.value();

        setting.Lmin = sqlite3_column_double(stmt, 1);
        setting.Lmax = sqlite3_column_double(stmt, 2);
        setting.Pmin = sqlite3_column_double(stmt, 3);
        setting.Pmax = sqlite3_column_double(stmt, 4);
        setting.f    = sqlite3_column_double(stmt, 5);
        settings.push_back( std::move(setting) );
    }
    rc = sqlite3_finalize(stmt);

    std::string qprof = "SELECT * FROM profiles_injection_w WHERE s_number = ?";
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
            s.time  = sqlite3_column_double(stmt, 1);
            s.value = sqlite3_column_double(stmt, 2);
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
    return 0;
}

} // namespace shimmer