#include <optional>
#include <sqlite3.h>
#include "errors.h"
#include "sqlite.hpp"

/* Database I/O functions for pipe type 'Compressor' */

namespace shimmer {

namespace pipe {

/* Define indices of columns in parameters table */
enum class params_col : int {
    s_from      = 1,
    s_to        = 2,
    diameter    = 3,
    length      = 4,
    roughness   = 5,
};

} // namespace pipe

int
network_database::import_pipe(std::vector<setting_pipe>& settings)
{
    using namespace pipe;

    char *zErrMsg = nullptr;
    sqlite3_stmt *stmt = nullptr;

    std::string qlim = "SELECT * FROM pipe_parameters";
    int rc = sqlite3_prepare_v2(db_, qlim.c_str(), qlim.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << qlim << "': " << zErrMsg << std::endl;
        sqlite3_free(zErrMsg);
        return SHIMMER_DATABASE_PROBLEM;
    }

    /* Import limits for all the stations */
    while (sqlite3_step(stmt) == SQLITE_ROW) {
        setting_pipe setting;
        int u_from = sqlite3_column_int(stmt, +params_col::s_from);
        int u_to = sqlite3_column_int(stmt, +params_col::s_to);
        
        auto i_sfrom_opt = s_u2i.at(u_from);
        auto i_sto_opt = s_u2i.at(u_to);
        if (not i_sfrom_opt or not i_sto_opt) {
            std::cerr << "s_u2i: invalid station numbers. Inconsistent data in DB?" << std::endl;
            return SHIMMER_DATABASE_PROBLEM;
        }

        setting.i_sfrom = i_sfrom_opt.value();
        setting.i_sto = i_sto_opt.value();
        setting.diameter = sqlite3_column_double(stmt, +params_col::diameter);
        setting.length = sqlite3_column_double(stmt, +params_col::length);
        setting.roughness = sqlite3_column_double(stmt, +params_col::roughness);
        settings.push_back( std::move(setting) );
    }
    rc = sqlite3_finalize(stmt);

    std::sort(settings.begin(), settings.end());
    return SHIMMER_SUCCESS;
}

} // namespace shimmer