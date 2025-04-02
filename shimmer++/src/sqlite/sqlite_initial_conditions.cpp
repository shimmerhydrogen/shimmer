#include <optional>
#include <sqlite3.h>
#include "errors.h"
#include "sqlite.hpp"

namespace shimmer {

int
network_database::import_station_initial_conditions()
{
    int rc;
    char *zErrMsg = nullptr;
    std::string zSql = "SELECT * FROM station_initial_conditions";

    sqlite3_stmt *stmt;
    rc = sqlite3_prepare_v2(db_, zSql.c_str(), zSql.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << zSql << "': " << zErrMsg << std::endl;
        sqlite3_free(zErrMsg);
        return SHIMMER_DATABASE_PROBLEM;
    }

    while (sqlite3_step(stmt) != SQLITE_DONE) {
        int u_snum = sqlite3_column_int(stmt, 0);
        int i_snum = s_u2vd.at(u_snum).value();
        station_initial_condition sic;
        sic.i_snum = i_snum;
        sic.init_P = sqlite3_column_double(stmt, 1);
        sic.init_L = sqlite3_column_double(stmt, 2);
        sics.push_back(sic);
    }

    sqlite3_finalize(stmt);
    return SHIMMER_SUCCESS;
}

int
network_database::import_pipe_initial_conditions()
{
    int rc;
    char *zErrMsg = nullptr;
    std::string zSql = "SELECT * FROM pipe_initial_conditions";

    sqlite3_stmt *stmt;
    rc = sqlite3_prepare_v2(db_, zSql.c_str(), zSql.length(), &stmt, nullptr);
    if (rc) {
        std::cerr << "SQL error on query '" << zSql << "': " << zErrMsg << std::endl;
        sqlite3_free(zErrMsg);
        return SHIMMER_DATABASE_PROBLEM;
    }

    while (sqlite3_step(stmt) != SQLITE_DONE) {
        std::string name = (char *) sqlite3_column_text(stmt, 0);
        int u_sfrom = sqlite3_column_int(stmt, 1);
        int i_sfrom = s_u2vd.at(u_sfrom).value();
        int u_sto = sqlite3_column_int(stmt, 2);
        int i_sto = s_u2vd.at(u_sto).value();

        pipe_initial_condition pic;
        pic.i_sfrom = i_sfrom;
        pic.i_sto = i_sto;
        pic.init_G = sqlite3_column_double(stmt, 3);
        pics.push_back(pic);
    }

    sqlite3_finalize(stmt);
    return SHIMMER_SUCCESS;
}

} //namespace shimmer